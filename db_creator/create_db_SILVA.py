import glob, re, gzip
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import zipfile, gzip, io
from io import StringIO
import subprocess, shutil, os, datetime
from tqdm import tqdm
from Bio import Entrez
from ete3 import NCBITaxa
from concurrent.futures import ThreadPoolExecutor, as_completed
import sqlite3

# Provide email address
Entrez.email = "macher@uni-trier.de"

# Initialize ete3 NCBI database
ncbi = NCBITaxa()

def create_acc2taxid_db(gz_path, db_path):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("DROP TABLE IF EXISTS acc2taxid")
    c.execute("CREATE TABLE acc2taxid (accession TEXT PRIMARY KEY, taxid INTEGER)")

    with gzip.open(gz_path, "rt") as f:
        next(f)  # skip header
        batch = []
        for i, line in enumerate(f):
            fields = line.strip().split("\t")
            if len(fields) >= 3:
                accession = fields[0]
                taxid = int(fields[2])
                batch.append((accession, taxid))
            if i % 100000 == 0:
                c.executemany("INSERT OR IGNORE INTO acc2taxid VALUES (?, ?)", batch)
                conn.commit()
                batch = []
        if batch:
            c.executemany("INSERT OR IGNORE INTO acc2taxid VALUES (?, ?)", batch)
            conn.commit()

    conn.close()

def get_taxid_from_sqlite(accession, db_path):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("SELECT taxid FROM acc2taxid WHERE accession=?", (accession,))
    result = c.fetchone()
    conn.close()
    return result[0] if result else None

def get_taxonomy_from_accession_ete3(accession):
    try:
        # Step 1: Try local lookup
        taxid = get_taxid_from_sqlite(accession, db_path)

        # Step 2: If not found, fallback to Entrez
        if not taxid:
            handle = Entrez.esummary(db="nuccore", id=accession)
            summary = Entrez.read(handle)
            handle.close()
            taxid = int(summary[0]["TaxId"])

        # Step 3: Use ete3 to get lineage
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)

        # Map desired ranks
        wanted_ranks = {
            'superkingdom': 'unclassified',
            'phylum': 'unclassified',
            'class': 'unclassified',
            'order': 'unclassified',
            'family': 'unclassified',
            'genus': 'unclassified',
            'species': 'unclassified'
        }

        for tid in lineage:
            rank = ranks.get(tid)
            name = names.get(tid)
            if rank in wanted_ranks:
                wanted_ranks[rank] = name

        wanted_ranks['Accession'] = accession
        return wanted_ranks

    except Exception as e:
        return {'Accession': accession, 'Error': str(e)}

def zip_to_gz(fasta_file):
    # Determine the output .gz file path
    gz_path = os.path.splitext(fasta_file)[0] + '.gz'

    temp_dir = 'temp_unzip'
    os.makedirs(temp_dir, exist_ok=True)

    # Unzip the contents
    with zipfile.ZipFile(fasta_file, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    # Gzip the contents
    with open(gz_path, 'wb') as gz_file:
        for root, _, files in os.walk(temp_dir):
            for file in files:
                file_path = os.path.join(root, file)
                with open(file_path, 'rb') as f_in:
                    with gzip.open(gz_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

    # Clean up the temporary directory
    shutil.rmtree(temp_dir)
    os.remove(fasta_file)

    return gz_path

def silva_taxonomy(fasta_file, taxonomy_xlsx):
    # Read existing taxonomy
    try:
        all_accession_numbers_df = pd.read_excel(taxonomy_xlsx)
    except FileNotFoundError:
        all_accession_numbers_df = pd.DataFrame(columns=['Accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

    known_accessions = set(all_accession_numbers_df['Accession'].values)
    all_accession_numbers_df = all_accession_numbers_df.drop_duplicates(subset='Accession', keep='last')
    accession_dict = all_accession_numbers_df.set_index('Accession').T.to_dict('list')

    # Read sequences
    with gzip.open(fasta_file, 'rt') as handle:
        sequences = list(SeqIO.parse(handle, 'fasta'))

    n_seqs = len(sequences)

    results = []

    def process_record(record):
        id = record.id
        accession = id.split('.')[0]

        if id in known_accessions:
            return [id] + accession_dict[id]
        else:
            tax = get_taxonomy_from_accession_ete3(accession)
            if 'Error' in tax:
                assignment = record.description.split(';')[-1]
                return [id] + ['unclassified']*6 + [assignment]
            else:
                return [id] + [str(i) for i in list(tax.values())[:-1]]

    # Use threads or processes depending on I/O vs CPU time of get_taxonomy
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(process_record, record): record for record in sequences}
        for i, future in enumerate(as_completed(futures), 1):
            print(f"Processed {i}/{n_seqs}")
            results.append(future.result())

    # Build output DataFrame
    accession_df = pd.DataFrame(results, columns=['Accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

    # Combine with old and deduplicate
    out_df = pd.concat([all_accession_numbers_df, accession_df], ignore_index=True).drop_duplicates(subset=['Accession'])

    # Save
    out_df.to_excel(taxonomy_xlsx, index=False)

    return accession_df

def run_silva(output_path, fasta_file, taxonomy_xlsx):

    print('{} : Starting to collect taxonomy from fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## collect accession numbers
    accession_df = silva_taxonomy(fasta_file, taxonomy_xlsx)

    print('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    # create database
    db_name = Path(fasta_file).name.replace('.fasta.gz', '')
    db_folder = Path(output_path).joinpath(f'db_{db_name}')
    if not os.path.isdir(db_folder):
        os.mkdir(db_folder)
    db_folder = db_folder.joinpath('db')
    command = f'zcat < {fasta_file} | makeblastdb -in - -title db -dbtype nucl -out {db_folder}'
    os.system(command)

    # write taxonomy file
    taxonomy_file_snappy = db_folder.parent.joinpath('db_taxonomy.parquet.snappy')
    accession_df.to_parquet(taxonomy_file_snappy)

    ## zip the folder
    output = Path(output_path).joinpath(f'db_{db_name}')
    shutil.make_archive(output, 'zip', output)

    print('{} : Finished to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ####################################################################################################################

def run_silva_filtered(output_path, fasta_file, taxonomy_xlsx):

    print('{} : Starting to collect taxonomy from fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## collect accession numbers
    accession_df = silva_taxonomy(fasta_file, taxonomy_xlsx)

    print('{} : Starting to filter fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    accession_df_filtered = accession_df[accession_df['phylum'] != 'unclassified']
    accession_keep = set(accession_df_filtered['Accession'])
    fasta_file_filtered = Path(str(Path(fasta_file)).replace('.fasta.gz', '_sp.fasta.gz'))

    with gzip.open(fasta_file, "rt") as in_handle, gzip.open(fasta_file_filtered, "wt") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            if record.id in accession_keep:
                SeqIO.write(record, out_handle, "fasta")

    print('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    # create database
    db_name = Path(fasta_file).name.replace('.fasta.gz', '_sp')
    db_folder = Path(output_path).joinpath(f'db_{db_name}')
    if not os.path.isdir(db_folder):
        os.mkdir(db_folder)
    db_folder = db_folder.joinpath('db')
    command = f'zcat < {fasta_file_filtered} | makeblastdb -in - -title db -dbtype nucl -out {db_folder}'
    os.system(command)

    # write taxonomy file
    taxonomy_file_snappy = db_folder.parent.joinpath('db_taxonomy.parquet.snappy')
    accession_df.to_parquet(taxonomy_file_snappy)

    ## zip the folder
    output = Path(output_path).joinpath(f'db_{db_name}')
    shutil.make_archive(output, 'zip', output)

    print('{} : Finished to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ####################################################################################################################

## Variables
output_path = '/Volumes/Coruscant/APSCALE_raw_databases/2025_07'
files = glob.glob('/Volumes/Coruscant/APSCALE_raw_databases/2025_07_fasta/SILVA/*.fasta*')
taxonomy_xlsx = '/Volumes/Coruscant/APSCALE_raw_databases/taxids.xlsx'
gz_path = "/Volumes/Coruscant/APSCALE_raw_databases/accession2taxid/nucl_gb.accession2taxid.gz"
db_path = '/Volumes/Coruscant/APSCALE_raw_databases/accession2taxid/db/nucl_gb.accession2taxid.db'

for fasta_file in files:
    ## Run
    if Path(fasta_file).suffix == '.zip':
        fasta_file = zip_to_gz(fasta_file)
    run_silva(output_path, fasta_file, taxonomy_xlsx)
    run_silva_filtered(output_path, fasta_file, taxonomy_xlsx)









