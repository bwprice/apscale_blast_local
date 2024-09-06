import glob, re, gzip
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import zipfile, gzip, io
from io import StringIO
import subprocess, shutil, os, datetime
from tqdm import tqdm

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

def silva_taxonomy(fasta_file):
    all_accession_numbers = []
    with gzip.open(fasta_file, 'rt') as myfile:
        data = myfile.read()
        sequences = SeqIO.parse(StringIO(data), 'fasta')
        for record in sequences:
            accession = record.id
            taxonomy = record.description.split(' ')[1].split(';')
            n_taxa = len(taxonomy)
            if n_taxa == 7:
                ## correct taxonomy. TAKE IT!
                all_accession_numbers.append([accession] + taxonomy)
            elif n_taxa < 7:
                ## add missing levels
                missing = 7 - n_taxa
                taxonomy = taxonomy + [''] * missing
                all_accession_numbers.append([accession] + taxonomy)
            else:
                taxonomy = taxonomy[:7]
                all_accession_numbers.append([accession] + taxonomy)

    accession_df = pd.DataFrame(all_accession_numbers, columns=['Accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

    return accession_df

def run_silva(output_path, fasta_file):

    print('{} : Starting to collect taxonomy from fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## collect accession numbers
    accession_df = silva_taxonomy(fasta_file)

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

## Variables
output_path = '/Volumes/Coruscant/APSCALE_raw_databases/2024_09'
files = glob.glob('/Volumes/Coruscant/APSCALE_raw_databases/2024_09_fasta/SILVA/*.fasta*')

for fasta_file in files:
    ## Run
    if Path(fasta_file).suffix == '.zip':
        fasta_file = zip_to_gz(fasta_file)
    run_silva(output_path, fasta_file)





























