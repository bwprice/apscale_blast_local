import glob, re
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import zipfile, gzip, io
from io import StringIO
import subprocess, shutil, os, datetime
from tqdm import tqdm

def run_diat_barcode(output_path, diat_barcode_xlsx):

    print('{} : Starting to collect accession numbers from xlsx file.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    diat_barcode_df = pd.read_excel(diat_barcode_xlsx, sheet_name='diatbarcode v12').fillna('')

    print('{} : Writing fasta file.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    fasta_file = diat_barcode_xlsx.replace('.xlsx', '.fasta.gz').replace(' ', '_')
    with gzip.open(fasta_file, 'wt') as f:
        for line in diat_barcode_df[['Sequence ID', 'Sequence']].values.tolist():
            if line[0] != '':
                f.write(f'>{line[0]}\n')
                f.write(f'{line[1]}\n')

    print('{} : Finished to convert fasta format.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    print('{} : Starting to generate taxonomy file.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    records = []
    for line in diat_barcode_df[["Species", "Genus", "Family (following Round, Crawford & Mann 1990)", "Order (following Round, Crawford & Mann 1990)", "Class (following Round, Crawford & Mann 1990)", "Phylum (following Algaebase 2018)", "Subkingdom (following Algaebase 2018)", "Sequence ID"]].values.tolist():
        if line[0] != '':
            records.append(line[::-1])
    records_df = pd.DataFrame(records, columns=['Accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

    print('{} : Finished to convert accession numbers to taxonomy.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    print('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    # create database
    db_name = Path(fasta_file).name.replace('.fasta.gz', '').replace('.', '_')
    db_folder = Path(output_path).joinpath(f'db_{db_name}')
    if not os.path.isdir(db_folder):
        os.mkdir(db_folder)
    db_folder = db_folder.joinpath('db')
    command = f'zcat < {fasta_file} | makeblastdb -in - -title db -dbtype nucl -out {db_folder}'
    os.system(command)

    # move taxonomy file
    taxonomy_file = db_folder.parent.joinpath('db_taxonomy.parquet.snappy')
    records_df.to_parquet(taxonomy_file)

    ## zip the folder
    output = Path(output_path).joinpath(f'db_{db_name}')
    shutil.make_archive(output, 'zip', output)

    print('{} : Finished to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

## Variables
diat_barcode_xlsx = '/Volumes/Coruscant/APSCALE_raw_databases/2024_09_fasta/DiatBarcode/2024-05-28-Diat.barcode_release-version 12.4.xlsx'
output_path = '/Volumes/Coruscant/APSCALE_raw_databases/2024_09'

# Run
run_diat_barcode(output_path, diat_barcode_xlsx)



