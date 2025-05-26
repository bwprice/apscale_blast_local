import glob, re
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import zipfile, gzip, io
from io import StringIO
import subprocess, shutil, os, datetime
from tqdm import tqdm
import zipfile
import gzip
import shutil
import os


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

def pr2_taxonomy(fasta_file):
    all_accession_numbers = []
    with gzip.open(fasta_file, 'rt') as myfile:
        data = myfile.read()
        sequences = SeqIO.parse(StringIO(data), 'fasta')
        for record in sequences:
            accession = record.id
            taxonomy_dict = {i.split(':')[0]:i.split(':')[1] for i in record.id.split(';')[1].split('=')[1].split(',')}
            taxonomy = []
            for taxon in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
                try:
                    taxonomy.append(taxonomy_dict[taxon].replace('_', ' '))
                except KeyError:
                    taxonomy.append('')
            all_accession_numbers.append([accession] + taxonomy)

    accession_df = pd.DataFrame(all_accession_numbers, columns=['Accession', 'superkingdom','phylum', 'class', 'order', 'family', 'genus', 'species'])
    return accession_df

def run_pr2(output_path, fasta_file):

    print('{} : Starting to collect taxonomy from fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## collect accession numbers
    accession_df = pr2_taxonomy(fasta_file)

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

## Variables
output_path = '/Volumes/Coruscant/APSCALE_raw_databases/2024_09'
files = glob.glob('/Volumes/Coruscant/APSCALE_raw_databases/2024_09_fasta/PR2/*.fasta*')

for fasta_file in files:
    ## Run
    if Path(fasta_file).suffix == '.zip':
        fasta_file = zip_to_gz(fasta_file)
    run_pr2(output_path, fasta_file)