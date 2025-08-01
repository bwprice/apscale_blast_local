import glob, re, gzip
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import zipfile, gzip, io
from io import StringIO
import subprocess, shutil, os, datetime
from tqdm import tqdm

def trnl_taxonomy(taxonomy_file):
    accession_df = pd.read_csv(taxonomy_file, sep='[;\t]', engine='python').fillna('')
    accession_df.columns = ['Accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    return accession_df

def run_trnl(output_path, fasta_file, taxonomy_file):

    print('{} : Starting to collect taxonomy from fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## collect accession numbers
    accession_df = trnl_taxonomy(taxonomy_file)

    print('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    # create database
    db_name = Path(fasta_file).name.replace('.fasta.gz', '')
    db_folder = Path(output_path).joinpath(f'db_{db_name}')
    os.makedirs(db_folder, exist_ok=True)
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
output_path = '/Users/tillmacher/Documents/GitHub/APSCALE_database_creator/db/trnl'
fasta_file = '/Users/tillmacher/Documents/GitHub/APSCALE_database_creator/db/trnL/trnL.fasta.gz'
taxonomy_file = '/Users/tillmacher/Documents/GitHub/APSCALE_database_creator/db/trnL/trnL_taxonomy.txt'
run_trnl(output_path, fasta_file, taxonomy_file)





























