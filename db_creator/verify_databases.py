import numpy as np
import gzip, glob, os
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import time, subprocess
from tqdm import tqdm
from plotly.subplots import make_subplots
import plotly.graph_objects as go

def create_test_fasta(output_path):
    # #create test folder
    test_folder = f'{output_path}_verification'
    os.makedirs(test_folder, exist_ok=True)

    ## collect all fasta files
    fasta_folder = f'{output_path}_fasta'
    files = glob.glob(f'{fasta_folder}/*/*.fasta.gz')

    ## create an new test fasta file
    entries = {}
    for file in files:
        name = Path(file).name.replace('.fasta.gz', '')
        with gzip.open(file, "rt") as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fasta")):
                if i < 10:
                    id = f'>{name}>>{record.id}\n'
                    seq = f'{str(record.seq)}\n'
                    entries[id] = seq
                else:
                    break

    ## write test fasta file
    with open(f'{output_path}_verification/test.fasta', 'w') as f:
        for key, values in entries.items():
            f.write(key)
            f.write(values)

def run_apscale_blast(output_path):

    ## collect files
    test_folder = f'{output_path}_verification'
    test_fasta = f'{output_path}_verification/test.fasta'
    databases = sorted([i.replace('.zip', '') for i in glob.glob(f'{output_path}/*.zip')])
    timings = {}

    # db = '/Volumes/Coruscant/APSCALE_raw_databases/2024_09/db_MIDORI2_UNIQ_NUC_SP_GB261_srRNA_BLAST'

    for db in tqdm(databases):
        ## time the blast search
        start_time = time.time()

        ## define output folder
        db_name = Path(db).name
        output_folder = Path(test_folder).joinpath(db_name)
        os.makedirs(output_folder, exist_ok=True)

        ## perform blastn
        command = f'apscale_blast blastn -database {db} -query_fasta {test_fasta} -out {output_folder} -subset_size 20'
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        # Print the output and error (if any)
        print("Output:", result.stdout)
        print("Error:", result.stderr)

        ## perform filtering
        command = f'apscale_blast filter -database {db} -blastn_folder {output_folder}'
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        # Print the output and error (if any)
        print("Output:", result.stdout)
        print("Error:", result.stderr)

        # Measure the execution time
        end_time = time.time()
        timings[db_name] = [round(end_time - start_time,2), round((end_time - start_time) / 60, 2)]

    return timings

def fasta_count(fasta_file):
    with gzip.open(fasta_file, "rt") as handle:
        return sum(1 for _ in SeqIO.parse(handle, "fasta"))

## Variables
output_path = '/Volumes/Coruscant/APSCALE_raw_databases/2024_09'

## Run
create_test_fasta(output_path)
timings = run_apscale_blast(output_path)

## File sizes
files = sorted(glob.glob('/Volumes/Coruscant/APSCALE_raw_databases/2024_09_fasta/*/*.fasta.gz'))
n_sequences = [fasta_count(fasta_file) for fasta_file in files]

## benchmark
x_values = ['diat.barcode 12.4',
 'MIDORI2_A6',
 'MIDORI2_A8',
 'MIDORI2_CO1',
 'MIDORI2_CO2',
 'MIDORI2_CO3',
 'MIDORI2_cytb',
 'MIDORI2_ND1',
 'MIDORI2_ND2',
 'MIDORI2_ND3',
 'MIDORI2_ND4L',
 'MIDORI2_ND4',
 'MIDORI2_ND5',
 'MIDORI2_ND6',
 'MIDORI2_lrRNA',
 'MIDORI2_srRNA',
 'SILVA_138.2_LSURef_NR99',
 'SILVA_138.2_SSURef_NR99',
 'PR2_5.0.0_SSU',
 'UNITE_all_eukaryotes',
 'UNITE_fungi']
y_values = [i[1] for i in list(timings.values())]

fig = make_subplots(rows=2, cols=1, subplot_titles=['Database size', 'Benchmark (100 sequences)'], shared_xaxes=True)

fig.add_trace(go.Bar(x=x_values, y=n_sequences), row=1,col=1)
fig.update_yaxes(title='sequences', row=1,col=1)

fig.add_trace(go.Bar(x=x_values, y=y_values), row=2,col=1)
fig.update_yaxes(title='run time (min)', row=2,col=1)

fig.update_layout(template='simple_white', width=600, height=600, showlegend=False)
fig.write_image(f'{output_path}_verification/timings.png')