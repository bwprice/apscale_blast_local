import argparse
import multiprocessing
import os
import sys
import datetime
from pathlib import Path
from apscale_blast.a_blastn import main as a_blastn
from apscale_blast.b_filter import main as b_filter

def main():
    """
    APSCALE BLASTn suite
    Command-line tool to run and filter BLASTn searches.
    """

    # Introductory message with usage examples
    message = """
    APSCALE blast command line tool - v1.3.0 (HPC-compatible)
    Example commands:
    $ apscale_blast -h
    $ apscale_blast -db ./MIDORI2_UNIQ_NUC_GB259_srRNA_BLAST -q ./12S_apscale_ESVs.fasta
    """
    print(message)

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description='APSCALE blast v1.3.0 (HPC-compatible)')

    # === Main settings ===
    main_settings = parser.add_argument_group("Main settings")
    main_settings.add_argument('-database', '-db', type=str, required=False, help='PATH to local database.')
    main_settings.add_argument('-query_fasta', '-q', type=str, help='PATH to fasta file.')
    main_settings.add_argument('-out', '-o', type=str, default='./', help='PATH to output directory. A new folder will be created here. [DEFAULT: ./]')

    # === General settings ===
    general_settings = parser.add_argument_group("General settings")
    general_settings.add_argument('-n_cores', type=int, default=multiprocessing.cpu_count() - 1, help='Number of CPU cores to use. [DEFAULT: CPU count - 1]')
    general_settings.add_argument('-task', type=str, default='blastn', help='Blastn task: blastn, megablast, or dc-megablast. [DEFAULT: blastn]')
    general_settings.add_argument('-subset_size', type=int, default=100, help='Number of sequences per query fasta subset. [DEFAULT: 100]')
    general_settings.add_argument('-max_target_seqs', type=int, default=20, help='Number of hits retained from the blast search. [DEFAULT: 20]')
    general_settings.add_argument('-thresholds', type=str, default='97,95,90,87,85', help='Taxonomy filter thresholds. [DEFAULT: 97,95,90,87,85]')

    # === Advanced settings ===
    advanced_settings = parser.add_argument_group("Advanced settings")
    advanced_settings.add_argument('-blastn_exe', type=str, default='blastn', help='PATH to blast executable. [DEFAULT: blastn]')
    advanced_settings.add_argument('-masking', type=str, default='Yes', help='Activate masking. [DEFAULT="Yes"]')
    advanced_settings.add_argument('-gui', action='store_true', help='Only required for Apscale-GUI.')

    # Parse the arguments
    args = parser.parse_args()

    # Handle missing arguments interactively for both commands
    if not args.database and not args.query_fasta:
        args.database = input("Please enter PATH to database: ").strip('"')
        args.query_fasta = input("Please enter PATH to query fasta: ").strip('"')

        # Set output directory if default value is used
        if args.out == './':
            args.out = str(args.query_fasta).replace('.fasta', '')
            if not os.path.isdir(args.out):
                os.mkdir(Path(args.out))  # Create the output directory

    ## CHECK IF FILES ALREADY EXIST
    project_folder = args.out  # Use the output directory specified by the user

    # Handle the 'blastn' command
    continue_blast = False

    # Convert db to Path
    database = Path(args.database.strip('"'))

    if args.query_fasta:
        # Run the BLASTn function
        continue_blast = a_blastn(args.blastn_exe,
                 args.query_fasta.strip('"'),
                 database,
                 project_folder,
                 args.n_cores,
                 args.task,
                 args.subset_size,
                 args.max_target_seqs,
                 args.masking,
                 False,  # headless (no longer used)
                 args.gui,
                 [],     # organism_mask (no longer used)
                 False   # include_uncultured (no longer used)
                                  )
    else:
        print('\nError: Please provide a fasta file!')
    # Handle the 'filter' command
    if continue_blast == False:
        print('\nNot all fasta subsets have been processed yet!')
    elif not os.path.isfile(Path(project_folder).joinpath('log.txt')):
        print('\nError: Could not find the BLAST results folder!')
    else:
        # Run the filter function
        b_filter(project_folder, database, args.thresholds, args.n_cores)

# Run the main function if script is called directly
if __name__ == "__main__":
    main()