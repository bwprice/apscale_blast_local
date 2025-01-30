import argparse
import multiprocessing
import os
import sys
from pathlib import Path
from apscale_blast.a_blastn import main as a_blastn
from apscale_blast.b_filter import main as b_filter
from ete3 import NCBITaxa, Tree

def main():
    """
    APSCALE BLASTn suite
    Command-line tool to run and filter BLASTn searches.
    """

    # Introductory message with usage examples
    message = """
    APSCALE blast command line tool - v1.2.0
    Example commands:
    $ apscale_blast -h
    $ apscale_blast -db ./MIDORI2_UNIQ_NUC_GB259_srRNA_BLAST -q ./12S_apscale_ESVs.fasta
    
    Remember to update your local ete3 NCBI taxonomy regularly, if using the "remote" blastn!
    This can be performed by running:
    $ apscale_blast -u
    """
    print(message)

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description='APSCALE blast v1.1.0')

    # Arguments for both blastn and filter
    parser.add_argument('-database', '-db', type=str, required=False, help=f'PATH to local database. Use "remote" to blast against the complete GenBank database (might be slow)')
    parser.add_argument('-blastn_exe', type=str, default='blastn', help='PATH to blast executable. [DEFAULT: blastn]')
    parser.add_argument('-query_fasta', '-q', type=str, help='PATH to fasta file.')
    parser.add_argument('-n_cores', type=int, default=multiprocessing.cpu_count() - 1, help='Number of CPU cores to use. [DEFAULT: CPU count - 1]')
    parser.add_argument('-task', type=str, default='blastn', help='Blastn task: blastn, megablast, or dc-megablast. [DEFAULT: blastn]')
    parser.add_argument('-out', '-o', type=str, default='./', help='PATH to output directory. A new folder will be created here. [DEFAULT: ./]')
    parser.add_argument('-subset_size', type=int, default=100, help='Number of sequences per query fasta subset. [DEFAULT: 100]')
    parser.add_argument('-max_target_seqs', type=int, default=20, help='Number of hits retained from the blast search. Larger values increase runtimes and storage needs. [DEFAULT: 20]')
    parser.add_argument('-masking', type=str, default='Yes', help='Activate masking [DEFAULT="Yes"]')
    parser.add_argument('-thresholds', type=str, default='97,95,90,87,85', help='Taxonomy filter thresholds. [DEFAULT: 97,95,90,87,85]')
    # parser.add_argument('-headless', type=str, default="True", help='Display the Chromium browser during the remote blast [DEFAULT: True')
    parser.add_argument('-update_taxids', '-u', action='store_true', help='Update NCBI taxid backbone')
    parser.add_argument('-gui', action='store_true', help='Only required for Apscale-GUI')

    # Parse the arguments
    args = parser.parse_args()

    # Handle taxonomy update
    if args.update_taxids:
        print("Updating NCBI taxonomy database...")
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()
        print("Taxonomy database updated successfully.")
        if Path('./taxdump.tar.gz').exists():
            os.remove('./taxdump.tar.gz')
            print('Removed taxdmup.tar.gz')
        return  # Exit after updating taxonomy

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
    headless = "True"
    continue_blast = False
    if args.query_fasta:
        # Run the BLASTn function
        continue_blast = a_blastn(args.blastn_exe,
                 args.query_fasta.strip('"'),
                 args.database.strip('"'),
                 project_folder,
                 args.n_cores,
                 args.task,
                 args.subset_size,
                 args.max_target_seqs,
                 args.masking,
                 headless,
                 args.gui)
    else:
        print('\nError: Please provide a fasta file!')

    # Handle the 'filter' command
    if continue_blast == False:
        print('\nNot all fasta subsets have been processed yet!')
    elif not os.path.isfile(Path(project_folder).joinpath('log.txt')):
        print('\nError: Could not find the BLAST results folder!')
    else:
        # Run the filter function
        b_filter(project_folder, args.database.strip('"'), args.thresholds, args.n_cores)

# Run the main function if script is called directly
if __name__ == "__main__":
    main()