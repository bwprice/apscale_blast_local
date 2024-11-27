import argparse
import multiprocessing
import os
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
    APSCALE blast command line tool - v1.0.2 - 09/09/2024
    Usage examples:
    $ apscale_blast blastn -h
    $ apscale_blast blastn -database ./MIDORI2_UNIQ_NUC_GB259_srRNA_BLAST -query_fasta ./12S_apscale_ESVs.fasta
    $ apscale_blast filter -h
    $ apscale_blast filter -database ./MIDORI2_UNIQ_NUC_GB259_srRNA_BLAST -blastn_folder ./12S_apscale_ESVs_blastn
    """
    print(message)

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description='APSCALE blast v1.0.2')

    # Creating subcommands (blastn and filter)
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # BLASTn subcommand parser
    blastn_parser = subparsers.add_parser('blastn', help='Perform blastn search on a selected fasta file.')
    # Filter subcommand parser
    filter_parser = subparsers.add_parser('filter', help='Filter blastn results based on thresholds.')

    # Arguments for both blastn and filter
    common_args = ['-database', '-apscale_gui']
    for arg in common_args:
        blastn_parser.add_argument(arg, type=str, required=False, help=f'PATH to {arg.lstrip("-")}.')
        filter_parser.add_argument(arg, type=str, required=False, help=f'PATH to {arg.lstrip("-")}.')

    # Arguments specific to BLASTn
    blastn_parser.add_argument('-blastn_exe', type=str, default='blastn',
                               help='PATH to blast executable. [DEFAULT: blastn]')
    blastn_parser.add_argument('-query_fasta', type=str, help='PATH to fasta file.')
    blastn_parser.add_argument('-n_cores', type=int, default=multiprocessing.cpu_count() - 1,
                               help='Number of CPU cores to use. [DEFAULT: CPU count - 1]')
    blastn_parser.add_argument('-task', type=str, default='blastn',
                               help='Blastn task: blastn, megablast, or dc-megablast. [DEFAULT: blastn]')
    blastn_parser.add_argument('-out', type=str, default='./',
                               help='PATH to output directory. A new folder will be created here. [DEFAULT: ./]')
    blastn_parser.add_argument('-subset_size', type=int, default=100,
                               help='Number of sequences per query fasta subset. [DEFAULT: 100]')
    blastn_parser.add_argument('-max_target_seqs', type=int, default=20,
                               help='Number of hits retained from the blast search. Larger values increase runtimes and storage needs. [DEFAULT: 20]')
    blastn_parser.add_argument('-masking', type=str, default='Yes',
                               help='Activate masking [DEFAULT="Yes"]')

    # Arguments specific to filter
    filter_parser.add_argument('-blastn_folder', type=str, help='PATH to blastn results folder for filtering.')
    filter_parser.add_argument('-thresholds', type=str, default='97,95,90,87,85',
                               help='Taxonomy filter thresholds. [DEFAULT: 97,95,90,87,85]')
    filter_parser.add_argument('-n_cores', type=int, default=multiprocessing.cpu_count() - 1,
                               help='Number of CPU cores to use. [DEFAULT: CPU count - 1]')

    # Parse the arguments
    args = parser.parse_args()

    # Handle missing arguments interactively for both commands
    if args.command == 'blastn' and not args.database and not args.query_fasta:
        args.database = input("Please enter PATH to database: ").strip('"')
        args.query_fasta = input("Please enter PATH to query fasta: ").strip('"')

        # Set output directory if default value is used
        if args.out == './':
            args.out = str(args.query_fasta).replace('.fasta', '')
            if not os.path.isdir(args.out):
                os.mkdir(Path(args.out))  # Create the output directory

    elif args.command == 'filter' and not args.database and not args.blastn_folder:
        args.database = input("Please enter PATH to database: ").strip('"')
        args.blastn_folder = input("Please enter PATH to blastn folder: ").strip('"')

    # Handle the 'blastn' command
    if args.command == 'blastn':
        if args.query_fasta:
            project_folder = args.out  # Use the output directory specified by the user
            # Run the BLASTn function
            a_blastn(args.blastn_exe,
                     args.query_fasta.strip('"'),
                     args.database.strip('"'),
                     project_folder,
                     args.n_cores,
                     args.task,
                     args.subset_size,
                     args.max_target_seqs,
                     args.masking,
                     args.apscale_gui)
        else:
            print('\nError: Please provide a fasta file!')

    # Handle the 'filter' command
    elif args.command == 'filter':
        if args.blastn_folder:
            # Run the filter function
            b_filter(args.blastn_folder.strip('"'), args.database.strip('"'), args.thresholds, args.n_cores)
        else:
            print('\nError: Please provide a blastn results file folder (.csv)!')


# Run the main function if script is called directly
if __name__ == "__main__":
    main()