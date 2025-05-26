import argparse
import os
from pathlib import Path
import datetime
import subprocess
import shutil
import time
import glob
import multiprocessing
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import threading
import pyarrow as pa
import pyarrow.parquet as pq
from ete3 import NCBITaxa, Tree
import requests
import xmltodict
from tqdm import tqdm

# Lock for thread-safe print statements
print_lock = threading.Lock()

# Ensure multiprocessing works when frozen (e.g. for executables)
multiprocessing.freeze_support()

def accession2taxonomy(df_1, taxid_dict, col_names_2, db_name):
    """
    Convert accession IDs to taxonomy information.
    """
    df_2_list = []
    for row in df_1.values.tolist():
        ID_name = row[0]
        accession = row[1]

        evalue = row[-1]
        similarity = row[-2]
        try:
            taxonomy = taxid_dict[accession]
        except KeyError:
            taxonomy = ['NoMatch'] * 7
        df_2_list.append([ID_name] + taxonomy + [similarity, evalue])
    df_2 = pd.DataFrame(df_2_list, columns=col_names_2)
    return df_2

def ncbi_taxid_request(taxid):
    desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    try:
        ncbi = NCBITaxa()
        lineage = ncbi.get_lineage(taxid)
        lineage2ranks = ncbi.get_rank(lineage)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        results = {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}
        taxids = [str(taxid) for taxid in list(results.values())]
        taxonomy = [list(ncbi.get_taxid_translator([taxid]).values())[0] for taxid in taxids]
        return taxonomy
    except ValueError:
        taxonomy_placeholder = f'Unknown taxid {str(taxid)}'
        return [taxonomy_placeholder] * 7

def filter_blastn_csvs(file, taxid_dict, i, n_subsets, thresholds, db_name):
    """
    Filter BLASTn results from CSV files.
    """

    ## load blast results
    col_names = ['unique ID', 'Sequence ID', 'Similarity', 'evalue']

    if not os.path.isfile(file):
        print('{}: Skipping missing subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))
        return
    if os.path.getsize(file) == 0:
        print('{}: Skipping empty subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))
        return

    csv_df = pd.read_csv(file, header=None, sep=';;', names=col_names, engine='python').fillna('NAN')
    csv_df['Similarity'] = [float(i) for i in csv_df['Similarity'].values.tolist()]

    if len(csv_df) == 0:
        print('{}: Error during filtering for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))
        return

    thresholds = thresholds.split(',')
    # Ensure correct number of thresholds
    if len(thresholds) != 5:
        print('Please provide 5 comma-separated threshold values!')
        print('Using default values...')
        thresholds = ['97', '95', '90', '87', '85']

    species_threshold = int(thresholds[0])
    genus_threshold = int(thresholds[1])
    family_threshold = int(thresholds[2])
    order_threshold = int(thresholds[3])
    class_threshold = int(thresholds[4])

    ## filter hits
    ID_set = csv_df['unique ID'].drop_duplicates().values.tolist()
    col_names_2 = ['unique ID', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity', 'evalue']
    taxonomy_df = pd.DataFrame()

    # loop through IDs
    # ID = ID_set[0]
    for ID in ID_set:
        ## filter by evalue
        df_0 = csv_df.loc[csv_df['unique ID'] == ID]
        max_sim = max(df_0['Similarity'])

        ##### 1) FILTERING BY SIMILARITY THEN BY EVALUE (or the other way round, let the user decide)
        df_1 = df_0.loc[df_0['Similarity'] == max_sim]
        max_e = max(df_1['evalue'])
        df_1 = df_1.loc[df_1['evalue'] == max_e]

        ############################################################################################################
        ## convert fasta headers to taxonomy
        if db_name != 'remote':
            df_2 = accession2taxonomy(df_1, taxid_dict, col_names_2, db_name)
        else:
            df_2_values = []
            for row in df_1.values.tolist():
                taxid = row[1]
                taxonomy = ncbi_taxid_request(taxid)
                df_2_values.append([row[0]] + taxonomy + row[2:])
            df_2 = pd.DataFrame(df_2_values, columns=col_names_2)

        ## Filter out missing taxids
        if df_2['Species'].str.contains(r'\bUnknown taxid\b', na=False).any():
            df_2_reduced = df_2.loc[~df_2['Species'].str.contains(r'\bUnknown taxid\b', na=False)]
            if len(df_2_reduced) != 0:
                df_2 = df_2_reduced.copy()

        ############################################################################################################

        ##### 2) ROBUSTNESS TRIMMING BY SIMILARITY
        n_hits = len(df_2)

        if max_sim >= species_threshold:
            pass
        elif max_sim < species_threshold and max_sim >= genus_threshold:
            df_2['Species'] = ['']*n_hits
        elif max_sim < genus_threshold and max_sim >= family_threshold:
            df_2['Species'] = ['']*n_hits
            df_2['Genus'] = ['']*n_hits
        elif max_sim < family_threshold and max_sim >= order_threshold:
            df_2['Species'] = ['']*n_hits
            df_2['Genus'] = ['']*n_hits
            df_2['Family'] = ['']*n_hits
        elif max_sim < order_threshold and max_sim >= class_threshold:
            df_2['Species'] = ['']*n_hits
            df_2['Genus'] = ['']*n_hits
            df_2['Family'] = ['']*n_hits
            df_2['Order'] = ['']*n_hits
        else:
            df_2['Species'] = ['']*n_hits
            df_2['Genus'] = ['']*n_hits
            df_2['Family'] = ['']*n_hits
            df_2['Order'] = ['']*n_hits
            df_2['Class'] = ['']*n_hits

        ##### 3) REMOVAL OF DUPLICATE HITS
        df_3 = df_2.drop_duplicates().copy()

        ##### 4) MORE THAN ONE TAXON REMAINING?
        if len(df_3) == 1:
            df_3['Flag'] = ['']*len(df_3)
            df_3['Ambigous taxa'] = [''] * len(df_3)
            taxonomy_df = pd.concat([taxonomy_df, df_3])

        ##### 5) SPECIES LEVEL REFERENCE?
        elif max_sim < species_threshold:
            # remove taxonomic levels until a single hit remains
            for level in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'][::-1]:
                n_hits = len(df_3)
                df_3.loc[:, level] = [''] * n_hits
                df_3 = df_3.drop_duplicates()
                if len(df_3) == 1:
                    break
            df_3['Flag'] = [''] * len(df_3)
            df_3['Ambigous taxa'] = [''] * len(df_3)
            taxonomy_df = pd.concat([taxonomy_df, df_3])

        ##### 6) AMBIGOUS SPECIES WORKFLOW (FLAGGING SYSTEM)
        else:
            ##### 7) DOMINANT SPECIES PRESENT? (F1)
            df_2['duplicate_count'] = df_2.groupby(df_2.columns.tolist()).transform('size')
            df_2_dominant = df_2.loc[df_2['duplicate_count'] == max(df_2['duplicate_count'])].drop_duplicates()

            if len(df_2_dominant) == 1:
                df_3 = df_2_dominant.drop(columns=['duplicate_count'])
                df_3['Flag'] = ['F1'] * len(df_3)
                df_3['Ambigous taxa'] = [', '.join(df_2['Species'].drop_duplicates().values.tolist())] * len(df_3)
                taxonomy_df = pd.concat([taxonomy_df, df_3])

            else:
                ##### 8) TWO SPECIES OF ONE GENUS? (F2)
                n_genera = len(set(df_3['Genus']))
                if n_genera == 1 and len(df_3) == 2:
                    genus = df_3['Genus'].values.tolist()[0]
                    species = '/'.join([i.replace(genus + ' ', '') for i in df_3['Species'].values.tolist()])
                    df_3['Species'] = ['{} {}'.format(genus, species)]*len(df_3)
                    df_3['Flag'] = ['F2'] * len(df_3)
                    df_3['Ambigous taxa'] = [', '.join(df_2['Species'].drop_duplicates().values.tolist())] * len(df_3)
                    taxonomy_df = pd.concat([taxonomy_df, df_3])

                ##### 9) MULTIPLE SPECIES OF ONE GENUS? (F3)
                elif n_genera == 1 and len(df_3) != 2:
                    genus = df_3['Genus'].values.tolist()[0]
                    df_3['Species'] = [genus + ' sp.']*len(df_3)
                    df_3['Flag'] = ['F3'] * len(df_3)
                    df_3['Ambigous taxa'] = [', '.join(df_2['Species'].drop_duplicates().values.tolist())] * len(df_3)
                    taxonomy_df = pd.concat([taxonomy_df, df_3])

                ##### 10) TRIMMING TO MOST RECENT COMMON TAXON (F4)
                else:
                    # remove taxonomic levels until a single hit remains
                    for level in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'][::-1]:
                        n_hits = len(df_3)
                        df_3.loc[:, level] = ['']*n_hits
                        df_3 = df_3.drop_duplicates()
                        if len(df_3) == 1:
                            break
                    df_3['Flag'] = ['F4'] * len(df_3)
                    df_3['Ambigous taxa'] = [', '.join(df_2['Species'].drop_duplicates().values.tolist())] * len(df_3)
                    taxonomy_df = pd.concat([taxonomy_df, df_3])

    # export dataframe
    taxonomy_df.columns = ['unique ID', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity', 'evalue', 'Flag', 'Ambigous taxa']
    blastn_filtered_xlsx = file.replace('.csv', '_filtered.xlsx')
    taxonomy_df.to_excel(blastn_filtered_xlsx, sheet_name='Taxonomy table', index=False)

    print('{}: Finished filtering for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))

def filter_blast_csvs_dbDNA(file, i, n_subsets, thresholds):
    """
    Filter BLASTn results from CSV files for the dbDNA database.
    """
    ## load results
    col_names = ['unique ID', 'Sequence ID', 'Similarity', 'evalue']
    blast_df = pd.read_csv(file, header=None, sep=';;', names=col_names, engine='python').fillna('NAN')
    blast_df['Similarity'] = [float(i) for i in blast_df['Similarity'].values.tolist()]
    all_OTUs = blast_df['unique ID'].drop_duplicates().tolist()
    rating_snappy = Path('/Volumes/Coruscant/dbDNA/FEI_genera_v2_BarCodeBank/3_BarCodeBank/FEI_genera_v2.BarCodeBank.parquet.snappy')
    rating_df = pd.read_parquet(rating_snappy)

    ## collect information about the hits
    reference_list = []
    for OTU in all_OTUs:
        tmp = blast_df.loc[blast_df['unique ID'] == OTU].sort_values('evalue', ascending=True)
        max_similarity = max(tmp['Similarity'])
        tmp = tmp.loc[tmp['Similarity'] == max_similarity]
        hits = tmp['Sequence ID'].values.tolist()[:20]

        for hit in hits:
            sequenceID = hit.split('__')[0].replace('>', '')
            processid = hit.split('__')[1]
            # species_name = hit.split('__')[2].replace('_', ' ')
            reference = rating_df.loc[(rating_df['sequenceID'] == sequenceID) & (rating_df['processid'] == processid)].values.tolist()
            if reference != []:
                reference_list.append([OTU, max_similarity] + reference[0])
            else:
                reference_list.append([OTU, max_similarity] + reference)

    # store intermediate results
    reference_df = pd.DataFrame(reference_list, columns=['unique ID', 'Similarity'] + rating_df.columns.tolist())

    # filter results
    gold_threshold = 40
    silver_threshold = 25
    bronze_threshold = 10
    species_threshold = 97
    genus_threshold = 94
    family_threshold = 91
    order_threshold = 88
    class_threshold = 85

    blast_filtered_list = []

    for OTU in all_OTUs:
        tmp = reference_df.loc[reference_df['unique ID'] == OTU].copy()
        species = tmp['species_name']
        n_hits = len(species)
        similarity = tmp['Similarity'].drop_duplicates().values.tolist()[0]

        # 1) trim taxonomy according to similarity
        if similarity >= species_threshold:
            pass
        elif similarity < species_threshold and similarity >= genus_threshold:
            tmp['species_name'] = ['']*n_hits
        elif similarity < genus_threshold and similarity >= family_threshold:
            tmp['species_name'] = ['']*n_hits
            tmp['genus_name'] = [''] * n_hits
        elif similarity < family_threshold and similarity >= order_threshold:
            tmp['species_name'] = ['']*n_hits
            tmp['genus_name'] = [''] * n_hits
            tmp['family_name'] = [''] * n_hits
        elif similarity < order_threshold and similarity >= class_threshold:
            tmp['species_name'] = ['']*n_hits
            tmp['genus_name'] = [''] * n_hits
            tmp['family_name'] = [''] * n_hits
            tmp['order_name'] = [''] * n_hits
        else:
            tmp['species_name'] = ['']*n_hits
            tmp['genus_name'] = [''] * n_hits
            tmp['family_name'] = [''] * n_hits
            tmp['order_name'] = [''] * n_hits
            tmp['class_name'] = [''] * n_hits

        # 2) trim taxonomy according to rating
        max_rating = max(tmp['rating'])
        if max_rating >= gold_threshold:
            tmp = tmp.loc[tmp['rating'] >= gold_threshold].copy()
            rating_str = 'A - Gold'
        elif max_rating >= silver_threshold:
            tmp = tmp.loc[tmp['rating'] >= silver_threshold].copy()
            rating_str = 'B - Silver'
        elif max_rating >= bronze_threshold:
            tmp = tmp.loc[tmp['rating'] >= bronze_threshold].copy()
            rating_str = 'C - Bronze'
        else:
            tmp = tmp.loc[tmp['rating'] < bronze_threshold].copy()
            rating_str = 'D - unreliable'

        # Export hit
        blast_hit = []
        relevant_columns = ['unique ID', 'Similarity',
                            'rating', 'bin_uri',
                            'phylum_name', 'class_name',
                            'order_name', 'family_name',
                            'genus_name', 'species_name',
                            'phylogeny', 'species_group',
                            'identification_by', 'institution_storing',
                            'country', 'province',
                            'region', 'exactsite',
                            'lifestage', 'sex']

        for col in relevant_columns:
            res = tmp[col].drop_duplicates().values.tolist()

            # calculate average rating
            if col == 'rating':
                blast_hit.append(np.average(res))

            # trim taxonomy if necessary and keep all information in separate cell
            elif '_name' in col:
                if len(res) != 1:
                    blast_hit.append('')
                    str_list = ', '.join(map(str, res))
                    blast_hit.append(str_list)
                else:
                    blast_hit.append(res[0])
                    blast_hit.append(res[0])
            ## merge all other information
            elif len(res) != 1:
                str_list = ', '.join(map(str, res))
                blast_hit.append(str_list)
            else:
                blast_hit.append(res[0])

        blast_filtered_list.append(blast_hit + [rating_str])

    # create a dataframe
    columns = ['unique ID', 'Similarity',
                'rating', 'bin_uri',
                'Phylum', 'all_phyla',
                'Class', 'all_classes',
                'Order', 'all_orders',
                'Family', 'all_families',
                'Genus', 'all_genera',
                'Species', 'all_species',
                'phylogeny', 'species_group',
                'identification_by', 'institution_storing',
                'country', 'province',
                'region', 'exactsite',
                'lifestage', 'sex', 'Standard']

    blast_filtered_df = pd.DataFrame(blast_filtered_list, columns=columns)

    # remove BINs of hit is not on species level
    blast_filtered_df.loc[blast_filtered_df['Species'] == '', 'bin_uri'] = ''

    # sort dataframe to TaXon table format compatible table
    sorted_columns = ['unique ID', 'Phylum',
                'Class', 'Order',
                'Family', 'Genus',
                'Species', 'Similarity',
                'rating', 'Standard', 'bin_uri',
                'all_phyla', 'all_classes',
                'all_orders', 'all_families',
                'all_genera', 'all_species',
                'phylogeny', 'species_group',
                'identification_by', 'institution_storing',
                'country', 'province',
                'region', 'exactsite',
                'lifestage', 'sex']

    ## sort df
    blast_filtered_df_sorted = blast_filtered_df[sorted_columns]

    # write dataframe
    blastn_filtered_xlsx = file.replace('.csv', '_filtered.xlsx')
    blast_filtered_df_sorted.to_excel(blastn_filtered_xlsx, sheet_name='Taxonomy table', index=False)

    ## finish command
    print('{}: Finished filtering for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))

def main(blastn_folder, blastn_database, thresholds, n_cores):
    """ Filter results according to Macher et al., 2023 (Fish Mock Community paper) """

    print('{}: Starting to filter blast results for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), blastn_folder))
    print('{}: Your database: {}'.format(datetime.datetime.now().strftime('%H:%M:%S'),  Path(blastn_database).stem))

    ## load blast results
    csv_files = glob.glob('{}/subsets/*.csv'.format(blastn_folder))

    ## check if a dbDNA database was used
    if "_dbDNA" in str(blastn_database):
        # PARALLEL FILTER COMMAND
        n_subsets = len(csv_files)
        Parallel(n_jobs = n_cores, backend='threading')(delayed(filter_blast_csvs_dbDNA)(file, i, n_subsets, thresholds) for i, file in enumerate(csv_files))

        ## also already define the no match row
        NoMatch = ["No Match"] * 6 + [0] * 2 + ['']*17
    else:

        if blastn_database != 'remote':
            ## load taxid table
            taxid_table = Path(blastn_database).joinpath('db_taxonomy.parquet.snappy')
            taxid_df = pd.read_parquet(taxid_table).fillna('')
            taxid_dict = {i[0]:i[1::] for i in taxid_df.values.tolist()}
            ## collect name of database
            db_name = Path(blastn_database).stem
        else:
            ## remote blast does not require the taxid_dict, but it must be defined anyways
            taxid_dict = {}
            ## collect name of databaseU
            db_name = "remote"

        # PARALLEL FILTER COMMAND
        n_subsets = len(csv_files)
        # file = csv_files[2]
        # i = 0
        if blastn_database != 'remote':
            Parallel(n_jobs = n_cores, backend='threading')(delayed(filter_blastn_csvs)(file, taxid_dict, i, n_subsets, thresholds, db_name) for i, file in enumerate(csv_files))
        else:
            [filter_blastn_csvs(file, taxid_dict, i, n_subsets, thresholds, db_name) for i, file in enumerate(csv_files)]

        ## also already define the no match row
        NoMatch = ['NoMatch'] * 7 + [0, 1, '', '']

    # Get a list of all the xlsx files
    xlsx_files = glob.glob('{}/subsets/*.xlsx'.format(blastn_folder))

    if len(xlsx_files) == 0:
        print('## Warning ##')
        print('{}: Error during filtering of the blast results for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), blastn_folder))
        print('{}: Please check your database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
        print('{}: Errors usually occur during unzipping of the database file.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
        print('## Warning ##')
    else:

        # Create a list to hold all the individual DataFrames
        df_list = []

        # Loop through the list of xlsx files
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Collecting {len(xlsx_files)} filtered subsets.')
        for file in tqdm(xlsx_files):
            # Read each xlsx file into a DataFrame
            df = pd.read_excel(file).fillna('')
            # Append the DataFrame to the list
            df_list.append(df)

        # Concatenate all the DataFrames in the list into one DataFrame
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Concatenating subsets.')
        merged_df = pd.concat(df_list, ignore_index=True)
        name = Path(blastn_folder).name
        blastn_filtered_xlsx = Path('{}/{}_taxonomy.xlsx'.format(blastn_folder, name))

        ## add OTUs without hit
        # Drop duplicates in the DataFrame
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Dropping duplicates.')
        merged_df = merged_df.drop_duplicates()
        output_df_list = []

        # Read the IDs from the file
        ID_list = Path(blastn_folder).joinpath('IDs.txt')
        IDs = [i.rstrip() for i in ID_list.open()]

        # Check if each ID is already in the DataFrame
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Constructing taxonomy table.')
        for ID in tqdm(IDs):
            if ID not in merged_df['unique ID'].values.tolist():
                # Create a new row with the ID and other relevant information
                row = [ID] + NoMatch
                output_df_list.append(row)
            else:
                row = merged_df.loc[merged_df['unique ID'] == ID].values.tolist()[0]
                output_df_list.append(row)

        ## sort table
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Saving taxonomy table...')

        merged_df.columns.tolist()
        output_df = pd.DataFrame(output_df_list, columns=merged_df.columns.tolist())
        output_df['Status'] = 'apscale blast'
        output_df.to_excel(blastn_filtered_xlsx, sheet_name='Taxonomy table', index=False)

        print('{}: Finished to filter blast results for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), blastn_folder))

if __name__ == '__main__':
    main()