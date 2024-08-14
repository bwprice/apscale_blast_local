import argparse
import os
from pathlib import Path
from Bio import SeqIO
import datetime
import subprocess
import os, shutil
import time
import glob
import multiprocessing
multiprocessing.freeze_support()
import pandas as pd
from pathlib import Path
import numpy as np
from joblib import Parallel, delayed
import threading
import pyarrow as pa
import pyarrow.parquet as pq

def fasta_subset(fasta_file, subset_size):

    print('{} : Creating subset(s) from fasta file.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    subset_size = int(subset_size)
    fasta_file = Path(fasta_file)

    ## create a new folder for the database
    subset_folder = Path(fasta_file.parent).joinpath('fasta_subsets')
    try:
        os.mkdir(subset_folder)
    except FileExistsError:
        ## delete existing subset files
        files = glob.glob(str(subset_folder / '*.fasta'))
        for f in files:
            os.remove(f)

    ## create batches from the main fasta file
    with open(fasta_file) as handle:
        i = 1
        n = 1
        chunk_fasta_files = []
        for record in SeqIO.parse(handle, 'fasta'):
            ## create a new fasta file for each chunk
            chunk_fasta = '{}/subset_{}.fasta'.format(subset_folder, i)
            ## save the name of all batches
            if chunk_fasta not in chunk_fasta_files:
                chunk_fasta_files.append(chunk_fasta)
            ## write the record to the respective file
            with open(chunk_fasta, 'a') as output_handle:
                SeqIO.write(record, output_handle, 'fasta')
            ## proceed to next chunk
            if n == subset_size:
                n = 1
                i += 1
            else:
                n += 1

    print('{} : Created {} subset(s) from fasta file.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i))

    return Path(fasta_file.parent).joinpath('fasta_subsets')

def accession2taxonomy(df_1, taxid_dict, col_names_2, db_name):
    df_2_list = []
    for row in df_1.values.tolist():
        ID_name = row[0]

        if 'db_MIDORI2' in db_name:
            # row = ['test', 'GU981053.1.<1.>964###root_1;Eukaryota_2759;Chordata_7711;Mammalia_40674;Eulipotyphla_9362;Soricidae_9376;Neomys_52813;Neomys_fodiens_62282', 100, 0.005]
            accession = row[1].split(".")[0] + '.' + row[1].split(".")[1]

        if 'db_PR2' in db_name:
            # row = ['test', "HQ446279.1.1754_U;tax=k:Eukaryota,d:TSAR,p:Alveolata-Ciliophora,c:Oligohymenophorea,o:Astomatia,f:Astomatida,g:Astomatida_X,s:Astomatida_X_sp.", 100, 0.005]
            accession = row[1].split('|')[1] + '.1'

        if 'db_SILVA' in db_name:
            # row = ['test', "HQ446279.1.1754_U;tax=k:Eukaryota,d:TSAR,p:Alveolata-Ciliophora,c:Oligohymenophorea,o:Astomatia,f:Astomatida,g:Astomatida_X,s:Astomatida_X_sp.", 100, 0.005]
            accession = row[1].split('|')[1] + '.1'

        if 'db_UNITE' in db_name:
            # row = 'Pezizaceae_sp|UDB04315424|SH0970857.10FU|reps|k__Fungi;p__Ascomycota;c__Pezizomycetes;o__Pezizales;f__Pezizaceae;g__Pezizaceae_gen_Incertae_sedis;s__Pezizaceae_sp'
            accession = row[1].split('|')[1]

        if 'db_DIATBARCODE' in db_name:
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

print_lock = threading.Lock()

def blastn_parallel(fasta_file, n_subsets, blastn_subset_folder, blastn_exe, db_folder, i, print_lock, task, max_target_seqs):
    """ Run blastn commands in parallel """

    ## create the output file
    blastn_csv = blastn_subset_folder.joinpath(Path(fasta_file).stem + '_' + task + '.csv')

    if os.path.isfile(blastn_csv):
        with print_lock:
            print('{}: Skipping {} (already exists).'.format(datetime.datetime.now().strftime('%H:%M:%S'), blastn_csv.stem))
        time.sleep(1)
    else:
        ## run blast search
        subprocess.call([blastn_exe, '-task', task, '-db', str(db_folder), '-query', str(fasta_file), '-num_threads', str(1), '-max_target_seqs', str(max_target_seqs),  '-outfmt', '6 delim=;; qseqid sseqid pident evalue', '-out', str(blastn_csv)])
        with print_lock:
            print('{}: Finished blastn for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i+1, n_subsets))

def filter_blastn_csvs(file, taxid_dict, i, n_subsets, thresholds, db_name):

    # file = '/Users/tillmacher/Desktop/APSCALE_projects/test_apscale/11_NCBI_BLAST/klimawerk_apscale_OTUs_filtered/subsets/subset_1_blastn.csv'

    ## load blast results
    col_names = ['unique ID', 'Sequence ID', 'Similarity', 'evalue']
    csv_df = pd.read_csv(file, header=None, sep=';;', names=col_names, engine='python').fillna('NAN')
    csv_df['Similarity'] = [float(i) for i in csv_df['Similarity'].values.tolist()]

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
    # ID ='OTU_99'
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
        df_2 = accession2taxonomy(df_1, taxid_dict, col_names_2, db_name)

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
                    genus = df_3['Genus'][0]
                    species = '/'.join([i.replace(genus + ' ', '') for i in df_3['Species'].values.tolist()])
                    df_3['Species'] = ['{} {}'.format(genus, species)]*len(df_3)
                    df_3['Flag'] = ['F2'] * len(df_3)
                    df_3['Ambigous taxa'] = [', '.join(df_2['Species'].drop_duplicates().values.tolist())] * len(df_3)
                    taxonomy_df = pd.concat([taxonomy_df, df_3])

                ##### 9) MULTIPLE SPECIES OF ONE GENUS? (F3)
                elif n_genera == 1 and len(df_3) != 2:
                    genus = df_3['Genus'][0]
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

def blastn_v2(blastn_exe, query_fasta, blastn_database, project_folder, n_cores, task, subset_size, max_target_seqs, apscale_gui):
    """ New version of blast that makes better use of multithreading by running multiple blast searches in parallel """

    # split fasta file
    subset_folder = fasta_subset(query_fasta, subset_size)

    ## convert to path
    project_folder = Path(project_folder)

    # run blastn command in parallel for all fasta subsets
    fasta_files = sorted(glob.glob(str(subset_folder) + '/*.fasta'))
    n_subsets = len(fasta_files)

    ## define blast task
    if task == 'Highly similar sequences (megablast)':
        task = 'megablast'
    elif task == 'More dissimilar sequences (discontiguous megablast)':
        task = 'dc-megablast'
    elif task == 'Somewhat similar sequences (blastn)':
        task = 'blastn'

    ## replace spaces and dots in the filename
    filename = Path(query_fasta).stem
    filename = filename.replace('.', '_').replace(' ', '_')

    print('{}: Starting {} for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), task, filename))
    print('{}: Your database: {}'.format(datetime.datetime.now().strftime('%H:%M:%S'),  Path(blastn_database).stem))

    ## collect files and folders
    db_folder = Path(blastn_database).joinpath('db')

    ## create a new folder for each blast search
    blastn_subset_folder = Path(project_folder).joinpath('subsets')
    if not os.path.isdir(blastn_subset_folder):
        os.mkdir(blastn_subset_folder)

    # PARALLEL BLASTN COMMAND
    Parallel(n_jobs = n_cores, backend='threading')(delayed(blastn_parallel)(fasta_file, n_subsets, blastn_subset_folder, blastn_exe, db_folder, i, print_lock, task, max_target_seqs) for i, fasta_file in enumerate(fasta_files))

    ## write the name of the database
    blastn_log = project_folder.joinpath('log.txt')
    f = open(blastn_log, 'w')
    f.write('Your database: {}\n'.format(Path(blastn_database).stem))
    f.write('Your task: {}\n'.format(task))
    f.close()

    ## write all OTUs to a table (to search OTUs without a match
    OTU_report = project_folder.joinpath('IDs.txt')
    f = open(OTU_report, 'w')
    with open(query_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            f.write(record.id + '\n')
    f.close()

    print('{}: Finished {} for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), task, filename))

    ## merge all .csv files to a single .snappy table
    # List of CSV file names
    csv_files = glob.glob('{}/*.csv'.format(str(blastn_subset_folder)))
    # Read and concatenate all CSV files into a single DataFrame
    col_names = ['unique ID', 'Sequence ID', 'Similarity', 'evalue']
    df = pd.concat((pd.read_csv(f, header=None, sep=';;', names=col_names, engine='python').fillna('NAN') for f in csv_files))
    # Convert the DataFrame to an Arrow Table
    table = pa.Table.from_pandas(df)
    # Write the Table to a Parquet file with Snappy compression
    blastn_snappy = '{}/{}.parquet.snappy'.format(project_folder, filename)
    pq.write_table(table, blastn_snappy, compression='snappy')

    ## remove subset fasta folder
    shutil.rmtree(subset_folder)

def blastn_filter(blastn_folder, blastn_database, thresholds, n_cores):
    """ Filter results according to Macher et al., 2023 (Fish Mock Community paper) """

    print('{}: Starting to filter blast results for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), blastn_folder))
    print('{}: Your database: {}'.format(datetime.datetime.now().strftime('%H:%M:%S'),  Path(blastn_database).stem))

    ## load blast results
    csv_files = glob.glob('{}/subsets/*.csv'.format(blastn_folder))

    ## collect name of database
    db_name = Path(blastn_database).stem

    ## check if a dbDNA database was used
    if "_dbDNA" in blastn_database:
        # PARALLEL FILTER COMMAND
        n_subsets = len(csv_files)
        Parallel(n_jobs = n_cores, backend='threading')(delayed(filter_blast_csvs_dbDNA)(file, i, n_subsets, thresholds) for i, file in enumerate(csv_files))

        ## also already define the no match row
        NoMatch = ["NoMatch"] * 6 + [0] * 2 + ['']*17
    else:
        ## load taxid table
        taxid_table = Path(blastn_database).joinpath('db_taxonomy.parquet.snappy')
        taxid_df = pd.read_parquet(taxid_table).fillna('')
        taxid_dict = {i[0]:i[1::] for i in taxid_df.values.tolist()}

        # PARALLEL FILTER COMMAND
        n_subsets = len(csv_files)
        # file = csv_files[0]
        # i = 1
        Parallel(n_jobs = n_cores, backend='threading')(delayed(filter_blastn_csvs)(file, taxid_dict, i, n_subsets, thresholds, db_name) for i, file in enumerate(csv_files))

        ## also already define the no match row
        NoMatch = ['NoMatch'] * 7 + [0, 1, '', '']

    # Get a list of all the xlsx files
    xlsx_files = glob.glob('{}/subsets/*.xlsx'.format(blastn_folder))

    # Create a list to hold all the individual DataFrames
    df_list = []

    # Loop through the list of xlsx files
    for file in xlsx_files:
        # Read each xlsx file into a DataFrame
        df = pd.read_excel(file).fillna('')
        # Append the DataFrame to the list
        df_list.append(df)

    # Concatenate all the DataFrames in the list into one DataFrame
    merged_df = pd.concat(df_list, ignore_index=True)
    name = Path(blastn_folder).name
    blastn_filtered_xlsx = Path('{}/{}_taxonomy.xlsx'.format(blastn_folder, name))

    ## add OTUs without hit
    # Drop duplicates in the DataFrame
    merged_df = merged_df.drop_duplicates()
    output_df_list = []

    # Read the IDs from the file
    ID_list = Path(blastn_folder).joinpath('IDs.txt')
    IDs = [i.rstrip() for i in ID_list.open()]

    # Check if each ID is already in the DataFrame
    for ID in IDs:
        if ID not in merged_df['unique ID'].values.tolist():
            # Create a new row with the ID and other relevant information
            row = [ID] + NoMatch
            output_df_list.append(row)
        else:
            row = merged_df.loc[merged_df['unique ID'] == ID].values.tolist()[0]
            output_df_list.append(row)

    ## sort table
    merged_df.columns.tolist()

    output_df = pd.DataFrame(output_df_list, columns=merged_df.columns.tolist())
    output_df['Status'] = 'apscale blast'
    output_df.to_excel(blastn_filtered_xlsx, sheet_name='Taxonomy table', index=False)

    print('{}: Finished to filter blast results for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), blastn_folder))

def main():
    """ BLASTn suite for apscale """

    message = """
    APSCALE blast command line tool - v1.0.0 - 01/07/2024
    Usage examples:
    $ apscale_blast blastn -h
    $ apscale_blast blastn -database ./MIDORI2_UNIQ_NUC_GB259_srRNA_BLAST -query_fasta ./12S_apscale_ESVs.fasta
    $ apscale_blast filter -h
    $ apscale_blast filter -database ./MIDORI2_UNIQ_NUC_GB259_srRNA_BLAST -blastn_folder ./12S_apscale_ESVs_blastn
    """

    print(message)

    # parse command line arguments
    parser = argparse.ArgumentParser(description='APSCALE blast v0.1.')

    # COMMANDS
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    blastn_parser = subparsers.add_parser('blastn', help='Perform blastn search on selected fasta file.')
    filter_parser = subparsers.add_parser('filter',  help='Perform blastn filter on selected blast result folder. Blastn must be run first!')

    # GENERAL VARIABLES
    blastn_parser.add_argument('-database', type=str, required=False, help='PATH to blastn database.')
    blastn_parser.add_argument('-apscale_gui', type=str, default=False, help='Can be ignored: Only required for APSCALE-GUI.')

    filter_parser.add_argument('-database', type=str, required=False, help='PATH to blastn database.')
    filter_parser.add_argument('-apscale_gui', type=str, default=False, help='Can be ignored: Only required for APSCALE-GUI.')

    ## VARIABLES FOR BLASTN
    blastn_parser.add_argument('-blastn_exe', type=str, default='blastn', help='PATH to blast executable. [DEFAULT: blastn]')
    blastn_parser.add_argument('-query_fasta', type=str, help='PATH to fasta file.')
    blastn_parser.add_argument('-n_cores', type=int, default=multiprocessing.cpu_count() - 1, help='Number of cores to use. [DEFAULT: CPU count - 1]')
    blastn_parser.add_argument('-task', type=str, default='blastn', help='Blastn task: blastn, megablast, or dc-megablast. [DEFAULT: blastn]')
    blastn_parser.add_argument('-out', type=str, default='./', help='PATH to output directory. A new folder will be created here. [DEFAULT: ./]')
    blastn_parser.add_argument('-subset_size', type=int, default=100, help='Number of sequences for each subset of the query fasta. [DEFAULT: 100]')
    blastn_parser.add_argument('-max_target_seqs', type=int, default=20, help='Number of hits retained from the blast search. Larger numbers will increase runtimes and required storage space [DEFAULT: 20]')

    ## VARIABLES FOR FILTER
    filter_parser.add_argument('-blastn_folder', type=str, help='PATH to blastn folder for filtering.')
    filter_parser.add_argument('-thresholds', type=str, default='97,95,90,87,85', help='Taxonomy filter thresholds. [DEFAULT: 97,95,90,87,85]')
    filter_parser.add_argument('-n_cores', type=int, default=multiprocessing.cpu_count() - 1, help='Number of cores to use. [DEFAULT: CPU count - 1]')

    args = parser.parse_args()

    ## if no args are given ask for database, query_fasta, or blastn_folder
    if args.command == 'blastn' and args.database == None and args.query_fasta == None:
        args.database = input("Please enter PATH to database: ")
        args.query_fasta = input("Please enter PATH to query fasta: ")

        args.database = args.database.strip('"')
        args.query_fasta = args.query_fasta.strip('"')

        print(args.query_fasta)

        if args.out == './':
            args.out = str(args.query_fasta).replace('.fasta', '')
            if not os.path.isdir(args.out):
                os.mkdir(Path(args.out))

    if args.command == 'filter' and args.database == None and args.blastn_folder == None:
        args.database = input("Please enter PATH to database: ")
        args.blastn_folder = input("Please enter PATH to blastn folder: ")

    ## BLASTN
    if args.command == 'blastn':
        if args.query_fasta:
            project_folder = args.out
            ## on windows the " must be removed
            args.database = args.database.strip('"')
            args.query_fasta = args.query_fasta.strip('"')
            blastn_v2(args.blastn_exe, args.query_fasta, args.database, project_folder, args.n_cores, args.task, args.subset_size, args.max_target_seqs, args.apscale_gui)
        else:
            print('\nError: Please provide a fasta file!')

    ## FILTER
    elif args.command == 'filter':
        if args.blastn_folder:
            ## on windows the " must be removed
            args.database = args.database.strip('"')
            args.blastn_folder = args.blastn_folder.strip('"')
            thresholds = args.thresholds.split(',')
            if len(thresholds) != 5:
                print('Please provide 5 comma separated values!')
                print('Using default values...')
                thresholds = '97,95,90,87,85'.split(',')
            blastn_filter(args.blastn_folder, args.database, thresholds, args.n_cores)
        else:
            print('\nError: Please provide a blast results file folder (.csv)!')

## run only if called as toplevel script
if __name__ == "__main__":
    main()

