#!/usr/bin/env python3
#This script is written for Illumina Pair-End reads
import subprocess as sp
from os import listdir, remove
from os.path import isfile, join, isdir
import re
import pandas as pd
import logging as lg
import multiprocessing
import os
from collections import defaultdict
import gzip
from functools import partial
import tempfile
import shutil
import json
import numpy as np
from tqdm import tqdm
import glob


#Number of processess for multiprocessing
n_proc = 4

lg.basicConfig(filename='log_arg_counter.log', filemode='w', format='%(asctime)s - %(levelname)s - %(message)s', level=lg.DEBUG)
lg.info('Starting ARG Counter Script')

#Test Protocol

TEST_PROTOCOL = False
TEST_SAMPLES = 2

#Quality check
step_1 = False
#Trim reads
step_2 = False
#Quality check after trimming
step_3 = False
#Create Databases
step_4 = False
#Create Resfinder DB
step_4_1 = False
#Run Diamond
step_5 = False
#Perform analys on count data
step_6 = True
# Taxonomic Profiling
step_7 = False


#Set the working directory
wdir = "/path/to/project/directory/"
out_dir = f"{wdir}data/output/"
data_dir = f"{wdir}data/dataDirectory/"
db_dir = f"{wdir}databases/"

# Holds file paths grouped by sample
samples = defaultdict(lambda: {'R1': [], 'R2': []})

# Recursively find all fastq.gz files
for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith(".fastq.gz"):
            match = re.match(r"(?P<sample>.+)_L00[567]_(?P<read>R[12])_001.fastq.gz", file)
            if match:
                sample = match.group("sample")
                read = match.group("read")
                samples[sample][read].append(os.path.join(root, file))

# Ensure consistency
for sample in samples:
    samples[sample]['R1'].sort()
    samples[sample]['R2'].sort()

lg.info(f"Discovered {len(samples)} samples")
for sample_id in sorted(samples):
    r1_count = len(samples[sample_id]['R1'])
    r2_count = len(samples[sample_id]['R2'])
    lg.info(f"  - {sample_id}: {r1_count} R1 files, {r2_count} R2 files")

def create_fasta_df(fasta_dir):
    sequences = pd.DataFrame()
    with open(fasta_dir) as f:
        ids = []
        seq = []
        sequence = ''
        for line in f:
            if line.startswith('>'):
                if len(ids) > 0:
                    seq.append(sequence)
                line = re.sub("\n", "", line)
                line = re.sub(">", "", line)
                ids.append(line)
                sequence = ''
            else:
                sequence += re.sub("\n", "", line)
        seq.append(sequence)

    sequences["IDs"] = ids
    sequences["Seq"] = seq
    return sequences


if step_1:
    lg.info(f"{5*'#'} Step 1 | Performing initial QC analysis {5*'#'}")
    os.makedirs(f"{out_dir}1_FastQC_Before", exist_ok=True)
    
    #First QC Analysis
    def fastq1_multiproc(paths):
        for fq in paths['R1'] + paths['R2']:
            try:
                lg.info(f'Performing QC: {fq}')
                sp.run(["fastqc", "-o", f"{out_dir}1_FastQC_Before/", fq])
            except Exception as e:
                lg.exception(f'Problem with QC on file: {fq}')


    with multiprocessing.Pool(processes=n_proc) as pool:
        pool.map(fastq1_multiproc, samples.values())
    try:
        lg.info(f'Performing MultiQC')
        sp.run(["multiqc", f"{out_dir}1_FastQC_Before/", "-o", f"{out_dir}1_FastQC_Before/"])
    except Exception as e:
        lg.exception(f'Problem with MultiQC')


#Quality Trimming
def trimming_multiproc(sample_tuple):
    sample, paths = sample_tuple
    try:
        out_r1 = f"{out_dir}2_Trimmed_Reads/{sample}_R1_001.fastq.gz"
        out_r2 = f"{out_dir}2_Trimmed_Reads/{sample}_R2_001.fastq.gz"

        if paths['R1'] and paths['R2']:
            sp.run([
                "bbduk.sh",
                f"in={','.join(paths['R1'])}",
                f"in2={','.join(paths['R2'])}",
                f"out={out_r1}",
                f"out2={out_r2}",
                "ref=adapters", "ktrim=r", "tbo", "qtrim=rl"
            ])
        else:
            sp.run([
                "bbduk.sh",
                f"in={','.join(paths['R1'])}",
                f"out={out_r1}",
                "ref=adapters", "ktrim=r", "tbo", "qtrim=rl"
            ])
        lg.info(f'Trimmed sample: {sample}')
    except Exception as e:
        lg.exception(f'Problem with trimming on sample: {sample} {e}')

if step_2:
    lg.info(f"{5*'#'} Step 2 | Starting quality trimming {5*'#'}")
    os.makedirs(f"{out_dir}2_Trimmed_Reads", exist_ok=True)

    sample_items = list(samples.items())
    if TEST_PROTOCOL:
        sample_items = sample_items[:TEST_SAMPLES]
        lg.info(f'TEST_PROTOCOL enabled: Trimming only first {TEST_SAMPLES} samples')

    for sample, paths in tqdm(sample_items, desc="Trimming samples", unit="sample"):
        try:
            out_r1 = f"{out_dir}2_Trimmed_Reads/{sample}_R1_001.fastq.gz"
            out_r2 = f"{out_dir}2_Trimmed_Reads/{sample}_R2_001.fastq.gz"

            r1_files = paths['R1']
            r2_files = paths['R2']

            with tempfile.TemporaryDirectory() as tempdir:
                # Concatenate R1 files
                r1_cat = os.path.join(tempdir, f"{sample}_R1_cat.fastq.gz")
                with open(r1_cat, 'wb') as wout:
                    for fq in r1_files:
                        with open(fq, 'rb') as fin:
                            shutil.copyfileobj(fin, wout)

                cmd = [
                    "bbduk.sh",
                    f"in={r1_cat}",
                    f"out={out_r1}",
                    "ref=adapters", "ktrim=r", "tbo", "qtrim=rl",
                    "threads=30"  # Set your desired threads here
                ]

                if r2_files:
                    r2_cat = os.path.join(tempdir, f"{sample}_R2_cat.fastq.gz")
                    with open(r2_cat, 'wb') as wout:
                        for fq in r2_files:
                            with open(fq, 'rb') as fin:
                                shutil.copyfileobj(fin, wout)
                    cmd += [f"in2={r2_cat}", f"out2={out_r2}"]

                lg.info(f"Running bbduk for sample: {sample}")
                sp.run(cmd, check=True)
                lg.info(f"Trimmed sample: {sample}")

        except Exception as e:
            lg.exception(f'Problem with trimming on sample: {sample} {e}')

#Second QC Analysis
if step_3:
    lg.info(f"{5*'#'} Step 3 | Performing post-trimming QC analysis {5*'#'}")
    os.makedirs(f"{out_dir}3_FastQC_After", exist_ok=True)

    # List all fastq.gz files in the trimmed reads directory
    all_trimmed_files = [
        f for f in listdir(f"{out_dir}2_Trimmed_Reads")
        if isfile(join(f"{out_dir}2_Trimmed_Reads", f)) and f.endswith(".fastq.gz")
    ]

    # If TEST_PROTOCOL is True, only use files from first 2 samples
    if TEST_PROTOCOL:
        sample_names = list(samples.keys())[:TEST_SAMPLES]
        all_trimmed_files = [
            f for f in all_trimmed_files
            if any(f.startswith(s) for s in sample_names)
        ]
        lg.info(f'TEST_PROTOCOL enabled: QC only on trimmed files from first 2 samples')

    def fastq2_multiproc(filename):
        fq_path = join(f"{out_dir}2_Trimmed_Reads", filename)
        try:
            lg.info(f'Performing Post-Trimming QC: {fq_path}')
            sp.run(["fastqc", "-o", f"{out_dir}3_FastQC_After/", fq_path], check=True)
        except Exception as e:
            lg.exception(f'Problem with Post-Trimming QC on file: {fq_path}')

    with multiprocessing.Pool(processes=n_proc) as pool:
        pool.map(fastq2_multiproc, all_trimmed_files)

    try:
        lg.info('Running MultiQC on post-trimming FastQC reports')
        sp.run(["multiqc", f"{out_dir}3_FastQC_After/", "-o", f"{out_dir}3_FastQC_After/"], check=True)
    except Exception as e:
        lg.exception(f'Problem with Post-Trimming MultiQC')


#Create Diamond Databases
if step_4:
    lg.info(f"{5*'#'} Step 4 | Creating Diamond Databases {5*'#'}")
    
    if TEST_PROTOCOL:
        lg.info(f'TEST_PROTOCOL enabled: Running Step 4 in test mode')

    def db_make(dir_name, input_file, format='Nucleotide'):
        # Validate input files exist
        if format == 'Nucleotide':
            if not os.path.exists(input_file):
                lg.error(f"Input file not found: {input_file}")
                return False
            
            #Cluster sequences to 100% similarity to remove duplicates
            lg.info('Clustering Sequences')
            try:
                result = sp.run(["vsearch", "--cluster_fast", input_file,
                                "--centroids", f"{out_dir}4_Diamond_Databases/{dir_name}/2_clstr_seqs.fasta",
                                "--id", "1.0", "--fasta_width", "0", "--notrunclabels"], check=True)
                if result.returncode != 0:
                    raise Exception("VSEARCH clustering failed")
            except Exception as e:
                lg.exception(f'Error in clustering: {e}')
                return False

            #Translate nucleotide reference sequence to protein sequence
            lg.info('Translating Nucleotide Sequences')
            try:
                result = sp.run(["transeq", "-sequence", f"{out_dir}4_Diamond_Databases/{dir_name}/2_clstr_seqs.fasta",
                               "-outseq", f"{out_dir}4_Diamond_Databases/{dir_name}/3_reference_sequences_prot.fasta",
                               "-frame", "6"], check=True)
                if result.returncode != 0:
                    raise Exception("Translation failed")
            except Exception as e:
                lg.exception(f'Error in translation: {e}')
                return False

        #Import protein fasta
        prot_file = f"{out_dir}4_Diamond_Databases/{dir_name}/3_reference_sequences_prot.fasta"
        if not os.path.exists(prot_file):
            lg.error(f"Protein file not found: {prot_file}")
            return False

        try:
            lg.info(f'Importing Sequences from {prot_file}')
            sequences = create_fasta_df(prot_file)
            if sequences.empty:
                raise Exception("No sequences found in FASTA file")
        except Exception as e:
            lg.exception(f'Could not import sequences: {e}')
            return False

        lg.info('Filtering bad ORFs')
        #Remove ORFs with more than 1 stop codon (with proper parentheses for operator precedence)
        good_sequences = sequences[
            (sequences['Seq'].str.count('\*') <= 1) & 
            (sequences['Seq'].str.endswith('*') | ~sequences['Seq'].str.contains('\*'))
        ]
        bad_sequences = sequences[~sequences.index.isin(good_sequences.index)]
        
        lg.info(f"Filtering results - Total: {len(sequences)}, Good: {len(good_sequences)}, Bad: {len(bad_sequences)}")

        try:
            # Write the sequence data to a FASTA file
            lg.info('Writing Filtered ORF FASTA File')
            output_file = f'{out_dir}4_Diamond_Databases/{dir_name}/4_filtered_ORFs.fasta'
            with open(output_file, 'w') as f:
                for i, row in good_sequences.iterrows():
                    f.write(f">{row['IDs']}\n")
                    f.write(f"{row['Seq']}\n")
        except Exception as e:
            lg.exception(f'Could not write new FASTA file: {e}')
            return False

        try:
            #Make Diamond DB from Protein reference file.
            lg.info('Making Diamond DB')
            result = sp.run(["diamond", "makedb", "--in", output_file,
                           "-d", f"{out_dir}4_Diamond_Databases/{dir_name}/5_diamond_db"], check=True)
            if result.returncode != 0:
                raise Exception("DIAMOND makedb failed")
            return True
        except Exception as e:
            lg.exception(f'Error while creating database: {e}')
            return False

    # Create main output directory
    os.makedirs(f"{out_dir}4_Diamond_Databases", exist_ok=True)

    if step_4_1 and (not TEST_PROTOCOL or TEST_SAMPLES >= 2):
        lg.info('Step 4.4 | Creating ResFinder Database')
        os.makedirs(f"{out_dir}4_Diamond_Databases/4_4_ResFinder_DB", exist_ok=True)
        
        # Verify input file exists (assuming it's already prepared)
        input_file = f"{db_dir}resfinder_db/1_reference_sequences_nuc.fasta"
        if not os.path.exists(input_file):
            lg.error(f"Input file not found: {input_file}")
            exit()
        
        if db_make('4_4_ResFinder_DB', input_file, 'Nucleotide'):
            lg.info('Step 4.4 | Completed')
        else:
            lg.error('Step 4.4 | Failed')

    lg.info(f"{5*'#'} Step 4 | Completed {5*'#'}")

#Run Diamond Blast
if step_5:
    lg.info(f"{5*'#'} Step 5 | Starting Diamond BLAST {5*'#'}")
    
    # Create output directory
    os.makedirs(f"{out_dir}5_Diamond_Results", exist_ok=True)
    
    # Get list of database directories
    db_dirs = [d for d in os.listdir(f"{out_dir}4_Diamond_Databases/") 
            if os.path.isdir(f"{out_dir}4_Diamond_Databases/{d}")]

    # Only keep ResFinder (step 4.4)
    db_dirs = [d for d in db_dirs if d == "4_4_ResFinder_DB"]
    
    # Get list of input files (handles both compressed and uncompressed FASTQ)
    raw_files = sorted([f for f in os.listdir(f"{out_dir}2_Trimmed_Reads") 
                if f.endswith((".fastq.gz", ".fastq"))])
    
    if TEST_PROTOCOL:
        raw_files = raw_files[:TEST_SAMPLES*2]  # Account for R1 and R2
    
    def diamond_multiproc(args):
        file, db = args
        try:
            sample_name = re.sub(r"\.fastq(\.gz)?$", "", file)
            output_file = f"{sample_name}_{db[4:]}_diamond_results.tsv"
            output_path = f"{out_dir}5_Diamond_Results/5{db[1:]}/{output_file}"
            db_path = f"{out_dir}4_Diamond_Databases/{db}/5_diamond_db.dmnd"
            input_path = f"{out_dir}2_Trimmed_Reads/{file}"
            
            # Validate inputs exist
            if not os.path.exists(db_path):
                lg.error(f"Diamond database not found: {db_path}")
                return
            if not os.path.exists(input_path):
                lg.error(f"Input file not found: {input_path}")
                return
            
            lg.info(f'Running Diamond BLASTx for {sample_name} vs {db[4:]}')
            
           
            # Run Diamond with optimized parameters
            result = sp.run([
                "diamond", "blastx",
                "-d", db_path,
                "-q", input_path,
                "-o", output_path,
                "-f", "6", "qseqid", "sseqid", "pident", "length", "mismatch", 
                         "gapopen", "qstart", "qend", "sstart", "send", 
                         "evalue", "bitscore",  # Extended output format
                "--id", "70",                  # Minimum identity percentage
                "--query-cover", "80",          # Minimum query coverage
                "--max-target-seqs", "1",       # Top hit only
                "--evalue", "1e-5",             # Strict e-value threshold
                "--more-sensitive",             # Increased sensitivity
                "--threads", "10",
                "--block-size", "8",            # Memory usage (GB)
                "--tmpdir", "/tmp",             # Explicit temp directory
            ], check=True, stderr=sp.PIPE)
            
            if result.returncode == 0:
                lg.info(f'Successfully completed {sample_name} vs {db[4:]}')
            else:
                lg.error(f'Failed: {sample_name} vs {db[4:]}')
                lg.error(result.stderr.decode())
                
        except sp.CalledProcessError as e:
            lg.exception(f'Diamond failed for {file}: {e.stderr.decode()}')
        except Exception as e:
            lg.exception(f'Unexpected error with {file}: {str(e)}')
    
    for db in db_dirs:
        db_output_dir = f"{out_dir}5_Diamond_Results/5{db[1:]}"
        os.makedirs(db_output_dir, exist_ok=True)
        
        lg.info(f'Processing database: {db[4:]}')
        
        # Prepare arguments for multiprocessing
        args_list = [(file, db) for file in raw_files]
        
        # Run with progress tracking
        with multiprocessing.Pool(processes=4) as pool:
            try:
                for _ in tqdm(pool.imap_unordered(diamond_multiproc, args_list),
                            total=len(args_list),
                            desc=f"BLASTing against {db[4:]}"):
                    pass
            except Exception as e:
                lg.exception(f"Multiprocessing error: {str(e)}")
                pool.terminate()
                raise
            finally:
                pool.close()
                pool.join()
        
        lg.info(f'Completed all samples for database {db[4:]}')
    
    lg.info(f"{5*'#'} Step 5 | Completed Diamond BLAST {5*'#'}")

#Perform some data screening, filtering and adding phenotype information. Counting the number of ARG occurances and grouping by various characteristics.
if step_6:
    lg.info(f"{5*'#'} Step 6 | Counting Resistance Gene Occurrences {5*'#'}")
    os.makedirs(f"{out_dir}6_ARG_Counts", exist_ok=True)
    
    # Get list of database result directories
    db_dirs = [d for d in os.listdir(f"{out_dir}5_Diamond_Results/") 
              if os.path.isdir(f"{out_dir}5_Diamond_Results/{d}")]
    
    # Only keep ResFinder (step 4.4)
    db_dirs = [d for d in db_dirs if d == "5_4_ResFinder_DB"]
    
    for db in db_dirs:
        db_name = db[4:]  # Remove '5_1_' prefix to get database name
        db_output_dir = f"{out_dir}6_ARG_Counts/6{db[1:]}"
        os.makedirs(db_output_dir, exist_ok=True)
        
        lg.info(f'Processing database: {db_name}')
        

        if db_name == "ResFinder_DB":
            phenotypes = pd.read_csv(
                f"{db_dir}resfinder_db/phenotypes.txt",
                delimiter='\t',
                skiprows=1,
                names=["ID", "Class", "Phenotype", "PMID", "Mechanism", "Notes", "Required_Genes"],
                dtype={'ID': str, 'Class': str, 'Phenotype': str}
            )
            merge_key = ('ID', 'ID')


        # Get all diamond result files for this database
        diamond_files = [f for f in os.listdir(f"{out_dir}5_Diamond_Results/{db}") 
                        if f.endswith(".tsv") and isfile(join(f"{out_dir}5_Diamond_Results/{db}", f))]

        # Extract unique sample names from the filenames
        sample_names = set()
        for f in diamond_files:
            match = re.match(r"(?P<sample>.+)_(R[12])_001_(.+)_diamond_results\.tsv", f)
            if match:
                sample_names.add(match.group("sample"))
        
        if TEST_PROTOCOL:
            sample_names = list(sample_names)[:TEST_SAMPLES]
            lg.info(f'TEST_PROTOCOL enabled: Processing only first {TEST_SAMPLES} samples')
        
        def process_sample(sample):
            lg.info(f'Processing sample: {sample}')
            
            # Check for both R1 and R2 files
            r1_file = f"{out_dir}5_Diamond_Results/{db}/{sample}_R1_001_{db_name}_diamond_results.tsv"
            r2_file = f"{out_dir}5_Diamond_Results/{db}/{sample}_R2_001_{db_name}_diamond_results.tsv"
            
            # Read and combine R1 and R2 files if they exist
            dfs = []
            if os.path.exists(r1_file):
                df_r1 = pd.read_csv(r1_file, sep='\t', 
                                  names=["Read", "ID", "pident", "length", "mismatch", "gapopen", 
                                         "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
                dfs.append(df_r1)
            if os.path.exists(r2_file):
                df_r2 = pd.read_csv(r2_file, sep='\t', 
                                  names=["Read", "ID", "pident", "length", "mismatch", "gapopen", 
                                         "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
                dfs.append(df_r2)
            
            if not dfs:
                lg.warning(f"No Diamond results found for sample {sample}")
                return None
            
            df = pd.concat(dfs)
            
            # Filter results
            mask = df['pident'] >= 90  # 90% identity threshold
            df_filtered = df[mask].sort_values(['bitscore', 'Read'], ascending=[False, False])
            
            # Remove duplicate reads (keep best hit per read)
            df_filtered = df_filtered.drop_duplicates(subset=["Read"], keep='first')
            
            # Save filtered results
            df_filtered.to_csv(f"{db_output_dir}/{sample}_{db_name}_filtered_counts.tsv", sep="\t", index=False)
            
            # MERGE WITH PHENOTYPE DATA
            try:
                # Get the correct merge columns for this database
                df_col, pheno_col = merge_key
                
                # Merge with phenotype information
                df_with_pheno = pd.merge(
                    df_filtered,
                    phenotypes,
                    left_on=df_col,
                    right_on=pheno_col,
                    how='left'
                )
                
                # Save phenotype-annotated results
                df_with_pheno.to_csv(f"{db_output_dir}/{sample}_{db_name}_phenotype_counts.tsv", 
                                   sep="\t", index=False)
                
                count_column = None
                if db_name == "ResFinder_DB":
                    count_column = "ID"
                
                if count_column and count_column in df_with_pheno.columns:
                    # Split combined categories (e.g., "CompoundA; CompoundB")
                    df_with_pheno = df_with_pheno.assign(
                        Class=df_with_pheno[count_column].str.split(r'[;,]\s*')
                    ).explode('Class')
                    
                    # Clean up class/compound names
                    df_with_pheno['Class'] = df_with_pheno['Class'].str.strip()
                    df_with_pheno = df_with_pheno[df_with_pheno['Class'] != '']
                    
                    # Count occurrences by class/compound
                    class_counts = df_with_pheno['Class'].value_counts().reset_index()
                    class_counts.columns = ['Class', 'Count']
                    class_counts['Sample'] = sample
                    
                    return class_counts
                
            except Exception as e:
                lg.exception(f"Error merging phenotype data for sample {sample}: {str(e)}")
                # Debug: Print sample IDs that failed to merge
                if db_name == "BacMet_DB":
                    lg.error(f"Sample IDs: {df_filtered['ID'].unique()[:5]}")
                    lg.error(f"Phenotype IDs: {phenotypes['GI_number'].unique()[:5]}")
            
            return None
        
        # Process all samples for this database
        all_class_counts = []
        with multiprocessing.Pool(processes=n_proc) as pool:
            results = pool.map(process_sample, sample_names)
            all_class_counts = [res for res in results if res is not None]
        
        # Combine all class/compound counts and save
        if all_class_counts:
            combined_counts = pd.concat(all_class_counts)
            output_file = f"{db_output_dir}/{db_name}_counts_by_group.tsv"
            combined_counts.to_csv(output_file, sep="\t", index=False)
            lg.info(f"Saved counts to {output_file}")
        
        lg.info(f'Completed processing for database: {db_name}')
    
    lg.info(f"{5*'#'} Step 6 | Completed Resistance Gene Counting {5*'#'}")

if step_7: #Database GTDB_r226 accessed 28/07/2025
    lg.info(f"{5*'#'} Step 8 | Metagenomic Profiling {5*'#'}")
    os.makedirs(f"{out_dir}8_Taxonomy", exist_ok=True)

    taxonomy_dir = f"{out_dir}8_Taxonomy/taxonomy_file_folder"
    if not os.path.exists(taxonomy_dir):
        os.makedirs(taxonomy_dir, exist_ok=True)
        sp.run(["sylph-tax", "download", "--download-to", f"{taxonomy_dir}"])

    r1_files = glob.glob(f"{out_dir}2_Trimmed_Reads/*_R1_001.fastq.gz")
    r2_files = glob.glob(f"{out_dir}2_Trimmed_Reads/*_R2_001.fastq.gz")

    # Ensure the files are sorted to match correctly (optional but often needed)
    r1_files.sort()
    r2_files.sort()

    sp.run([
        "sylph", "profile",
        f"{db_dir}sylph_db/gtdb-r226-c200-dbv1.syldb",
        "-1", *r1_files,
        "-2", *r2_files,
        "-t", "40",
        "-o", f"{out_dir}8_Taxonomy/results.tsv"
    ])
    sp.run(["sylph-tax", "taxprof", f"{out_dir}8_Taxonomy/results.tsv", "-t", "GTDB_r226", "-o", f"{out_dir}/8_Taxonomy/taxonomy_"])
    sylphmpa_files = glob.glob(f"{out_dir}8_Taxonomy/*.sylphmpa")
    sp.run(["sylph-tax", "merge", *sylphmpa_files, "--column", "relative_abundance", "-o", f"{out_dir}8_Taxonomy/merged_taxonomy_abundance.tsv"])
