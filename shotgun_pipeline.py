#!/usr/bin/env python

import numpy as np
import pandas as pd
from collections import defaultdict
import subprocess
import os
from tqdm import tqdm
from loguru import logger
import csv
import argparse
import sys

__version__ = "2026.01.30"

def get_sample(sample_id, sra_path='~/bin/sratoolkit.3.2.0-centos_linux64/bin', skip_if_exists=True, paired=False, log_file='process.log'):
    '''Download a single sample from SRA given its SRA ID

    Parameters
    ----------
    sample_id: str
            SRA sample ID (SRRxxxxxx)
    sra_path: str, optional
            path to the sra-toolkit binaries
    skip_if_exists: bool, optional
            if true, skip downloading if the sample already exists
    paired: bool, optional
            whether the input data is paired-end (True) or single-end (False)
    '''
    sra_path = os.path.expanduser(sra_path)
    logger.info(f"Downloading sample {sample_id}")
    if skip_if_exists:
        if os.path.exists(f"{sample_id}_1.fastq"):
                logger.info(f"Sample {sample_id} already exists, skipping download")
                return
    # prefetch
    params = [os.path.join(sra_path, 'prefetch'), sample_id]
    logger.debug(f"Running command: {' '.join(params)}")
    subprocess.call(params)
    logger.debug(f"Prefetched sample {sample_id}")
    # fasterq-dump
    params = [os.path.join(sra_path, 'fasterq-dump'), './' + sample_id, '--threads', '4']
    if paired:
        # if paired, we create a single file with both f and r reads (we don't care about separating them since we will be aligning to uniref and not doing assembly)
        params.append('--split-spot')
    else:
        # not pair - we will use only the first read, so we will ignore the _2 file if it is created
        params.remove('--split-files')
    logger.debug(f"Running command: {' '.join(params)}")
    with open(log_file, 'a') as logfile:
        res = subprocess.call(params, stdout=logfile, stderr=logfile)
    if res != 0:
        logger.error(f"fasterq-dump failed with return code {res}")
        raise RuntimeError("fasterq-dump execution failed")

    # check is split to only one file, rename it to standard name
    if not os.path.exists(f"{sample_id}_1.fastq"):
        if os.path.exists(f"{sample_id}.fastq"):
                os.rename(f"{sample_id}.fastq", f"{sample_id}_1.fastq")
                logger.info(f'Renamed {sample_id}.fastq to {sample_id}_1.fastq')

    logger.info(f"Converted sample {sample_id} to fastq")
    return


def clean_sample(sample_id, fastp_path='~/bin/fastp', skip_if_exists=True, log_file='process.log'):
    '''Clean a single sample using fastp

    Parameters
    ----------
    sample_id: str
            SRA sample ID (SRRxxxxxx)
    fastp_path: str, optional
            path to the fastp binary
    skip_if_exists: bool, optional
            if true, skip cleaning if the cleaned sample already exists
    '''
    logger.info(f"Cleaning sample {sample_id}")
    fastp_path = os.path.expanduser(fastp_path)
    if skip_if_exists:
        if os.path.exists(f"{sample_id}-1.clean.fastq"):
                logger.info(f"Cleaned sample {sample_id} already exists, skipping cleaning")
                return
    # run fastp
    input_r1 = f"{sample_id}_1.fastq"
    output_r1 = f"{sample_id}-1.clean.fastq"
    params = [fastp_path, '-i', input_r1, '-o', output_r1, '--length_required', '50', '--qualified_quality_phred', '20', '--unqualified_percent_limit', '1']
    logger.debug(f"Running command: {' '.join(params)}")
    with open(log_file, 'a') as logfile:
        res = subprocess.call(params, stdout=logfile, stderr=logfile)
    if res != 0:
        logger.error(f"fastp failed with return code {res}")
        raise RuntimeError("fastp execution failed")
    logger.debug(f"Cleaned sample {sample_id}")
    return


def convert_to_fasta(sample_id, seqtk_path='seqtk', skip_if_exists=True, log_file='process.log'):
    '''Convert cleaned fastq files to fasta using seqtk

    Parameters
    ----------
    sample_id: str
            SRA sample ID (SRRxxxxxx)
    seqtk_path: str, optional
            path to the seqtk binary
    skip_if_exists: bool, optional
            if true, skip conversion if the fasta file already exists
    '''
    seqtk_path = os.path.expanduser(seqtk_path)
    logger.info(f"Converting sample {sample_id} to fasta")
    if skip_if_exists:
        if os.path.exists(f"{sample_id}-1.clean.fasta"):
                logger.info(f"Fasta file for sample {sample_id} already exists, skipping conversion")
                return
    input_fastq = f"{sample_id}-1.clean.fastq"
    output_fasta = f"{sample_id}-1.clean.fasta"
    params = [seqtk_path, 'seq', '-a', input_fastq]
    logger.debug(f"Running command: {' '.join(params)}")
    with open(output_fasta, 'w') as outfile, open(log_file, 'a') as logfile:
        res = subprocess.call(params, stdout=outfile, stderr=logfile)
    if res != 0:
        logger.error(f"seqtk failed with return code {res}")
        raise RuntimeError("seqtk execution failed")
    logger.debug(f"Converted sample {sample_id} to fasta")
    return

def rarify_fasta(sample_id, depth, seqtk_path='seqtk', skip_if_exists=True, random_seed=2026, log_file='process.log'):
    '''Rarify a fasta file to a specified depth by randomly subsampling readsm using seqtk
    if depth is None, no rarification is performed and the original fasta file is copied to the rarified fasta file (if it doesn't already exist)
    
    Parameters
    ----------
    sample_id: str
        SRA sample ID (SRRxxxxxx)
    depth: int
        number of reads to rarify to
    seqtk_path: str, optional
        path to the seqtk binary
    skip_if_exists: bool, optional
        if true, skip rarification if the rarified fasta file already exists
    log_file: str, optional
        log file path
    '''
    if depth is None:
        logger.info(f"No rarification depth specified, copying original fasta file for sample {sample_id} to rarified fasta file")
        input_fasta = f"{sample_id}-1.clean.fasta"
        output_fasta = f"{sample_id}-1.clean.rarified.fasta"
        if skip_if_exists:
            if os.path.exists(output_fasta):
                    logger.info(f"Rarified fasta file for sample {sample_id} already exists, skipping copy")
                    return
        subprocess.call(['cp', input_fasta, output_fasta])
        logger.debug(f"Copied original fasta file for sample {sample_id} to rarified fasta file")
        return
    logger.info(f"Rarifying sample {sample_id} to depth {depth}")
    seqtk_path = os.path.expanduser(seqtk_path)
    input_fasta = f"{sample_id}-1.clean.fasta"
    output_fasta = f"{sample_id}-1.clean.rarified.fasta"
    if skip_if_exists:
        if os.path.exists(output_fasta):
                logger.info(f"Rarified fasta file for sample {sample_id} already exists, skipping rarification")
                return
    # run seqtk sample to rarify the fasta file
    params = [seqtk_path, 'sample', '-s', str(random_seed), input_fasta, str(depth)]
    logger.debug(f"Running command: {' '.join(params)}")
    with open(output_fasta, 'w') as outfile, open(log_file, 'a') as logfile:
        res = subprocess.call(params, stdout=outfile, stderr=logfile)
    if res != 0:
        logger.error(f"seqtk sample failed with return code {res}")
        raise RuntimeError("seqtk sample execution failed")
    logger.debug(f"Rarified sample {sample_id} to depth {depth}")
    return

def align_to_uniref(sample_id, diamond_db='~/databases/uniref/db-uniref90.dmnd', diamond_path='~/bin/diamond', skip_if_exists=True, log_file='process.log',sensitivity='fast', threads='10', iterate=False, tmp_dir=None):
    '''Align input fasta file to UniRef database using DIAMOND
    Parameters
    ----------
    input_fasta: str
        path to input fasta file
    output_file: str
        path to output file
    diamond_db: str, optional
        path to DIAMOND UniRef database
    diamond_path: str, optional
        path to DIAMOND binary
    skip_if_exists: bool, optional
        if true, skip alignment if the output file already exists
    sensitivity: str, optional
        sensitivity mode for DIAMOND (fast, sensitive, more-sensitive)
    iterate: bool, optional
        whether to use the --iterate flag for DIAMOND alignment (for iterative searches, recommended for more sensitive modes)
    threads: str, optional
        number of threads to use for DIAMOND alignment
    tmp_dir: str, optional
        path to temporary directory to use for DIAMOND (if not set, will use current directory)
    '''
    output_file = f"{sample_id}-aligned.txt"
    input_fasta = f"{sample_id}-1.clean.rarified.fasta"
    if skip_if_exists:
        if os.path.exists(output_file):
                logger.info(f"Alignment output file {output_file} already exists, skipping alignment")
                return
    logger.info(f'Aligning {input_fasta} to UniRef database using DIAMOND')
    diamond_path = os.path.expanduser(diamond_path)
    diamond_db = os.path.expanduser(diamond_db)
    diamond_output_columns = ['sseqid','qseqid', 'qstart', 'qend', 'qframe', 'qstrand', 'sstart', 'send', 'qseq', 'qseq_translated', 'length', 'evalue', 'bitscore', 'gapopen', 'cigar', 'full_qseq']
    command = [
        diamond_path,
        'blastx',
        '--db', diamond_db,
        '--out', output_file,
        '--' + sensitivity,
        '--min-orf', '1',
        '--query', input_fasta,
        '--strand', 'both',
        '--max-target-seqs', '1',
        '--un', output_file + '.unmatched.fasta',
        '--threads', threads]
    # if using tmp dir, add it to the command
    if tmp_dir is not None:
        tmp_dir = os.path.expanduser(tmp_dir)
        tmp_dir = os.path.join(tmp_dir, sample_id)
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        command.extend(['--tmpdir', tmp_dir])
    # add --iterate flag for more sensitive modes to improve performance
    if iterate:
        command.append('--iterate')
    # set the output format to tabular with specific columns
    command.extend(['--outfmt', '6'])
    command.extend(diamond_output_columns)
    with open(log_file, 'a') as logfile:
        res = subprocess.run(command, stdout=logfile, stderr=logfile)
    if res.returncode != 0:
        logger.error(f"DIAMOND failed with return code {res.returncode}")
        raise RuntimeError("DIAMOND execution failed")
    logger.info("DIAMOND completed successfully")
    return


def split_to_uniref(sample_id, skip_if_exists=True, min_keep=50, log_file='process.log', buffer_threshold=1000):
    '''Split DIAMOND output file into separate files per UniRef ID

    Parameters
    ----------
    sample_id: str
        SRA sample ID (SRRxxxxxx)
    skip_if_exists: bool, optional
        if true, skip splitting if the output directory already exists and is not empty
    min_keep: int, optional
        minimum number of reads to keep a uniref file
    log_file: str, optional
        log file path
    buffer_threshold: int, optional
        number of reads to buffer before flushing to file
    '''
    diamond_output = f"{sample_id}-aligned.txt"
    out_dir = f"{sample_id}-splits"
    logger.info(f"Splitting DIAMOND output {diamond_output} into directory {out_dir}")
    if skip_if_exists:
        if os.path.exists(out_dir) and len(os.listdir(out_dir)) > 0:
                logger.info(f"Output directory {out_dir} already exists and is not empty, skipping splitting")
                return
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # Dictionary to store buffers and counts for each uniref_id
    buffers = defaultdict(list)
    uniref_counts = defaultdict(int)
    
    def flush_buffer(uniref_id, force=False):
        """Flush buffer for a specific uniref_id to file"""
        if uniref_id in buffers and len(buffers[uniref_id]) > 0:
            # During final flush, only write if we have enough reads
            if force and uniref_counts[uniref_id] < min_keep:
                buffers[uniref_id].clear()
                return
            
            with open(os.path.join(out_dir, f"{uniref_id}.txt"), 'a') as out_f:
                out_f.writelines(buffers[uniref_id])
            buffers[uniref_id].clear()
    
    with open(diamond_output, 'r') as f:
        for line in f:
            parts = line.strip().split('\t', maxsplit=1)
            uniref_id = parts[0]
            uniref_counts[uniref_id] += 1
            buffers[uniref_id].append(line)
            
            # Flush if buffer reaches threshold
            if len(buffers[uniref_id]) >= buffer_threshold:
                flush_buffer(uniref_id)
    
    # Final flush of all remaining buffers
    kept_uniref_ids = 0
    for uniref_id in list(buffers.keys()):
        if uniref_counts[uniref_id] >= min_keep:
            flush_buffer(uniref_id, force=True)
            kept_uniref_ids += 1
        else:
            buffers[uniref_id].clear()  # Discard buffers below threshold
    
    logger.info(f"Completed splitting DIAMOND output into {kept_uniref_ids} uniref ids (kept), total lines: {sum(count for uid, count in uniref_counts.items() if count >= min_keep)}")
    return


def sample_pipeline(sample_id, skip_if_exists=True, start_step=0, database='~/databases/uniref/db-uniref50.dmnd', sensitivity='fast', threads='10', iterate=False, paired=False, depth=None, tmp_dir=None):
    '''Process a single sample given its SRA ID
    Steps:
    1. Download the sample using sra-toolkit prefetch+fasterq-dump
    2. Quality control using fastp
    3. convert to fasta using seqtk
    4. Align to UniRef using diamond
    5. split to per-unirefID files

    Parameters
    ----------
    sample_id: str
            SRA sample ID (SRRxxxxxx)
    skip_if_exists: bool, optional
            if true, skip each processing step if the relevant output file already exists
    start_step: int, optional
            step to start from (0: download, 1: clean, 2: convert to fasta, 3: align, 4: split)
    database: str, optional
        location of the diamond uniref database to use for alignment
    sensitivity: str, optional
        sensitivity mode for DIAMOND (fast, sensitive, more-sensitive)
    paired: bool, optional
        whether the input data is paired-end (True) or single-end (False)
    depth: int, optional
        if set, rarify each sample to this depth
    threads: str, optional
        number of threads to use for diamond alignment
    iterate: bool, optional
        whether to use the --iterate flag for diamond alignment (for iterative searches, recommended for more sensitive modes)
    tmp_dir: str, optional
        path to temporary directory to use for DIAMOND (if not set, will use current directory)
    '''
    log_file = f'process-{sample_id}.log'
    logger.info(f"Processing sample {sample_id}")
    if start_step <= 0:
        # Step 0: Download the sample
        get_sample(sample_id, skip_if_exists=skip_if_exists, paired=paired, log_file=log_file)
    if start_step <= 1:
        # Step 1: Clean the sample
        clean_sample(sample_id, skip_if_exists=skip_if_exists, log_file=log_file)
    if start_step <= 2:
        # Step 2: Convert to fasta
        convert_to_fasta(sample_id, skip_if_exists=skip_if_exists, log_file=log_file)
    if start_step <= 3:
        # Step 3: Rarify to specified depth (if depth is None, no rarification is performed and the original fasta file is copied to the rarified fasta file)
        rarify_fasta(sample_id, depth=depth, skip_if_exists=skip_if_exists, log_file=log_file)
    if start_step <= 4:
        # Step 3: Align to UniRef
        align_to_uniref(sample_id, skip_if_exists=skip_if_exists, log_file=log_file, diamond_db=database, sensitivity=sensitivity, threads=threads, iterate=iterate, tmp_dir=tmp_dir)
    if start_step <= 5:
        # Step 4: Split to per-UniRef ID files
        split_to_uniref(sample_id, skip_if_exists=skip_if_exists, log_file=log_file)
    logger.info(f"Finished processing sample {sample_id}")
    return


def main(argv):
    parser = argparse.ArgumentParser(description='Shotgun pipeline version ' + __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--accession', help='SRA sample accession to process', required=True)
    parser.add_argument('--skip-if-exists', action='store_true', help='Skip processing steps if output files already exist', default=True)
    parser.add_argument('--start-step', type=int, help='Step to start from (0: download, 1: clean, 2: convert to fasta, 3: rarify, 4: align, 5: split)', default=0)
    parser.add_argument('--database', type=str, help='Path to the database to use for alignment', default='~/databases/uniref/db-uniref50.dmnd')
    parser.add_argument('--sensitivity', type=str, help='Sensitivity mode for DIAMOND (fast, sensitive, more-sensitive)', default='fast')
    parser.add_argument('--threads', type=str, help='Number of threads to use for diamond alignment', default='10')
    parser.add_argument('--type', type=str, help='if "uniref50" or "uniref90" use relevant defaults (database, sensitivity, iterate)', default=None)
    parser.add_argument('--iterate', action='store_true', help='diamond --iterate flag (for iterative searches)', default=False)
    parser.add_argument('--paired', action='store_true', help='Whether to use paired-end data (True) or only forward reads (False)', default=False)
    parser.add_argument('--depth', type=int, help='If set, rarify each sample to this depth', default=None)
    parser.add_argument('--tmp-dir', type=str, help='Path to temporary directory to use for DIAMOND (if not set, will use current directory)', default=None)
    args = parser.parse_args(sys.argv[1:])
    # add file logging
    logger.add("shotgun_pipeline.log", rotation="10 MB")
    logger.info("Starting shotgun pipeline")
    sample_pipeline(args.accession, skip_if_exists=args.skip_if_exists, start_step=args.start_step, database=args.database, sensitivity=args.sensitivity, threads=args.threads, iterate=args.iterate, paired=args.paired, depth=args.depth, tmp_dir=args.tmp_dir)
    logger.info("Shotgun pipeline finished")


if __name__ == "__main__":
	main(sys.argv[1:])
