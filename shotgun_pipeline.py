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

def get_sample(sample_id, sra_path='~/bin/sratoolkit.3.2.0-centos_linux64/bin', skip_if_exists=True, log_file='process.log'):
    '''Download a single sample from SRA given its SRA ID

    Parameters
    ----------
    sample_id: str
            SRA sample ID (SRRxxxxxx)
    sra_path: str, optional
            path to the sra-toolkit binaries
    skip_if_exists: bool, optional
            if true, skip downloading if the sample already exists
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
    params = [os.path.join(sra_path, 'fasterq-dump'), './' + sample_id, '--split-files', '--threads', '4']
    logger.debug(f"Running command: {' '.join(params)}")
    with open(log_file, 'a') as logfile:
        res = subprocess.call(params, stdout=logfile, stderr=logfile)
    if res != 0:
        logger.error(f"fasterq-dump failed with return code {res}")
        raise RuntimeError("fasterq-dump execution failed")
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
    params = [fastp_path, '-i', input_r1, '-o', output_r1, '--length_required', '50', '--qualified_quality_phred', '25']
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


def align_to_uniref(sample_id, diamond_db='~/databases/uniref/db-uniref90.dmnd', diamond_path='~/bin/diamond', skip_if_exists=True, log_file='process.log'):
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
    '''
    output_file = f"{sample_id}-aligned.txt"
    input_fasta = f"{sample_id}-1.clean.fasta"
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
        '--fast',
        '--min-orf', '1',
        '--query', input_fasta,
        '--strand', 'both',
        '--max-target-seqs', '1',
        '--un', output_file + '.unmatched.fasta',
        '--threads', '10',
        '--outfmt', '6'
    ]
    command.extend(diamond_output_columns)
    with open(log_file, 'a') as logfile:
        res = subprocess.run(command, stdout=logfile, stderr=logfile)
    if res.returncode != 0:
        logger.error(f"DIAMOND failed with return code {res.returncode}")
        raise RuntimeError("DIAMOND execution failed")
    logger.info("DIAMOND completed successfully")
    return


def split_to_uniref(sample_id, skip_if_exists=True, log_file='process.log'):
    '''Split DIAMOND output file into separate files per UniRef ID

    Parameters
    ----------
    diamond_output: str
        path to DIAMOND output file
    out_dir: str
        output directory to store per-UniRef ID files
    skip_if_exists: bool, optional
        if true, skip splitting if the output directory already exists and is not empty
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
    uniref_ids = defaultdict(int)
    with open(diamond_output, 'r') as f:
        for line in f:
            parts = line.strip().split('\t', maxsplit=1)
            uniref_id = parts[0]
            uniref_ids[uniref_id] += 1
            with open(os.path.join(out_dir, f"{uniref_id}.txt"), 'a') as out_f:
                out_f.write(line)
    logger.info(f"Completed splitting DIAMOND output into {len(uniref_ids)} uniref ids, total lines: {sum(uniref_ids.values())}")
    return


def sample_pipeline(sample_id, skip_if_exists=True):
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
    '''
    log_file = f'process-{sample_id}.log'
    logger.info(f"Processing sample {sample_id}")
    # Step 1: Download the sample
    get_sample(sample_id, skip_if_exists=skip_if_exists, log_file=log_file)
    # Step 2: Clean the sample
    clean_sample(sample_id, skip_if_exists=skip_if_exists, log_file=log_file)
    # Step 3: Convert to fasta
    convert_to_fasta(sample_id, skip_if_exists=skip_if_exists, log_file=log_file)
    # Step 4: Align to UniRef
    align_to_uniref(sample_id, skip_if_exists=skip_if_exists, log_file=log_file)
    # Step 5: Split to per-UniRef ID files
    split_to_uniref(sample_id, skip_if_exists=skip_if_exists, log_file=log_file)
    logger.info(f"Finished processing sample {sample_id}")
    return


def main(argv):
    parser = argparse.ArgumentParser(description='Shotgun pipeline version ' + __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--accession', help='SRA sample accession to process')

    args = parser.parse_args(sys.argv[1:])
    # add file logging
    logger.add("shotgun_pipeline.log", rotation="10 MB")
    logger.info("Starting shotgun pipeline")
    if args.accession:
        sample_pipeline(args.accession)
    logger.info("Shotgun pipeline finished")


if __name__ == "__main__":
	main(sys.argv[1:])
