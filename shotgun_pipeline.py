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

def get_sample(sample_id, sra_path='~/bin/sratoolkit.3.2.0-centos_linux64/bin'):
    '''Download a single sample from SRA given its SRA ID

    Parameters
    ----------
    sample_id: str
            SRA sample ID (SRRxxxxxx)
    sra_path: str, optional
            path to the sra-toolkit binaries
    '''
    sra_path = os.path.expanduser(sra_path)
    logger.info(f"Downloading sample {sample_id}")
    # prefetch
    params = [os.path.join(sra_path, 'prefetch'), sample_id]
    logger.debug(f"Running command: {' '.join(params)}")
    subprocess.call(params)
    logger.debug(f"Prefetched sample {sample_id}")
    # fasterq-dump
    params = [os.path.join(sra_path, 'fasterq-dump'), './' + sample_id, '--split-files', '--threads', '4']
    logger.debug(f"Running command: {' '.join(params)}")
    subprocess.call(params)
    logger.info(f"Converted sample {sample_id} to fastq")
    return


def clean_sample(sample_id, fastp_path='~/bin/fastp'):
    '''Clean a single sample using fastp

    Parameters
    ----------
    sample_id: str
            SRA sample ID (SRRxxxxxx)
    fastp_path: str, optional
            path to the fastp binary
    '''
    logger.info(f"Cleaning sample {sample_id}")
    fastp_path = os.path.expanduser(fastp_path)
    # run fastp
    input_r1 = f"{sample_id}_1.fastq"
    output_r1 = f"{sample_id}-1.clean.fastq"
    params = [fastp_path, '-i', input_r1, '-o', output_r1, '--length_required', '50', '--qualified_quality_phred', '25']
    logger.debug(f"Running command: {' '.join(params)}")
    subprocess.call(params)
    logger.debug(f"Cleaned sample {sample_id}")
    return


def convert_to_fasta(sample_id, seqtk_path='~/bin/seqtk'):
    '''Convert cleaned fastq files to fasta using seqtk

    Parameters
    ----------
    sample_id: str
            SRA sample ID (SRRxxxxxx)
    seqtk_path: str, optional
            path to the seqtk binary
    '''
    seqtk_path = os.path.expanduser(seqtk_path)
    logger.info(f"Converting sample {sample_id} to fasta")
    input_fastq = f"{sample_id}-1.clean.fastq"
    output_fasta = f"{sample_id}-1.clean.fasta"
    params = [seqtk_path, 'seq', '-a', input_fastq]
    logger.debug(f"Running command: {' '.join(params)}")
    with open(output_fasta, 'w') as outfile:
        subprocess.call(params, stdout=outfile)
    logger.debug(f"Converted sample {sample_id} to fasta")
    return


def align_to_uniref(input_fasta, output_file, diamond_db='~/databases/uniref/db-uniref90.dmnd', diamond_path='~/bin/diamond'):
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
    '''
    logger.info(f"Running DIAMOND on {input_fasta}, outputting to {output_file} using DB {diamond_db}")
    diamond_path = os.path.expanduser(diamond_path)
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
    res = subprocess.run(command)
    if res.returncode != 0:
        logger.error(f"DIAMOND failed with return code {res.returncode}")
        raise RuntimeError("DIAMOND execution failed")
    logger.info("DIAMOND completed successfully")
    return

def sample_pipeline(sample_id):
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
    '''
    logger.info(f"Processing sample {sample_id}")
    # Step 1: Download the sample
    get_sample(sample_id)
    # Step 2: Clean the sample
    clean_sample(sample_id)
    # Step 3: Convert to fasta
    convert_to_fasta(sample_id)
    # Step 4: Align to UniRef
    align_to_uniref(sample_id)
    logger.info(f"Finished processing sample {sample_id}")
    return


def GetSRA(inputname, path, skipifthere=False, fastq=False, delimiter=None, outdir='fasta', split_files=False):
        '''Get all the samples from the SRA. Using the Metadata (runinfo metadata) table SraRunTable.txt from the run browser

        Parameters
        ----------
        inputname: str
                the SraRunInfo.txt file. A table containing a column with Run_s/Run/acc column that contains the SRR accession numbers.
        path: str
                path to the SraToolKit binary directory
        skipifthere: bool, optional
                if true, do not download files that already exist
        fastq: bool, optional
                if true, download fastq instead of fasta
        delimiter: str or None, optional
                delimiter for the table. If none, autodetect
        outdir: str, optional
                name of the output directory for the downloads
        split_files: bool, optional
                if True, split the samples into forward and reverse reads

        Returns
        -------
        num_files: int
                number of files downloaded
        '''
        logger.info('Starting GetSRA with input file %s' % inputname)
        # create output directory
        if not os.path.exists(outdir):
                os.makedirs(outdir)

        # get the delimiter if not provided
        if delimiter is None:
                with open(inputname) as csvfile:
                        xx = csv.Sniffer()
                        res = xx.sniff(csvfile.readline(), delimiters=',\t')
                        delimiter = res.delimiter
                        logger.info('Detected delimiter %s' % delimiter)

        ifile = csv.DictReader(open(inputname, 'r'), delimiter=delimiter)
        num_files = 0
        num_skipped = 0
        for cline in ifile:
                if 'Run_s' in cline:
                        csamp = cline['Run_s']
                elif 'Run' in cline:
                        csamp = cline['Run']
                elif 'acc' in cline:
                        csamp = cline['acc']
                num_files += 1
                
                if skipifthere:
                        if os.path.isfile(os.path.join(outdir, csamp) + '.fasta'):
                                logger.info("skipping sample %s. file exists" % csamp)
                                continue

                logger.info("getting file %s" % csamp)
                params = [os.path.join(path, 'fastq-dump'), '--disable-multithreading']
                params += ['--outdir', outdir]
                if split_files:
                        params += ['--split-files']
                if not fastq:
                        params += ['--fasta', '0']
                params += [csamp]
                logger.info(f"Running command: {' '.join(params)}")
                subprocess.call(params)
                print("got file %s" % csamp)
        print('got %d files.' % num_files)
        return num_files


def main(argv):
    parser = argparse.ArgumentParser(description='Shotgun pipeline version ' + __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--accession', help='SRA sample accession to process')

    args = parser.parse_args(sys.argv[1:])
    logger.info("Starting shotgun pipeline")
    if args.accession:
        sample_pipeline(args.accession)
    logger.info("Finished shotgun pipeline")


if __name__ == "__main__":
	main(sys.argv[1:])
