#!/usr/bin/env python
import os
import sys
import csv
import argparse
import subprocess
from loguru import logger


def run_pipeline_on_sra_table(inputname, parallel=True, pipeline_script='~/git/shotgun/shotgun_pipeline.py', skip_if_exists=True, start_step=0, database='~/databases/uniref/db-uniref50.dmnd'):
        '''Run the sample pipeline on all samples listed in the SRA metadata table
        
        Parameters
        ----------
        inputname: str
                the SraRunInfo.txt file. A table containing a column with Run_s/Run/acc column that contains the SRR accession numbers.
        parallel: bool, optional
                if true, run samples in parallel
        pipeline_script: str, optional
                path to the shotgun pipeline script (shotgun_pipeline.py)
        skip_if_exists: bool, optional
                if true, skip each processing step if the relevant output file already exists
        start_step: int, optional
                step to start from (0: download, 1: clean, 2: convert to fasta, 3: align, 4: split)
        database: str, optional
                location of the diamond uniref database to use for alignment
        '''
        pipeline_script = os.path.expanduser(pipeline_script)
        logger.info(f"Running pipeline on SRA table {inputname} with parallel={parallel}")
        # get the delimiter
        with open(inputname) as csvfile:
                xx = csv.Sniffer()
                res = xx.sniff(csvfile.readline(), delimiters=',\t')
                delimiter = res.delimiter
                logger.info(f"Detected delimiter {delimiter}")
        
        ifile = csv.DictReader(open(inputname, 'r'), delimiter=delimiter)
        for cline in ifile:
                if 'Run_s' in cline:
                        csamp = cline['Run_s']
                elif 'Run' in cline:
                        csamp = cline['Run']
                elif 'acc' in cline:
                        csamp = cline['acc']
                logger.info(f"Processing sample {csamp} from SRA table")
                cmd = [sys.executable, pipeline_script, '-a', csamp, '--start-step', str(start_step), '--database', database]
                if skip_if_exists:
                    cmd += ['--skip-if-exists']
                if parallel:
                    subprocess.Popen(cmd)
                else:
                    subprocess.call(cmd)
                
        logger.info("Finished processing all samples from SRA table")
        return


def main(argv):
    parser = argparse.ArgumentParser(description='process_sra_table', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', help='SRA table file name')
    parser.add_argument('-p', '--parallel', help='Process samples in parallel', action='store_true', default=False)
    parser.add_argument('--start-step', type=int, help='Step to start from (0: download, 1: clean, 2: convert to fasta, 3: align, 4: split)', default=0)
    parser.add_argument('--skip-if-exists', action='store_true', help='Skip processing steps if output files already exist', default=True)
    parser.add_argument('--pipeline-script', type=str, help='Path to the shotgun pipeline script', default='~/git/shotgun/shotgun_pipeline.py')
    parser.add_argument('--database', type=str, help='Path to the database to use for alignment', default='~/databases/uniref/db-uniref50.dmnd')

    args = parser.parse_args(sys.argv[1:])
    # add file logging
    logger.add("shotgun_pipeline.log", rotation="10 MB")
    logger.info("Starting shotgun pipeline")
    if args.input:
        run_pipeline_on_sra_table(args.input, parallel=args.parallel, skip_if_exists=args.skip_if_exists, start_step=args.start_step, pipeline_script=args.pipeline_script, database=args.database)
    logger.info("Shotgun pipeline finished")


if __name__ == "__main__":
	main(sys.argv[1:])
