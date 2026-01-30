#!/usr/bin/env python
import os
import sys
import csv
import argparse
import subprocess
from loguru import logger


def run_pipeline_on_sra_table(inputname, parallel=True):
        '''Run the sample pipeline on all samples listed in the SRA metadata table
        
        Parameters
        ----------
        inputname: str
                the SraRunInfo.txt file. A table containing a column with Run_s/Run/acc column that contains the SRR accession numbers.
        parallel: bool, optional
                if true, run samples in parallel
        '''
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
                if parallel:
                    subprocess.Popen([sys.executable, __file__, '-a', csamp])
                else:
                    subprocess.call([sys.executable, __file__, '-a', csamp])
                
        logger.info("Finished processing all samples from SRA table")
        return


def main(argv):
    parser = argparse.ArgumentParser(description='process_sra_table', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', help='SRA table file name')
    parser.add_argument('-p', '--parallel', help='Process samples in parallel', action='store_true', default=False)

    args = parser.parse_args(sys.argv[1:])
    # add file logging
    logger.add("shotgun_pipeline.log", rotation="10 MB")
    logger.info("Starting shotgun pipeline")
    if args.input:
        run_pipeline_on_sra_table(args.input, parallel=args.parallel)
    logger.info("Shotgun pipeline finished")


if __name__ == "__main__":
	main(sys.argv[1:])
