#!/usr/bin/env python

__version__ = '2026.06.15'

import os
import glob
from collections import defaultdict
import argparse
import sys
import subprocess
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
import sqlite3

import pandas as pd
import numpy as np
import scipy as sp

try:
    from loguru import logger
except ImportError:
    import logging
    logger = logging.getLogger(__name__)


def uniref_to_db(uniref_file, database_file, ids_file=None, delete_database=False, expected_lines_in_file=6146909076):
    '''
    Convert a UniRef XML file to a SQLite database.

    Parameters:
    - uniref_file: Path to the UniRef XML file to process.
    - database_file: Path to the SQLite database file to create or update.
    - ids_file: Path to a file containing UniRef IDs to extract (one ID per line). If not provided, all entries will be processed.
    - delete_database: If True, delete the existing database before processing.
    - expected_lines_in_file: Expected number of lines in the UniRef XML file (used for progress reporting).

    Stores the following fields in the database:
    - uniref_id: The UniRef ID of the entry.
    - taxons: Comma-separated list of taxons associated with the entry.
    - names: Comma-separated list of names associated with the entry.
    - go_ids: Comma-separated list of GO IDs associated with the entry.
    - accessions: Comma-separated list of UniProtKB accessions associated with the entry
    - avg_length: Average length of the sequences in the entry.
    '''
    cid = 'na'
    num_ids = 0
    taxons = set()
    names = set()
    go_ids = set()
    accessions = set()
    processing_entry = False
    processing_relevant = False
    found_count = 0
    num_lines = 0
    lengths = []

    # load the ids_file if provided (we expect it to be a tsv file with first column being the UniRef IDs to extract)
    if ids_file is not None:
        logger.info(f"Loading target IDs from file: {ids_file}")
        target_ids = set()
        with open(ids_file, "r", encoding="utf-8") as f:
            for line in f:
                target_ids.add(line.strip().split('\t')[0])
        logger.info(f"Loaded {len(target_ids)} target IDs from file.")
    else:
        target_ids = None
        logger.debug("No target IDs file provided; processing all entries.")

    # Connect to database
    logger.debug(f"Connecting to database: {database_file}")
    conn = sqlite3.connect(database_file)
    cursor = conn.cursor()

    if delete_database:
        logger.info(f"Deleting existing database: {database_file}")
        cursor.execute("DROP TABLE IF EXISTS entries;")
        conn.commit()
    # create the table
    logger.debug('creating table entries')
    cursor.execute("CREATE TABLE IF NOT EXISTS entries (uniref_id TEXT, taxons TEXT, names TEXT, go_ids TEXT, accessions TEXT, avg_length REAL);")

    batch_size = 100_000
    batch = []

    #Apply High-Performance Speed Settings
    cursor.execute("PRAGMA journal_mode = WAL;")  # WAL mode is highly efficient for 60M scales
    cursor.execute("PRAGMA synchronous = NORMAL;") # Balances incredible speed with safety
    cursor.execute("PRAGMA cache_size = -1000000;") # Allocates 1GB of RAM cache

    logger.info(f"Processing UniRef XML file: {uniref_file}")
    with open(uniref_file, "r", encoding="utf-8") as f:
        for cline in f:
            num_lines += 1
            if num_lines % 10000000 == 0:
                print(f"Processed {num_lines} lines ({100*num_lines/expected_lines_in_file:.2f}%), found {found_count} target entries so far...")
            # if a new entry starts, extract its ID and continue to the next line
            if cline.startswith('<entry id'):
                processing_entry = True
                cid = cline.split('"',maxsplit=2)[1]
                if target_ids is not None:
                    processing_relevant = (cid in target_ids)
                else:
                    processing_relevant = True
                if processing_relevant:
                    found_count += 1
                    num_ids += 1
                    lengths = []
                    taxons = set()
                    names = set()
                    go_ids = set()
                    accessions = set()
                    continue
            if processing_entry and cline.startswith('</entry>'):
                if processing_relevant:
                    batch.append((cid, ','.join(list(taxons)), ','.join(list(names)), ','.join(list(go_ids)), ','.join(list(accessions)), np.mean(lengths) if lengths else 0))
                    if len(batch) >= batch_size:
                        print('adding batch of', len(batch), 'entries to database...')
                        print(batch[0])
                        cursor.executemany("INSERT INTO entries VALUES (?, ?, ?, ?, ?, ?);", batch)
                        conn.commit()
                        batch.clear()

                processing_entry = False
                processing_relevant = False
                if target_ids is not None:
                    if found_count >= len(target_ids):
                        print("All target entries found. Stopping early.")
                        break
                continue
            if not processing_relevant:
                continue
            if cline.startswith('<property type="common taxon"'):
                ctaxon = cline.split('"',maxsplit=4)[3]
                taxons.add(ctaxon)
                continue
            if cline.startswith('<property type="source organism"'):
                ctaxon = cline.split('"',maxsplit=4)[3]
                taxons.add(ctaxon)
                continue
            if cline.startswith('<name>'):
                cname = cline.split('>',maxsplit=1)[1].split('<',maxsplit=1)[0]
                names.add(cname)
                continue
            if cline.startswith('<property type="protein name"'):
                cname = cline.split('"',maxsplit=4)[3]
                names.add(cname)
                continue
            if cline.startswith('<property type="GO'):
                cgo_id = cline.split('"',maxsplit=4)[3]
                go_ids.add(cgo_id)
                continue
            if cline.startswith('<property type="UniProtKB accession'):
                caccession = cline.split('"',maxsplit=4)[3]
                accessions.add(caccession)
                continue
            if cline.startswith('<property type="length"'):
                clength = int(cline.split('"',maxsplit=4)[3])
                lengths.append(clength)
                continue

    # final batch insert for any remaining entries
    if len(batch) > 0:
        cursor.executemany("INSERT INTO entries VALUES (?, ?, ?, ?, ?, ?);", batch)
        conn.commit()
    # Create index after all insertions for optimal performance
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_uniref_id ON entries (uniref_id);")



def main(argv):
    parser = argparse.ArgumentParser(description='uniref-to-db.py version ' + __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('--uniref-file', type=str, required=True, help='Path to the UniRef XML file to process.')
    parser.add_argument('--database-file', type=str, default='uniref.db', help='Path to the SQLite database file to create or update.')
    parser.add_argument('--ids-file', type=str, default=None, help='Path to a file containing UniRef IDs to extract (one ID per line). If not provided, all entries will be processed.')
    parser.add_argument('--delete-database', action='store_true', help='Delete the existing database before processing.')
    parser.add_argument('--expected-lines-in-file', type=int, default=6146909076, help='Expected number of lines in the UniRef XML file (used for progress reporting).')
    parser.add_argument('--log-level', type=str, default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level.')

    args = parser.parse_args(sys.argv[1:])
    logger.remove()  # Remove default logger
    logger.add(sys.stderr, level=args.log_level)  # Add new logger with specified log level

    logger.info("Starting split-to-uniref pipeline")
    uniref_to_db(uniref_file=args.uniref_file, ids_file=args.ids_file, database_file=args.database_file, delete_database=args.delete_database, expected_lines_in_file=args.expected_lines_in_file)
    logger.info("Split-to-uniref pipeline finished")


if __name__ == "__main__":
	main(sys.argv[1:])
