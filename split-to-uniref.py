#!/usr/bin/env python

__version__ = '2026.05.13'

import os
import glob
from collections import defaultdict
import argparse
import sys
import subprocess

import pandas as pd
import numpy as np

try:
    from loguru import logger
except ImportError:
    import logging
    logger = logging.getLogger(__name__)


def rev_comp(seq):
    '''Return the reverse complement of a DNA sequence.'''
    comp_dict = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp_dict)[::-1]


def align(ref_reads):
    '''Create the aligned sequences for all reads, by adding '-' to the left and right of the read sequence, according to the position of the read on the reference.
    Handles both forward and reverse reads

    Parameters:
    ref_reads (pd.DataFrame): the dataframe containing the reads aligned to the reference, with the following columns: 'qseqid', 'qstart', 'qend', 'qframe', 'qstrand', 'sstart', 'send', 'qseq', 'qseq_translated', 'length', 'evalue', 'bitscore', 'gapopen', 'cigar', 'full_qseq'

    Returns:
    all_prot (list): list of all aligned sequences, with '-' added to the left and right of the read sequence, according to the position of the read on the reference
    all_ids (list): list of all read ids
    all_info (list): list of all read info, including the strand and the aligned sequence
    '''
    all_prot = []
    all_ids = []
    all_info = []
    max_pos = ref_reads['sstart'].max()*3 + 150
    for cpos,crow in ref_reads.iterrows():
        cseq = crow['full_qseq']
        cframe = crow['qframe']
        if crow['qstrand']=='-':
            cseq = rev_comp(cseq[crow['qend']:crow['qstart']])
            cfull_seq = cseq
            cfull_seq = '-'*((crow['sstart'])*3) + cfull_seq
            # complete to the right with '-' up to the max position
            cfull_seq = cfull_seq + '-'*(max_pos - len(cfull_seq))
            # continue
        else:
            start_pos = crow['qstart']-1
            end_pos = crow['qend']-1
            cfull_seq = cseq[start_pos:end_pos]
            cfull_seq = '-'*(crow['sstart']*3) + cfull_seq
            cfull_seq = cfull_seq + '-'*(max_pos - len(cfull_seq))
        all_prot.append(cfull_seq)
        all_ids.append(crow['qseqid'])
        all_info.append([crow['qstrand'], cfull_seq])
    return all_prot, all_ids, all_info


def calc_windows(all_prot, window_size=25):
    '''Count the number of occurances of each sequence of length window_size at each position in the aligned sequences.

    Parameters:
    -----------
    all_prot (list): list of all aligned sequences, with '-' added to the left and right of the read sequence, according to the position of the read on the reference
    window_size (int): the size of the window to count the sequences

    Returns:
    --------
    pos_nums (defaultdict): a nested dictionary where the keys are the positions in the aligned sequences and the values are dictionaries where the keys are the sequences of length window_size and the values are the number of occurances of each sequence at that position
    '''
    pos_nums = defaultdict(lambda: defaultdict(int))
    for cseq in all_prot:
        for i in range(len(cseq)-window_size+1):
            tseq = cseq[i:i+window_size]
            if '-' in tseq:
                continue
            pos_nums[i][tseq] += 1
    return pos_nums


def denoise_reads(pos_nums, noise_level=0.05):
    '''Denoise the reads by removing sequences that are likely to be errors, based on their abundance and their similarity to more abundant sequences.
    Denoising is done on the window length sequences, by comparing each sequence to more abundant sequences at the same position and reducing their count if they are similar (hamming distance <= 2) to a more abundant sequence, by a factor of noise_level.

    Parameters:
    -----------
    pos_nums (defaultdict): a nested dictionary where the keys are the positions in the aligned sequences and the values are dictionaries where the keys are the sequences of length window_size and the values are the number of occurances of each sequence at that position
    noise_level (float): the threshold for considering a sequence as noise, as a fraction of the count of the most abundant sequence at that position

    Returns:
    --------
    out_pos_nums (defaultdict): a nested dictionary where the keys are the positions in the aligned sequences and the values are dictionaries where the keys are the denoised sequences of length window_size and the values are the number of occurances of each denoised sequence at that position
    '''
    import numpy as np
    
    out_pos_nums = defaultdict(lambda: defaultdict(int))
    
    for cpos, cseqs in pos_nums.items():
        if len(cseqs) == 0:
            continue
            
        # Convert to lists sorted by count (descending)
        sorted_items = sorted(cseqs.items(), key=lambda item: item[1], reverse=True)
        sequences = [item[0] for item in sorted_items]
        counts = np.array([item[1] for item in sorted_items], dtype=float)
        
        # Convert sequences to numpy character arrays for vectorized operations
        seq_arrays = np.array([list(seq) for seq in sequences])
        n_seqs = len(sequences)
        
        # Process each sequence in order of abundance
        for i in range(n_seqs):
            current_count = counts[i]
            if current_count <= 0:
                continue
                
            # Vectorized hamming distance calculation for all other sequences
            current_seq = seq_arrays[i]
            hamming_distances = np.sum(seq_arrays != current_seq, axis=1)
            
            # Find sequences with hamming distance <= 2 (excluding self)
            similar_mask = (hamming_distances <= 2) & (np.arange(n_seqs) != i)
            
            # Apply noise reduction to similar sequences
            reduction = noise_level * current_count
            counts[similar_mask] = np.maximum(0, counts[similar_mask] - reduction)
        
        # Convert back to integer dictionary, keeping only non-zero counts
        result_dict = {}
        for seq, count in zip(sequences, counts):
            int_count = int(count)
            if int_count > 0:
                result_dict[seq] = int_count
                
        out_pos_nums[cpos] = result_dict
    
    return out_pos_nums


def parse_results(res_dir, uniref_id, window_size=25, min_reads=0, only_one=False, min_files=5, percentile=95):
    '''Count the number of variants at each position for a given uniref_id in the result files in the specified directory, and optionally plot the distribution of the number of variants at each position

    Parameters
    ----------
    res_dir: str
        the directory containing the per-sample splits/ids.txt files (from the shotgun_pipeline.py)
    unired_id: str
        the uniref id of the sequences to check (with .txt appended)
    window_size: int
        the size of the window to count the sequences (default: 25)
    min_reads: int
        the minimum number of reads required to consider a sequence as valid (default: 0 since we use denoising instead)
    only_one: bool
        whether to only parse one file (for testing purposes, default: False)
    min_files: int
        the minimum number of files (samples) in which the uniref_id must be present to consider the results valid (default: 5)
    percentile: int, optional
        the percentile of the number of reads per position to return (use 100 to get the maximal number of variants)
    
    Returns
    -------
    all_num: defaultdict(float)
        a dictionary where the keys are the sample names and the values are the specified percentile of the number of variants at each position for that sample for the uniref_id
    all_tot_reads: defaultdict(int)
        a dictionary where the keys are the sample names and the values are the total number of reads for that sample mapped to the uniref_id
    id_len: int
        the length of the uniref_id sequence (in nucleotides)
    '''
    result_files = glob.glob(os.path.join(res_dir, f"*-splits/{uniref_id}"))
    all_num = defaultdict(float)
    all_tot_reads = defaultdict(int)
    if len(result_files) < min_files:
        logger.info(f'Found only {len(result_files)} result files for uniref_id {uniref_id} in directory {res_dir}. Expected at least {min_files} files. Check if the files are in the correct directory and have the correct naming convention.')
        return all_num, all_tot_reads,0
    id_len = 0
    columns = 'sseqid qseqid qstart qend qframe qstrand sstart send qseq qseq_translated length evalue bitscore gapopen cigar full_qseq'.split(' ')
    for res_file in result_files:
        logger.debug(f'parsing file {res_file}')
        res_df = pd.read_csv(res_file, sep='\t', header=None)
        res_df.columns = columns
        logger.debug(f'found {len(res_df)} matching reads in file {res_file}')
        all_prot, all_ids, all_info = align(res_df)
        pos_nums = calc_windows(all_prot, window_size=window_size)
        # denoise the reads to remove sequences that are read errors
        pos_nums = denoise_reads(pos_nums)

        num_diff = [np.sum(np.array(list(pos_nums[i].values()))>min_reads) for i in range(len(all_prot[0]))]
        num_unique = np.percentile(num_diff, percentile, axis=0)
        id_len = len(all_prot[0])
        logger.debug(f'sample {os.path.basename(res_file)}: {percentile}% unique: {num_unique}, total reads: {len(res_df)}')
        sample_name = res_file.split('/')[-2]
        all_num[sample_name] = num_unique
        all_tot_reads[sample_name] = len(res_df)
        if only_one:
            break
    return all_num, all_tot_reads, id_len


def create_reads_table(base_dir, output_file='reads.txt'):
    '''Create a table with the number of reads for each sample, by counting the number of sequences in the .clean.fasta files in the sample splits directories, and write it to a file called reads.txt in the base directory. The table will have two columns: sampleid and reads, where sampleid is the name of the sample (without the -1.clean.fasta suffix) and reads is the number of reads in that sample.
    for speed, runs grep -c ">" *.clean.fasta > reads.txt

    Parameters
    ----------
    base_dir : str
        The base directory where the sample splits are located
    output_file : str
        The name of the output file to write the reads table to (default: 'reads.txt')

    Returns
    -------
    str        The path to the output file containing the reads table
    '''
    logger.info(f'Creating reads table {output_file} for base directory {base_dir}')
    with open(output_file, 'w') as f:
        dir_name = os.path.join(base_dir, '*.clean.fasta')
        subprocess.run(f'grep -c ">" {dir_name} > {output_file}', shell=True)
        # read the output file into a pandas dataframe
        reads = pd.read_csv(output_file, sep=':', header=None, names=['sampleid', 'reads'])
        # remove the path and the -1.clean.fasta suffix from the sampleid column
        reads['sampleid'] = reads['sampleid'].apply(lambda x: os.path.basename(x))
        # write the dataframe back to the output file
        reads.to_csv(output_file, sep='\t', index=False)

    logger.debug(f'Reads table {output_file} created')
    return output_file


def get_ids_list(base_dir):
    '''Get the uniref ids and their counts across the samples
    
    Parameters
    ----------
    base_dir : str
        The base directory where the sample splits are located

    Returns
    -------
    dict of uniref ids appearing in the sample splits (key) and their counts (value)
    '''
    uniref_ids = defaultdict(int)
    for name in glob.glob(f'{base_dir}/*-splits/UniRef50_*.txt'):
        uniref_ids[name.split('/')[-1]] += 1

    # sort the dict
    uniref_ids = dict(sorted(uniref_ids.items(), key=lambda item: item[1], reverse=True))
    logger.info(f'Found {len(uniref_ids)} unique uniref ids across the samples')
    return uniref_ids

def split_to_uniref(base_dir='.', metadata_file='map.txt',reads_file=None, output_file='uniref-table.txt', min_samples_per_uniref=5, window_size=50, num_ids=None, normalize_reads_per_gene=False):
    if reads_file is None:
        reads_file = create_reads_table(base_dir)

    # metadata = pd.read_csv(base_dir+'/table.csv',sep=',', index_col='Run')
    metadata = pd.read_csv(os.path.join(base_dir, metadata_file), sep='\t', index_col='Run')
    reads = pd.read_csv(reads_file, sep='\t', header=None, names=['sampleid', 'reads'])

    # get the list of all the sample files
    all_samples = []
    split_dirs = glob.glob(os.path.join(base_dir, f"*-splits"))
    logger.info(f'Found {len(split_dirs)} sample split directories in base directory {base_dir}')
    for cname in split_dirs:
        baseid = os.path.basename(cname).split('-splits')[0]
        sid = f'{baseid}-1.clean.fasta'
        if sid in reads['sampleid'].values:
            all_samples.append(baseid)
        else:
            logger.warning(f'Sample {sid} not found in reads.txt, skipping it.')
    logger.info(f'Number of samples found: {len(all_samples)} out of {len(metadata)} samples in metadata')

    uniref_ids = get_ids_list(base_dir)
    ids_list = list(uniref_ids.keys())
    # permute the uniref ids list to get a random order of checking the uniref ids so we get better estimation of performance
    np.random.shuffle(ids_list)

    logger.info(f'Processing {len(ids_list)} uniref ids across {len(all_samples)} samples')
    with open(output_file, 'w') as f:
        f.write('uniref_id')
        for csample in all_samples:
            f.write(f'\t{csample}')
        f.write('\n')

        for cidx, uniref_id in enumerate(ids_list):
            if num_ids is not None and cidx >= num_ids:
                logger.info(f'Processed {cidx} uniref ids, stopping as num_ids is set to {num_ids}')
                break
            logger.debug(f'Processing uniref_id {cidx+1}/{len(ids_list)}: {uniref_id}')
            all_num, all_tot_reads, all_len = parse_results(base_dir, uniref_id, window_size=window_size, min_files=min_samples_per_uniref)
            if len(all_num) == 0:
                continue
            f.write(uniref_id)
            for csample in all_samples:
                csample_id = f'{csample}-splits'
                if csample_id in all_num:
                    num_variants = all_num[csample_id]
                    if normalize_reads_per_gene:
                        if csample_id in all_tot_reads and all_tot_reads[csample_id] > 0:
                            num_variants = num_variants / all_tot_reads[csample_id]
                        else:
                            logger.warning(f'No reads found for sample {csample_id} in uniref_id {uniref_id}, cannot normalize by reads, keeping original variant count')
                            num_variants = 0
                    f.write(f'\t{num_variants}')
                else:
                    f.write('\t0')
            f.write('\n')

def main(argv):
    parser = argparse.ArgumentParser(description='split-to-uniref version ' + __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('-b', '--base-dir', help='Base directory where the sample splits are located', default='.')
    parser.add_argument('-m', '--metadata-file', help='Metadata file (tab-separated) with sample information, must contain a "Run" column with sample ids', default='map.txt')
    parser.add_argument('-r', '--reads-file', help='File (tab-separated) with sample ids and their read counts, must contain "sampleid" and "reads" columns. tsv. If not provided, will auto-generate', default=None)
    parser.add_argument('-o', '--output-file', help='Output file to write the uniref table to', default='uniref-table.txt')
    parser.add_argument('--min-samples-per-uniref', type=int, help='Minimum number of samples a uniref_ID must be present in to be included in output table', default=5)
    parser.add_argument('--window-size', type=int, help='Window size for counting sequence variants', default=50)
    parser.add_argument('--log-level', type=str, help='Logging level (DEBUG, INFO, WARNING, ERROR)', default='INFO')
    parser.add_argument('--num-ids', type=int, help='Number of uniref ids to process (for testing purposes, None to process all)', default=None)
    parser.add_argument('--normalize-reads-per-gene', action='store_true', help='Whether to normalize each variant count by the total number of reads mapped to the gene (in each sample)')

    args = parser.parse_args(sys.argv[1:])
    logger.remove()  # Remove default logger
    logger.add(sys.stderr, level=args.log_level)  # Add new logger with specified log level

    logger.info("Starting split-to-uniref pipeline")
    split_to_uniref(base_dir=args.base_dir, metadata_file=args.metadata_file, reads_file=args.reads_file, output_file=args.output_file, min_samples_per_uniref=args.min_samples_per_uniref, window_size=args.window_size, num_ids=args.num_ids, normalize_reads_per_gene=args.normalize_reads_per_gene)
    logger.info("Split-to-uniref pipeline finished")


if __name__ == "__main__":
	main(sys.argv[1:])
