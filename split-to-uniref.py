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

import pandas as pd
import numpy as np
import scipy as sp

try:
    from loguru import logger
except ImportError:
    import logging
    logger = logging.getLogger(__name__)


def rev_comp(seq):
    '''Return the reverse complement of a DNA sequence.'''
    comp_dict = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp_dict)[::-1]


def count_window(ref_reads, window_size=25):
    '''Count the number of occurances of each sequence of length window_size at each position in the aligned sequences.

    Parameters:
    ref_reads (pd.DataFrame): the dataframe containing the reads aligned to the reference, with the following columns: 'qseqid', 'qstart', 'qend', 'qframe', 'qstrand', 'sstart', 'send', 'qseq', 'qseq_translated', 'length', 'evalue', 'bitscore', 'gapopen', 'cigar', 'full_qseq'
    window_size (int): the size of the window to count the sequences

    Returns:
    pos_nums (defaultdict): a nested dictionary where the keys are the positions in the aligned sequences and the values are dictionaries where the keys are the sequences of length window_size and the values are the number of occurances of each sequence at that position
    '''
    pos_nums = defaultdict(lambda: defaultdict(int))

    if ref_reads.empty:
        return pos_nums

    for crow in ref_reads.itertuples(index=False):
        qstart = int(crow.qstart)
        qend = int(crow.qend)
        sstart = int(crow.sstart)
        cseq = crow.full_qseq
        if crow.qstrand == '-':
            cseq = rev_comp(cseq[qend:qstart])
        else:
            cseq = cseq[qstart-1:qend]
        
        for i in range(len(cseq) - window_size + 1):
            pos_nums[sstart * 3 + i][cseq[i:i+window_size]] += 1
    return pos_nums




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
    # Some input files may include a header-like row or non-numeric values.
    # Coerce coordinate columns to numeric and skip invalid rows.
    ref_reads = ref_reads.copy()
    for col in ['qstart', 'qend', 'sstart']:
        ref_reads[col] = pd.to_numeric(ref_reads[col], errors='coerce')
    ref_reads = ref_reads.dropna(subset=['qstart', 'qend', 'sstart'])

    # Convert only the scalar max value (not the whole Series) to int.
    max_sstart = ref_reads['sstart'].max()
    max_pos = int(max_sstart) * 3 + 150 if pd.notna(max_sstart) else 150
    for crow in ref_reads.itertuples(index=False):
        qstart = int(crow.qstart)
        qend = int(crow.qend)
        sstart = int(crow.sstart)
        cseq = crow.full_qseq
        if crow.qstrand == '-':
            cseq = rev_comp(cseq[qend:qstart])
            cfull_seq = '-' * (sstart * 3) + cseq
            cfull_seq = cfull_seq + '-' * (max_pos - len(cfull_seq))
        else:
            start_pos = qstart - 1
            end_pos = qend - 1
            cfull_seq = cseq[start_pos:end_pos]
            cfull_seq = '-' * (sstart * 3) + cfull_seq
            cfull_seq = cfull_seq + '-' * (max_pos - len(cfull_seq))
        all_prot.append(cfull_seq)
        all_ids.append(crow.qseqid)
        all_info.append([crow.qstrand, cfull_seq])
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
        i = 0
        seq_len = len(cseq)
        while i < seq_len:
            if cseq[i] == '-':
                i += 1
                continue
            j = i
            while j < seq_len and cseq[j] != '-':
                j += 1
            seg_len = j - i
            if seg_len >= window_size:
                segment = cseq[i:j]
                max_local = seg_len - window_size
                for local in range(max_local + 1):
                    pos_nums[i + local][segment[local:local + window_size]] += 1
            i = j
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
                
        out_pos_nums[cpos] = defaultdict(int, result_dict)
    
    return out_pos_nums


def parse_results(res_dir, uniref_id, methods=['denoise_mean', 'num_reads', 'denoise_deblur', 'num_variants', 'entropy'],
                  window_size=25, min_reads=0, only_one=False, min_files=5, percentile=95, result_files=None, mean_error=0.005):
    '''Count the number of variants at each position for a given uniref_id in the result files in the specified directory, and optionally plot the distribution of the number of variants at each position

    Parameters
    ----------
    res_dir: str
        the directory containing the per-sample splits/ids.txt files (from the shotgun_pipeline.py)
    unired_id: str
        the uniref id of the sequences to check (with .txt appended)
    methods: list of str, optional
        the stats to calculate. default is all supported methods
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
    results_files: list or None, optional
        if not None, parse only the sample split files listed (i.e. 'SRRXXXX-splits/uniref50-YYY.txt')
        if None, parse all -splits directories in res_dir
    
    Returns
    -------
    results: dict of:
        'denoise_mean': defaultdict(float)
            a dictionary where the keys are the sample names and the values are the specified percentile of the number of variants at each position for that sample for the uniref_id
            using the mean noise denoising method (based on mean_error parameter)
        'denoise_deblur': defaultdict(float)
            a dictionary where the keys are the sample names and the values are the specified percentile of the number of variants at each position for that sample for the uniref_id
            using the deblur-like denoising method (based on mean_error parameter)
        'num_variants': defaultdict(float)
            a dictionary where the keys are the sample names and the values are the specified percentile of the number of variants at each position for that sample for the uniref_id
            without denoising
        'num_reads': defaultdict(float)
            a dictionary where the keys are the sample names and the values are the total number of reads for that sample mapped to the uniref_id        
        'entropy': defaultdict(float)
            a dictionary where the keys are the sample names and the values are the entropy of the number of variants at each position for that sample for the uniref_id
    id_len: int
        the length of the uniref_id sequence (in nucleotides)
    '''
    results = {}
    for cmethod in methods:
        results[cmethod] = defaultdict(float)
    if result_files is None:
        result_files = glob.glob(os.path.join(res_dir, f"*-splits/{uniref_id}"))
    all_num = defaultdict(float)
    all_tot_reads = defaultdict(int)
    if len(result_files) < min_files:
        logger.warning(f'Found only {len(result_files)} result files for uniref_id {uniref_id} in directory {res_dir}. Expected at least {min_files} files. Check if the files are in the correct directory and have the correct naming convention.')
        return results, 0
    id_len = 0
    columns = ['qseqid', 'qstart', 'qend', 'qstrand', 'sstart', 'full_qseq']
    usecols = [1, 2, 3, 5, 6, 15]
    for res_file in result_files:
        logger.debug(f'parsing file {res_file}')
        res_df = pd.read_csv(res_file, sep='\t', header=None, usecols=usecols, names=columns)
        # Handle accidental header rows / malformed rows in result files.
        for col in ['qstart', 'qend', 'sstart']:
            res_df[col] = pd.to_numeric(res_df[col], errors='coerce')
        res_df = res_df.dropna(subset=['qstart', 'qend', 'sstart', 'full_qseq'])
        logger.debug(f'found {len(res_df)} matching reads in file {res_file}')
        id_len = len(res_df['full_qseq'].iloc[0]) if len(res_df) > 0 else 0
        if 'num_reads' in methods:
            num_reads = len(res_df)
        if 'num_variants' in methods or 'denoise_deblur' in methods or 'denoise_mean' in methods or 'entropy' in methods:
            # all_prot, all_ids, all_info = align(res_df)
            # pos_nums = calc_windows(all_prot, window_size=window_size)
            pos_nums = count_window(res_df, window_size=window_size)
            # denoise the reads to remove sequences that are read errors
            if 'denoise_deblur' in methods:
                deblur_pos_nums = denoise_reads(pos_nums)
                # id_len = len(all_prot[0])
                num_diff = np.zeros(np.max(list(deblur_pos_nums.keys()))+1, dtype=np.int32)
                for pos, seq_counts in deblur_pos_nums.items():
                    num_diff[pos] = sum(1 for val in seq_counts.values() if val > min_reads)
                num_denoise_deblur = float(np.percentile(num_diff, percentile, axis=0))
            if 'num_variants' in methods:
                # id_len = len(all_prot[0])
                num_diff = np.zeros(np.max(list(pos_nums.keys()))+1, dtype=np.int32)
                for pos, seq_counts in pos_nums.items():
                    num_diff[pos] = sum(1 for val in seq_counts.values())
                num_variants = float(np.percentile(num_diff, percentile, axis=0))
            if 'denoise_mean' in methods:
                # id_len = len(all_prot[0])
                num_diff = np.zeros(np.max(list(pos_nums.keys()))+1, dtype=np.int32)
                for pos, seq_counts in pos_nums.items():
                    num_diff[pos] = len(seq_counts.values())
                    num_diff[pos] = num_diff[pos] - len(seq_counts) * (1-np.pow((1-mean_error),window_size))
                num_denoise_mean = float(np.percentile(num_diff, percentile, axis=0))
            if 'entropy' in methods:
                entropy_values = np.zeros(np.max(list(pos_nums.keys()))+1, dtype=np.float32)
                for pos, seq_counts in pos_nums.items():
                    counts = np.array(list(seq_counts.values()), dtype=np.float32)
                    entropy_values[pos] = sp.stats.entropy(counts)
                entropy = float(np.percentile(entropy_values, percentile, axis=0))

        # logger.debug(f'sample {os.path.basename(res_file)}: {percentile}% unique: {num_unique}, total reads: {len(res_df)}')
        logger.debug(f'sample {os.path.basename(res_file)}: {percentile}% unique:')
        sample_name = res_file.split('/')[-2]
        if 'denoise_mean' in methods:
            results['denoise_mean'][sample_name] = num_denoise_mean
        if 'denoise_deblur' in methods:
            results['denoise_deblur'][sample_name] = num_denoise_deblur
        if 'num_reads' in methods:
            results['num_reads'][sample_name] = num_reads
            logger.info(f'sample {os.path.basename(res_file)}: {num_reads} total reads')
        if 'num_variants' in methods:
            results['num_variants'][sample_name] = num_variants
        if 'entropy' in methods:
            results['entropy'][sample_name] = entropy
        if only_one:
            break
    return results, id_len


def process_uniref_task(task):
    '''Process a single UniRef ID for optional multiprocessing.'''
    uniref_id, base_dir, methods, window_size, min_reads, mean_noise, result_files, min_files = task
    result, id_len = parse_results(
        base_dir,
        uniref_id,
        methods=methods,
        window_size=window_size,
        min_reads=min_reads,
        mean_error=mean_noise,
        result_files=result_files,
        min_files=min_files,
    )
    return uniref_id, result, id_len


def build_uniref_index(base_dir):
    '''Build an index of UniRef ids to matching result files in one filesystem pass.'''
    uniref_to_files = defaultdict(list)
    split_dirs = glob.glob(os.path.join(base_dir, '*-splits'))
    for split_dir in split_dirs:
        with os.scandir(split_dir) as it:
            for entry in it:
                if entry.is_file() and entry.name.startswith('UniRef50_') and entry.name.endswith('.txt'):
                    uniref_to_files[entry.name].append(entry.path)
    return dict(sorted(uniref_to_files.items(), key=lambda item: len(item[1]), reverse=True))


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
    logger.info(f'Building UniRef index for base directory {base_dir}')
    uniref_ids = {uid: len(paths) for uid, paths in build_uniref_index(base_dir).items()}
    logger.info(f'Found {len(uniref_ids)} unique uniref ids across the samples')
    return uniref_ids

def split_to_uniref(base_dir='.', metadata_file='map.txt',reads_file=None, min_samples_per_uniref=5,
                    window_size=50, num_ids=None, min_reads_for_norm=20,
                    num_reads_out=None, denoise_mean_out=None, mean_noise=0.005,
                    denoise_deblur_out=None, variants_out=None, entropy_out=None,
                    jobs=1):

    # if reads_file is None:
    #     reads_file = create_reads_table(base_dir)

    # metadata = pd.read_csv(base_dir+'/table.csv',sep=',', index_col='Run')
    metadata = pd.read_csv(os.path.join(base_dir, metadata_file), sep='\t', index_col='Run')
    # reads = pd.read_csv(reads_file, sep='\t', header=None, names=['sampleid', 'reads'])

    # get the list of all the sample files
    all_samples = []
    split_dirs = glob.glob(os.path.join(base_dir, f"*-splits"))
    logger.info(f'Found {len(split_dirs)} sample split directories in base directory {base_dir}')
    for cname in split_dirs:
        baseid = os.path.basename(cname).split('-splits')[0]
        sid = f'{baseid}-1.clean.fasta'

        all_samples.append(baseid)
        # if sid in reads['sampleid'].values:
        #     all_samples.append(baseid)
        # else:
        #     logger.warning(f'Sample {sid} not found in reads.txt, skipping it.')
    logger.info(f'Number of samples found: {len(all_samples)} out of {len(metadata)} samples in metadata')

    uniref_index = build_uniref_index(base_dir)
    logger.info(f'Found {len(uniref_index)} unique uniref ids across the samples')
    ids_list = list(uniref_index.keys())
    # permute the uniref ids list to get a random order of checking the uniref ids so we get better estimation of performance
    np.random.shuffle(ids_list)

    logger.info(f'Processing {len(ids_list)} uniref ids across {len(all_samples)} samples')
    
    files=[]
    methods = []
    with ExitStack() as stack:
        if num_reads_out is not None:
            reads_file = stack.enter_context(open(num_reads_out, 'w', buffering=1))
            files.append(reads_file)
            methods.append('num_reads')
        if denoise_mean_out is not None:
            denoise_mean_file = stack.enter_context(open(denoise_mean_out, 'w', buffering=1))
            files.append(denoise_mean_file)
            methods.append('denoise_mean')
        if denoise_deblur_out is not None:
            denoise_deblur_file = stack.enter_context(open(denoise_deblur_out, 'w', buffering=1))
            files.append(denoise_deblur_file)
            methods.append('denoise_deblur')
        if variants_out is not None:
            variants_file = stack.enter_context(open(variants_out, 'w', buffering=1))
            files.append(variants_file)
            methods.append('num_variants')
        if entropy_out is not None:
            entropy_file = stack.enter_context(open(entropy_out, 'w', buffering=1))
            files.append(entropy_file)
            methods.append('entropy')
        # write the headers for the tables
        for cfile in files:       
            cfile.write('uniref_id')
            for csample in all_samples:
                cfile.write(f'\t{csample}')
            cfile.write('\n')

        selected_ids = ids_list[:num_ids] if num_ids is not None else ids_list

        if jobs > 1:
            task_iter = [
                (
                    uniref_id,
                    base_dir,
                    methods,
                    window_size,
                    min_reads_for_norm,
                    mean_noise,
                    uniref_index[uniref_id],
                    min_samples_per_uniref,
                )
                for uniref_id in selected_ids
            ]
            executor = stack.enter_context(ProcessPoolExecutor(max_workers=jobs))
            result_iter = executor.map(process_uniref_task, task_iter)
        else:
            result_iter = (
                process_uniref_task((
                    uniref_id,
                    base_dir,
                    methods,
                    window_size,
                    min_reads_for_norm,
                    mean_noise,
                    uniref_index[uniref_id],
                    min_samples_per_uniref,
                ))
                for uniref_id in selected_ids
            )

        for cidx, (uniref_id, result, id_len) in enumerate(result_iter, start=1):
            logger.debug(f'Processing uniref_id {cidx}/{len(selected_ids)}: {uniref_id}')
            if id_len == 0:
                logger.warning(f'No valid results found for uniref_id {uniref_id}, skipping it.')
                continue
            for cfile in files:
                cfile.write(uniref_id)
            for csample in all_samples:
                csample_id = f'{csample}-splits'
                if variants_out is not None:
                    num_variants = result['num_variants'].get(csample_id, 0)
                    variants_file.write(f'\t{num_variants}')
                if num_reads_out is not None:
                    num_reads = result['num_reads'].get(csample_id, 0)
                    reads_file.write(f'\t{num_reads}')
                if entropy_out is not None:
                    entropy = result['entropy'].get(csample_id, 0)
                    entropy_file.write(f'\t{entropy}')
                if denoise_mean_out is not None:
                    num_denoise_mean = result['denoise_mean'].get(csample_id, 0)
                    denoise_mean_file.write(f'\t{num_denoise_mean}')
                if denoise_deblur_out is not None:
                    num_denoise_deblur = result['denoise_deblur'].get(csample_id, 0)
                    denoise_deblur_file.write(f'\t{num_denoise_deblur}')
            for cfile in files:
                cfile.write('\n')
            logger.info(f'Finished processing uniref_id {cidx}/{len(selected_ids)}: {uniref_id}')

def main(argv):
    parser = argparse.ArgumentParser(description='split-to-uniref version ' + __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('-b', '--base-dir', help='Base directory where the sample splits are located', default='.')
    parser.add_argument('-m', '--metadata-file', help='Metadata file (tab-separated) with sample information, must contain a "Run" column with sample ids', default='map.txt')
    parser.add_argument('-r', '--reads-file', help='File (tab-separated) with sample ids and their read counts, must contain "sampleid" and "reads" columns. tsv. If not provided, will auto-generate', default=None)
    parser.add_argument('--min-samples-per-uniref', type=int, help='Minimum number of samples a uniref_ID must be present in to be included in output table', default=5)
    parser.add_argument('--window-size', type=int, help='Window size for counting sequence variants', default=50)
    parser.add_argument('--log-level', type=str, help='Logging level (DEBUG, INFO, WARNING, ERROR)', default='INFO')
    parser.add_argument('--num-ids', type=int, help='Number of uniref ids to process (for testing purposes, None to process all)', default=None)
    parser.add_argument('--min-reads-for-norm', type=int, help='Minimum number of reads for normalization by reads (lower read number are updated to it), to avoid inflating the variant count for lowly covered genes (default: 20)', default=20)
    parser.add_argument('--num-reads-out', help='If provided, file name for the number of reads/unirefid table')
    parser.add_argument('--denoise-mean-out', help='If provided, file name for the num variants/unirefid after mean error rate subtraction (see --mean-noise)')
    parser.add_argument('--mean-noise', type=float, help='Mean per-nucleotide error rate for the mean-out denoising', default=0.005)
    parser.add_argument('--denoise-deblur-out', help='If provided, file name for the deblur style denoised num variants/unirefid table')
    parser.add_argument('--variants-out', help='If provided, file name for the non-denoised variants/unirefid table')
    parser.add_argument('--entropy-out', help='If provided, file name for the entropy/unirefid table')
    parser.add_argument('--jobs', type=int, help='Number of worker processes to use for per-UniRef processing', default=1)

    args = parser.parse_args(sys.argv[1:])
    logger.remove()  # Remove default logger
    logger.add(sys.stderr, level=args.log_level)  # Add new logger with specified log level

    logger.info("Starting split-to-uniref pipeline")
    split_to_uniref(base_dir=args.base_dir, metadata_file=args.metadata_file, reads_file=args.reads_file, min_samples_per_uniref=args.min_samples_per_uniref,
                    window_size=args.window_size, num_ids=args.num_ids, min_reads_for_norm=args.min_reads_for_norm,
                    num_reads_out=args.num_reads_out, denoise_mean_out=args.denoise_mean_out, mean_noise=args.mean_noise,
                    denoise_deblur_out=args.denoise_deblur_out, variants_out=args.variants_out, entropy_out=args.entropy_out,
                    jobs=args.jobs)
    logger.info("Split-to-uniref pipeline finished")


if __name__ == "__main__":
	main(sys.argv[1:])
