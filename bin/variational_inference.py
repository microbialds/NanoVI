# inference.py

import math
import numpy as np
import pysam
from operator import add, mul
from multiprocessing import Pool
from sys import stdout
from flatten_dict import unflatten
from Bio.Seq import Seq

CIGAR_OPS = [1, 2, 4, 10]
CIGAR_OPS_ALL = [0, 1, 2, 4]

def get_align_stats(alignment):
    cigar_stats = alignment.get_cigar_stats()[0]
    n_mismatch = cigar_stats[10] - cigar_stats[1] - cigar_stats[2]
    return [cigar_stats[1], cigar_stats[2], cigar_stats[4], n_mismatch]
    pass

def get_align_len(alignment):
    return sum(alignment.get_cigar_stats()[0][cigar_op] for cigar_op in CIGAR_OPS_ALL)
    pass

def extract_alignment_info(alignment):
    cigar_stats = get_align_stats(alignment)
    align_len = get_align_len(alignment)
    query_name = alignment.query_name
    return query_name, align_len, cigar_stats
    pass

def process_alignment_chunk(alignment_chunk):
    cigar_stats_primary = [0] * len(CIGAR_OPS)
    dict_longest_align = {}
    
    for alignment_info in alignment_chunk:
        query_name, align_len, cigar_stats = alignment_info
        if align_len not in dict_longest_align:
            dict_longest_align[query_name] = align_len
        cigar_stats_primary = list(map(add, cigar_stats_primary, cigar_stats))
        if dict_longest_align[query_name] < align_len:
            dict_longest_align[query_name] = align_len
    
    return cigar_stats_primary, dict_longest_align
    pass

def get_cigar_op_log_probabilities(sam_path, threads):
    cigar_stats_primary = [0] * len(CIGAR_OPS) 
    dict_longest_align = {}

    sam_pysam = pysam.AlignmentFile(sam_path)
    
    # Extract necessary information from alignments
    alignment_info = [extract_alignment_info(alignment) for alignment in sam_pysam.fetch()]

    # Split alignment info into chunks for parallel processing
    chunk_size = len(alignment_info) // threads
    alignment_chunks = [alignment_info[i:i + chunk_size] for i in range(0, len(alignment_info), chunk_size)]
    
    with Pool(threads) as pool:
        results = pool.map(process_alignment_chunk, alignment_chunks)

    # Combine results from all processes
    for result in results:
        stats, longest_align = result
        cigar_stats_primary = list(map(add, cigar_stats_primary, stats))
        dict_longest_align.update(longest_align)

    zero_locs = [i for i, e in enumerate(cigar_stats_primary) if e == 0]
    if zero_locs:
        for i in sorted(zero_locs, reverse=True):
            del cigar_stats_primary[i]
    
    n_char = sum(cigar_stats_primary)
    return [math.log(x) for x in np.array(cigar_stats_primary) / n_char], zero_locs, dict_longest_align
    pass

def compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op, dict_longest_align, align_len):
    ref_name, query_name = alignment.reference_name, alignment.query_name
    log_score = sum(list(map(mul, log_p_cigar_op, cigar_stats))) * (dict_longest_align[query_name] / align_len)
    species_tid = int(ref_name.split(":")[0])
    return log_score, query_name, species_tid
    pass

def log_prob_rgs_dict(sam_path, log_p_cigar_op, dict_longest_align, p_cigar_op_zero_locs=None):
    log_p_rgs, unassigned_set = {}, set()
    sam_filename = pysam.AlignmentFile(sam_path, 'rb')

    if not p_cigar_op_zero_locs:
        for alignment in sam_filename.fetch():
            align_len = get_align_len(alignment)
            if alignment.reference_name and align_len:
                cigar_stats = get_align_stats(alignment)
                log_score, query_name, species_tid = compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op, dict_longest_align, align_len)
                if query_name not in log_p_rgs:
                    log_p_rgs[query_name] = ([species_tid], [log_score])
                elif query_name in log_p_rgs:
                    if species_tid not in log_p_rgs[query_name][0]:
                        log_p_rgs[query_name] = (log_p_rgs[query_name][0] + [species_tid], log_p_rgs[query_name][1] + [log_score])
                    else:
                        logprgs_idx = log_p_rgs[query_name][0].index(species_tid)
                        if log_p_rgs[query_name][1][logprgs_idx] < log_score:
                            log_p_rgs[query_name][1][logprgs_idx] = log_score
            else:
                unassigned_set.add(alignment.query_name)
    else:
        for alignment in sam_filename.fetch():
            align_len = get_align_len(alignment)
            if alignment.reference_name and align_len:
                cigar_stats = get_align_stats(alignment)
                if sum(cigar_stats[x] for x in p_cigar_op_zero_locs) == 0:
                    for i in sorted(p_cigar_op_zero_locs, reverse=True):
                        del cigar_stats[i]
                    log_score, query_name, species_tid = compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op, dict_longest_align, align_len)
                    if query_name not in log_p_rgs:
                        log_p_rgs[query_name] = ([species_tid], [log_score])
                    elif query_name in log_p_rgs and species_tid not in log_p_rgs[query_name][0]:
                        log_p_rgs[query_name] = (log_p_rgs[query_name][0] + [species_tid], log_p_rgs[query_name][1] + [log_score])
                    else:
                        logprgs_idx = log_p_rgs[query_name][0].index(species_tid)
                        if log_p_rgs[query_name][1][logprgs_idx] < log_score:
                            log_p_rgs[query_name][1][logprgs_idx] = log_score
            else:
                unassigned_set.add(alignment.query_name)

    assigned_reads = set(log_p_rgs.keys())
    unassigned_reads = unassigned_set - assigned_reads
    unassigned_count = len(unassigned_reads)
    stdout.write(f"Unassigned read count: {unassigned_count}\n")

    return log_p_rgs, unassigned_count, len(assigned_reads)
    pass

def variational_inference(log_p_rgs, freq, max_iterations=20, tolerance=5e-4):
    for iteration in range(max_iterations):
        q_dist, elbo = update_q_distribution(log_p_rgs, freq)
        new_freq = maximize_elbo(q_dist)
        for key in freq:
            if key not in new_freq:
                new_freq[key] = 0.0
        diff = sum(abs(new_freq[key] - freq.get(key, 0.0)) for key in new_freq)
        if diff < tolerance:
            print(f"Converged at iteration {iteration + 1} with ELBO = {elbo}")
            break
        freq = new_freq
        print(f"Iteration {iteration + 1}: ELBO = {elbo}")
    return freq, elbo
    pass

def update_q_distribution(log_p_rgs, freq):
    q_dist = {}
    elbo = 0
    for read in log_p_rgs:
        valid_seqs, log_p_rns = [], []
        for seq in range(len(log_p_rgs[read][0])):
            s_val = log_p_rgs[read][0][seq]
            if s_val in freq and freq[s_val] != 0:
                logprns_val = log_p_rgs[read][1][seq] + math.log(freq[s_val])
                valid_seqs.append(s_val)
                log_p_rns.append(logprns_val)
        if len(valid_seqs) != 0:
            logc = -np.max(log_p_rns)
            prnsc = np.exp(log_p_rns + logc)
            prc = np.sum(prnsc)
            elbo += (np.log(prc) - logc)
            for seq in enumerate(valid_seqs):
                q_dist[(seq[1], read)] = prnsc[seq[0]] / prc
    return unflatten(q_dist), elbo
    pass

def maximize_elbo(q_dist):
    freq = {}
    for tax_id, read_id in q_dist.items():
        freq[tax_id] = sum(read_id.values())
    total = sum(freq.values())
    for tax_id in freq:
        freq[tax_id] /= total
    return freq
    pass

def variational_inference_iterations(log_p_rgs, db_ids, lli_thresh, input_threshold):
    n_db = len(db_ids)
    n_reads = len(log_p_rgs)
    stdout.write(f"Assigned read count: {n_reads}\n")
    if n_reads == 0:
        raise ValueError("0 reads assigned")
    freq, counter = dict.fromkeys(db_ids, 1 / n_db), 1
    freq_thresh = 1 / n_reads
    if n_reads > 1000:
        freq_thresh = 10 / n_reads
    total_elbo = -math.inf
    while True:
        freq, _ = variational_inference(log_p_rgs, freq)
        _, new_elbo = update_q_distribution(log_p_rgs, freq)
        elbo_diff = new_elbo - total_elbo
        total_elbo = new_elbo
        if elbo_diff < 0:
            raise ValueError("ELBO decreased from prior iteration")
        if elbo_diff < lli_thresh:
            stdout.write(f"Number of VI iterations: {counter}\n")
            freq = {k: v for k, v in freq.items() if v >= freq_thresh}
            freq_full, _ = variational_inference(log_p_rgs, freq)
            freq_set_thresh = None
            if freq_thresh < input_threshold:
                freq = {k: v for k, v in freq_full.items() if v >= input_threshold}
                freq_set_thresh, _ = variational_inference(log_p_rgs, freq)
            return freq_full, freq_set_thresh, None
        counter += 1

    pass

def output_read_assignments(p_sgr, tsv_output_path):
    dist_df = pd.DataFrame(p_sgr)
    dist_df.to_csv("{}.tsv".format(tsv_output_path), sep='\t')
    return dist_df
    pass
