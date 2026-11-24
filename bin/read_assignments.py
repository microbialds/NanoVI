"""
Per-read species assignments (F55, D37).

estimate_abundance.py --read-assignments writes, for every read of the sample,
the table rows whose probability for the read is at least a threshold:

    read_id    species    probability

The probability of a row is the read's responsibility summed over the row's tax
IDs, from the final fit (φ for VB, EM's posterior for --estimator em), so the
rows are those of the abundance table (taxonomy.species_row_keys). A read's
probability on species that are not reported (all its species pruned, or
dropped by --min-abundance) is its 'Below detection' row; reads with no scored
alignment are 'Unassigned' with probability 1. Rows below the threshold are
left out, so a read's listed probabilities sum to 1 less what they leave out.
Two table rows with the same species name and different lineages (possible in
a custom taxonomy) are listed under that one name, as the table's species
column shows them.
"""

import gzip

import numpy as np

from taxonomy import BELOW_DETECTION_LABEL, UNASSIGNED_LABEL
from variational_inference_bayesian import SHARED_READ_MIN_SPLIT

# The share of a read that shared_read_fraction counts as a split (D37)
DEFAULT_MIN_PROBABILITY = SHARED_READ_MIN_SPLIT

# Lines gzipped at a time, so a large sample's file is not built in memory at once
CHUNK_LINES = 100_000


def read_assignment_rows(idx, phi, row_of_column, n_rows, min_probability=DEFAULT_MIN_PROBABILITY):
    """The (read, row, probability) lines of the scored reads, in file order.

    idx, phi: the padded responsibilities, idx the entry's column among the
    reported species or -1 (padding, or a species not reported). row_of_column:
    the table row (0 .. n_rows - 1) of each reported column. Row n_rows is
    'Below detection': 1 minus the read's responsibility on reported species,
    exactly 1 only when that responsibility is 0 (the table's Below detection
    reads), so a read with a tiny share on a reported species keeps it.
    Lines with a probability of 0 or below min_probability are left out. Reads
    keep their order; within a read, rows go by decreasing probability, then
    Below detection.
    """
    R = idx.shape[0]
    on = idx >= 0
    reads, cols = np.nonzero(on)
    weights = phi[reads, cols]
    rows = np.asarray(row_of_column, dtype=np.int64)[idx[reads, cols]]
    # A read's responsibility summed over the tax IDs of each row
    keys, inverse = np.unique(reads.astype(np.int64) * (n_rows + 1) + rows, return_inverse=True)
    prob = np.bincount(inverse, weights=weights, minlength=len(keys))
    line_reads, line_rows = keys // (n_rows + 1), keys % (n_rows + 1)

    reported = np.bincount(reads, weights=weights, minlength=R)
    below = np.where(reported > 0, np.minimum(1.0 - reported, np.nextafter(1.0, 0.0)), 1.0)
    line_reads = np.concatenate((line_reads, np.arange(R)))
    line_rows = np.concatenate((line_rows, np.full(R, n_rows)))
    prob = np.concatenate((prob, np.maximum(below, 0.0)))

    keep = (prob > 0) & (prob >= min_probability)
    line_reads, line_rows, prob = line_reads[keep], line_rows[keep], prob[keep]
    order = np.lexsort((line_rows, -prob, line_rows == n_rows, line_reads))
    return line_reads[order], line_rows[order], prob[order]


def write_read_assignments(path, read_names, unassigned_names, idx, phi, reported_ids, row_keys,
                           min_probability=DEFAULT_MIN_PROBABILITY):
    """Write <sample>_read_assignments.tsv.gz.

    read_names: one per row of idx/phi (the read-scores npz order);
    unassigned_names: the reads with no scored alignment, listed last as
    Unassigned. reported_ids: the species of the reported columns; row_keys
    maps each of them to its table row (taxonomy.species_row_keys), whose first
    element is the species name the table shows. Probabilities are written at
    full precision. The gzip header carries no name or time, so the same
    assignments give the same bytes.
    """
    row_number, labels = {}, []
    row_of_column = []
    for sid in reported_ids:
        key = row_keys[sid]
        if key not in row_number:
            row_number[key] = len(labels)
            labels.append(key[0])
        row_of_column.append(row_number[key])
    labels.append(BELOW_DETECTION_LABEL)

    line_reads, line_rows, prob = read_assignment_rows(idx, phi, row_of_column, len(labels) - 1, min_probability)
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as f:
        f.write(b"read_id\tspecies\tprobability\n")
        for start in range(0, len(prob), CHUNK_LINES):
            end = start + CHUNK_LINES
            f.write("".join(
                f"{read_names[r]}\t{labels[g]}\t{p!r}\n"
                for r, g, p in zip(line_reads[start:end].tolist(), line_rows[start:end].tolist(),
                                   prob[start:end].tolist())
            ).encode())
        for start in range(0, len(unassigned_names), CHUNK_LINES):
            f.write("".join(f"{name}\t{UNASSIGNED_LABEL}\t1.0\n"
                            for name in unassigned_names[start:start + CHUNK_LINES]).encode())
