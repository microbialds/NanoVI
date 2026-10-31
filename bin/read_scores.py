"""
The read-scores layout shared by LIKELIHOOD and RUN_VI (F24, F52).

likelihood.py writes each read's candidate species and scores as flat arrays
(the read-scores npz); variational_inference_bayesian.py builds the same arrays
from the dict form. Both use read_scores_to_arrays from this module, which
does not import pysam, so the estimators do not need it.

With --read-names, read_likelihood.py also stores the read names (F55):
"read_names", the name of each CSR row, and "unassigned_names", the reads with
no scored alignment, both in the order of the alignment file. Files without
them still load, without names.
"""

import numpy as np


def read_scores_to_arrays(log_p_rgs):
    """Per-read lists as flat arrays, in read order (the read-scores npz, F24).

    log_p_rgs maps each read to (species_ids, log_scores), as
    likelihood.score_reads returns it. Returns (offsets, tax_ids, scores):
    read r's candidates are tax_ids[offsets[r]:offsets[r + 1]], with their scores.
    """
    lengths = [len(sids) for sids, _ in log_p_rgs.values()]
    offsets = np.zeros(len(lengths) + 1, dtype=np.int64)
    np.cumsum(lengths, out=offsets[1:])
    tax_ids = np.fromiter((sid for sids, _ in log_p_rgs.values() for sid in sids),
                          dtype=np.int64, count=int(offsets[-1]))
    scores = np.fromiter((score for _, scores in log_p_rgs.values() for score in scores),
                         dtype=np.float64, count=int(offsets[-1]))
    return offsets, tax_ids, scores


def names_to_array(names):
    """Read names as one uint8 array: their UTF-8 bytes joined by newlines.

    A FASTQ read name holds no newline, and unlike a fixed-width bytes array one
    long name does not widen every other (np.load reads either without pickle).
    """
    return np.frombuffer("\n".join(names).encode(), dtype=np.uint8)


def names_from_array(array):
    """The read names that names_to_array stored."""
    text = array.tobytes().decode()
    return text.split("\n") if text else []


def save_read_scores(f, offsets, tax_ids, scores, read_names=None, unassigned_names=None):
    """Write the read-scores npz; read_names (one per row) and unassigned_names are optional."""
    names = {}
    if read_names is not None:
        if len(read_names) != len(offsets) - 1:
            raise ValueError(f"{len(read_names)} read names for {len(offsets) - 1} reads")
        names["read_names"] = names_to_array(read_names)
    if unassigned_names is not None:
        names["unassigned_names"] = names_to_array(unassigned_names)
    np.savez_compressed(f, offsets=offsets, tax_ids=tax_ids, scores=scores, **names)


def load_read_scores(path, names=True):
    """Read a read-scores npz: (offsets, tax_ids, scores, read_names, unassigned_names).

    The names are lists of str, or None for a file written without them; with
    names=False they are not read (None), which saves the memory they take.
    """
    with np.load(path) as data:
        offsets, tax_ids, scores = data["offsets"], data["tax_ids"], data["scores"]
        read_names, unassigned_names = (
            names_from_array(data[key]) if names and key in data.files else None
            for key in ("read_names", "unassigned_names")
        )
    return offsets, tax_ids, scores, read_names, unassigned_names
