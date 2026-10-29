"""
The read-scores layout shared by LIKELIHOOD and RUN_VI (F24, F52).

likelihood.py writes each read's candidate species and scores as flat arrays
(the read-scores npz); variational_inference_bayesian.py builds the same arrays
from the dict form. Both use read_scores_to_arrays from this module, which
does not import pysam, so the estimators do not need it.
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
