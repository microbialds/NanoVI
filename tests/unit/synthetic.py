"""Synthetic per-read scores shared by the estimator tests."""
import numpy as np


def synthetic_reads(counts, seed=0, base=-40.0, max_decoys=2, mean_gap=3.0):
    """Reads with a known source species and lower-scoring decoy alignments.

    counts: {species_id: number of reads drawn from it}. Each read aligns to its
    source with score base + noise and to up to max_decoys other species with a
    score lower by an exponential gap.
    """
    rng = np.random.default_rng(seed)
    species = sorted(counts)
    log_p_rgs = {}
    n = 0
    for sid in species:
        for _ in range(counts[sid]):
            own = base + rng.normal(0.0, 2.0)
            others = [s for s in species if s != sid]
            k = min(int(rng.integers(0, max_decoys + 1)), len(others))
            decoys = list(rng.choice(others, size=k, replace=False)) if k else []
            sids = [sid] + [int(d) for d in decoys]
            scores = [own] + [own - float(rng.exponential(mean_gap)) for _ in decoys]
            log_p_rgs[f'read{n}'] = (sids, scores)
            n += 1
    return log_p_rgs
