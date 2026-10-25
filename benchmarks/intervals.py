#!/usr/bin/env python3
"""Credible interval coverage (Phase 7a, item 7; M2, D7).

Three 95% intervals for each species' abundance, on the reference setting's
scores (exact database, -N 3, primary error profile, Emu score):

  mean-field  the Phase 5 columns: Beta(α_s, Σα − α_s) over the reported species
              (bayesian_vi with posterior_summary, as the pipeline's
              --credible_intervals; α summed over tax IDs that share a name)
  bootstrap   reads resampled with replacement (100 resamples, seeded), the
              pipeline's VB rerun on each; the 2.5% and 97.5% percentiles of the
              reported abundance (0 where a species is not reported)
  gibbs       a collapsed Gibbs sampler on the same sparse likelihood and
              Dirichlet(α0/S) prior, over every species with data (no pruning):
              z_r is resampled with p(z_r = s | z_-r) ∝ L_rs (n_s^-r + α0/S);
              after burn-in, each sweep draws π ~ Dirichlet(α0/S + n); the 2.5%
              and 97.5% percentiles of those draws, renormalised over the species
              the VB reports so that all three methods describe the same quantity

Coverage is the share of (sample, species) pairs whose interval contains the
truth. Truth is the design abundance renormalised over the reported species
(the model's π generates the reads); the realised read share is also recorded.
Writes <bench>/results/intervals.tsv (one row per sample, method and species)
and intervals_runtime.tsv.
"""
import argparse
import time
import zlib
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from numba import njit
from scipy.stats import beta as beta_dist

from common import BENCH, write_tsv
from evaluate import load_truth, taxonomy
from variational_inference_bayesian import bayesian_vi, sparse_likelihood_from_csr

N_BOOT = 100
GIBBS_BURN_IN = 200
GIBBS_SAMPLES = 1000
REF_NPZ = "lik/exact_N3_primary_emu.npz"
DB = BENCH / "db" / "exact"


def load_scores(sample):
    with np.load(BENCH / "runs" / sample / REF_NPZ) as d:
        offsets, tax_ids, scores = d["offsets"], d["tax_ids"], d["scores"]
    db_ids = np.array(sorted({int(t) for t in tax_ids}), dtype=np.int64)
    idx, logL = sparse_likelihood_from_csr(offsets, tax_ids, scores, db_ids)
    return idx, logL, db_ids


def run_vb(idx, logL, db_ids, summary=False):
    res = bayesian_vi(idx, logL, list(db_ids), abundance_threshold=1e-3, posterior_summary=summary)
    freq_full, freq_thresh = res[0], res[1]
    return (freq_thresh or freq_full), (res[4] if summary else None)


def by_name(values, tax):
    out = {}
    for tid, v in values.items():
        name = tax.loc[int(tid), "species"]
        out[name] = out.get(name, 0.0) + v
    return out


def mean_field(idx, logL, db_ids, tax):
    abund, post = run_vb(idx, logL, db_ids, summary=True)
    alpha, expected, shared = by_name({t: p["alpha"] for t, p in post.items()}, tax), \
        by_name({t: p["expected_reads"] for t, p in post.items()}, tax), \
        by_name({t: p["shared_reads"] for t, p in post.items()}, tax)
    A = sum(alpha.values())
    out = {}
    for name, a in alpha.items():
        rest = A - a
        lo, hi = beta_dist.ppf([0.025, 0.975], a, rest) if rest > 0 else (1.0, 1.0)
        out[name] = (float(lo), float(hi), shared[name] / expected[name] if expected[name] else 0.0)
    return by_name(abund, tax), out


def bootstrap(idx, logL, db_ids, tax, reported, seed):
    rng = np.random.default_rng(seed)
    R = idx.shape[0]
    draws = np.zeros((N_BOOT, len(reported)))
    col = {n: i for i, n in enumerate(reported)}
    for b in range(N_BOOT):
        rows = rng.integers(0, R, R)
        abund, _ = run_vb(idx[rows], logL[rows], db_ids)
        for name, v in by_name(abund, tax).items():
            if name in col:
                draws[b, col[name]] += v
    return {n: (float(np.percentile(draws[:, i], 2.5)), float(np.percentile(draws[:, i], 97.5)))
            for n, i in col.items()}


@njit(cache=True)
def _gibbs(idx, lik, n_cand, n_cols, prior, burn_in, n_samples, seed, z):
    np.random.seed(seed)
    R, K = idx.shape
    counts = np.zeros(n_cols)
    for r in range(R):
        counts[idx[r, z[r]]] += 1
    draws = np.zeros((n_samples, n_cols))
    w = np.zeros(K)
    for sweep in range(burn_in + n_samples):
        for r in range(R):
            if n_cand[r] == 1:  # one candidate: its assignment never changes
                continue
            counts[idx[r, z[r]]] -= 1
            total = 0.0
            for k in range(n_cand[r]):
                w[k] = lik[r, k] * (counts[idx[r, k]] + prior)
                total += w[k]
            u = np.random.random() * total
            k = 0
            acc = w[0]
            while acc < u and k < n_cand[r] - 1:
                k += 1
                acc += w[k]
            z[r] = k
            counts[idx[r, k]] += 1
        if sweep >= burn_in:
            total = 0.0
            for j in range(n_cols):  # π ~ Dirichlet(counts + prior), through gamma draws
                draws[sweep - burn_in, j] = np.random.standard_gamma(counts[j] + prior)
                total += draws[sweep - burn_in, j]
            draws[sweep - burn_in] /= total
    return draws


def gibbs(idx, logL, db_ids, tax, reported, seed, alpha_prior=1.0, renormalise=True):
    # Reads with at least one finite score; candidates sorted so the best is first
    keep = np.isfinite(logL).any(axis=1)
    idx, logL = idx[keep], logL[keep]
    order = np.argsort(-logL, axis=1)
    idx = np.take_along_axis(idx, order, axis=1)
    logL = np.take_along_axis(logL, order, axis=1)
    lik = np.exp(logL - logL[:, :1])  # 0 for padding
    cols = np.unique(idx[np.isfinite(logL)])
    remap = np.zeros(len(db_ids), dtype=np.int64)
    remap[cols] = np.arange(len(cols))
    idx_c = remap[idx]
    prior = alpha_prior / len(cols)
    n_cand = np.isfinite(logL).sum(axis=1)  # finite entries come first after the sort
    z = np.zeros(idx.shape[0], dtype=np.int64)  # start at each read's best candidate
    draws = _gibbs(idx_c, lik, n_cand, len(cols), prior, GIBBS_BURN_IN, GIBBS_SAMPLES, seed, z)
    names = tax.loc[db_ids[cols], "species"].to_numpy()
    col = {n: i for i, n in enumerate(reported)}
    sub = np.zeros((GIBBS_SAMPLES, len(reported)))
    for j, n in enumerate(names):
        if n in col:
            sub[:, col[n]] += draws[:, j]
    if renormalise:
        sub /= sub.sum(axis=1, keepdims=True)  # over the species the VB reports, as the other two
    return {n: (float(np.percentile(sub[:, i], 2.5)), float(np.percentile(sub[:, i], 97.5))) for n, i in col.items()}


def one_sample(sample):
    tax = taxonomy(str(DB))
    idx, logL, db_ids = load_scores(sample)
    truth = load_truth(BENCH / "sims" / sample)
    rows, times = [], {"sample": sample, "reads": idx.shape[0]}
    t0 = time.time()
    abund, mf = mean_field(idx, logL, db_ids, tax)
    times["mean-field"] = time.time() - t0
    reported = sorted(abund)
    t0 = time.time()
    bs = bootstrap(idx, logL, db_ids, tax, reported, seed=zlib.crc32(sample.encode()))
    times["bootstrap"] = time.time() - t0
    t0 = time.time()
    gb = gibbs(idx, logL, db_ids, tax, reported, seed=zlib.crc32(sample.encode()) % 2**31)
    times["gibbs"] = time.time() - t0
    design = dict(zip(truth["species"], truth["design_abundance"]))
    realised = dict(zip(truth["species"], truth["realised"]))
    d_total = sum(design.get(n, 0.0) for n in reported)
    r_total = sum(realised.get(n, 0.0) for n in reported)
    for method, ints in (("mean-field", {n: v[:2] for n, v in mf.items()}), ("bootstrap", bs), ("gibbs", gb)):
        for n in reported:
            lo, hi = ints[n]
            rows.append({"sample": sample, "method": method, "species": n, "estimate": abund[n],
                         "lower": lo, "upper": hi,
                         "truth_design": design.get(n, 0.0) / d_total if d_total else 0.0,
                         "truth_realised": realised.get(n, 0.0) / r_total if r_total else 0.0,
                         "shared_read_fraction": mf[n][2], "in_truth": n in design})
    return rows, times


def toy():
    """What the model can identify: species C (1,000 reads of its own) and two siblings A
    and B whose 300 shared reads score the same against both, B with 0, 5, 20 or 100
    reads of its own and A with none. Mean-field and Gibbs 95% intervals (Gibbs over
    every species, not renormalised). Writes docs/benchmarks/data/intervals_toy.tsv."""
    from common import DOCS_DATA
    rows = []
    for unique_b in (0, 5, 20, 100):
        R = 1000 + unique_b + 300
        idx = np.zeros((R, 3), dtype=np.int32)
        logL = np.full((R, 3), -np.inf)
        logL[:1000, 0] = -1.0
        idx[1000:1000 + unique_b, 0], logL[1000:1000 + unique_b, 0] = 2, -1.0
        idx[1000 + unique_b:, 0], idx[1000 + unique_b:, 1] = 1, 2
        logL[1000 + unique_b:, :2] = -1.0
        db_ids = np.array([10, 20, 30])
        tax = pd.DataFrame({"species": ["C", "A", "B"]}, index=db_ids)
        abund, mf = mean_field(idx, logL, db_ids, tax)
        gb = gibbs(idx, logL, db_ids, tax, ["A", "B", "C"], seed=1, renormalise=False)
        for sp in ("A", "B", "C"):
            rows.append({"unique_reads_B": unique_b, "species": sp, "vb_mean": abund.get(sp, 0.0),
                         "mf_lower": mf[sp][0] if sp in mf else float("nan"),
                         "mf_upper": mf[sp][1] if sp in mf else float("nan"),
                         "gibbs_lower": gb[sp][0], "gibbs_upper": gb[sp][1]})
    df = pd.DataFrame(rows)
    write_tsv(df, DOCS_DATA / "intervals_toy.tsv")
    print(df.round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--samples", help="comma-separated (default: every simulated sample with reference scores)")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--toy", action="store_true", help="only the identifiability example (intervals_toy.tsv)")
    args = ap.parse_args()
    if args.toy:
        return toy()
    samples = args.samples.split(",") if args.samples else sorted(
        p.name for p in (BENCH / "runs").iterdir() if (p / REF_NPZ).exists())
    rows, times = [], []
    with ProcessPoolExecutor(args.threads) as pool:
        for r, t in pool.map(one_sample, samples):
            rows += r
            times.append(t)
    df = pd.DataFrame(rows)
    for col in ("truth_design", "truth_realised"):
        df[f"covers_{col[6:]}"] = (df["lower"] <= df[col]) & (df[col] <= df["upper"])
    df["width"] = df["upper"] - df["lower"]
    write_tsv(df, BENCH / "results" / "intervals.tsv")
    write_tsv(pd.DataFrame(times), BENCH / "results" / "intervals_runtime.tsv")
    print(df.groupby("method")[["covers_design", "covers_realised", "width"]].mean())


if __name__ == "__main__":
    main()
