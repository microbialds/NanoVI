"""Metrics against the truth of a simulated sample (Phase 7a, item 5).

An estimate is {species name: abundance} over the reported species (tax IDs
mapped to names through the database's taxonomy.tsv). The truth is the realised
read share of each community species (truth.tsv); the design abundances are kept
for the interval coverage test.

Metrics, at species and genus level, and with species merged into the groups that
share identical 16S sequences within a genus (db_audit.py; 'group' level):
precision, recall and F1 of the reported set against the species with at least
one read; L1 distance (sum of absolute differences, 0-2) and Bray-Curtis
dissimilarity (L1 / 2 for two compositions summing to 1) against the realised
shares; false-positive count and abundance (species reported but absent).
"""
import json
from functools import lru_cache
from pathlib import Path

import pandas as pd

from common import BENCH, read_taxonomy
from estimate_abundance import UNKNOWN_TAX_ID

UNKNOWN = "Unknown (not in database)"


@lru_cache(maxsize=None)
def taxonomy(db_dir):
    return read_taxonomy(db_dir)


@lru_cache(maxsize=None)
def identical_groups():
    """Species -> group label for species of one genus that share an identical 16S sequence."""
    path = BENCH / "audit" / "ssu_all" / "identical_groups.tsv"
    if not path.exists():
        return {}
    g = pd.read_csv(path, sep="\t")
    return dict(zip(g["species"], g["group"]))


def load_estimate(json_path, db_dir):
    """Reported abundances by species name (tax IDs sharing a name are summed)."""
    data = json.loads(Path(json_path).read_text())
    tax = taxonomy(str(db_dir))
    est = {}
    for tid, value in data["abundance"].items():
        tid = int(tid)
        name = UNKNOWN if tid == UNKNOWN_TAX_ID else tax.loc[tid, "species"]
        est[name] = est.get(name, 0.0) + value
    return est, data.get("below_detection_count", 0)


def load_truth(sample_dir):
    t = pd.read_csv(Path(sample_dir) / "truth.tsv", sep="\t")
    t["realised"] = t["reads"] / t["reads"].sum()
    return t


def _genus(name):
    return name.split(" ")[0] if name != UNKNOWN else UNKNOWN


def collapse(abund, level):
    if level == "species":
        return dict(abund)
    groups = identical_groups()
    key = _genus if level == "genus" else (lambda s: groups.get(s, s))
    out = {}
    for name, v in abund.items():
        out[key(name)] = out.get(key(name), 0.0) + v
    return out


def compare(est, truth_abund, level="species", exclude_unknown=True):
    """Metrics of one estimate against one truth ({name: share})."""
    if exclude_unknown:
        est = {k: v for k, v in est.items() if k != UNKNOWN}
        total = sum(est.values())
        est = {k: v / total for k, v in est.items()} if total > 0 else est
    e, t = collapse(est, level), collapse(truth_abund, level)
    present = {k for k, v in t.items() if v > 0}
    reported = {k for k, v in e.items() if v > 0}
    tp = len(present & reported)
    fp = reported - present
    precision = tp / len(reported) if reported else float("nan")
    recall = tp / len(present) if present else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if tp else 0.0
    keys = set(e) | set(t)
    l1 = sum(abs(e.get(k, 0.0) - t.get(k, 0.0)) for k in keys)
    return {
        f"{level}_precision": precision, f"{level}_recall": recall, f"{level}_f1": f1,
        f"{level}_l1": l1, f"{level}_bray_curtis": l1 / 2,
        f"{level}_fp_count": len(fp), f"{level}_fp_abundance": sum(e[k] for k in fp),
        f"{level}_tp": tp, f"{level}_reported": len(reported),
    }


def all_metrics(est, truth):
    t = dict(zip(truth["species"], truth["realised"]))
    out = {}
    for level in ("species", "genus", "group"):
        out.update(compare(est, t, level))
    out["unknown_abundance"] = est.get(UNKNOWN, 0.0)
    return out


def per_species(est, truth):
    """One row per truth or reported species: estimate, realised share, design abundance."""
    rows = []
    t = truth.set_index("species")
    total = sum(v for k, v in est.items() if k != UNKNOWN)
    for name in sorted(set(t.index) | set(est)):
        rows.append({
            "species": name,
            "estimate": est.get(name, 0.0),
            "estimate_known": est.get(name, 0.0) / total if name != UNKNOWN and total else float("nan"),
            "realised": float(t.loc[name, "realised"]) if name in t.index else 0.0,
            "design": float(t.loc[name, "design_abundance"]) if name in t.index else 0.0,
            "reads": int(t.loc[name, "reads"]) if name in t.index else 0,
            "in_truth": name in t.index,
        })
    return rows


# ─── Real mock communities ──────────────────────────────────────────────────

TRUTH_DIR = Path(__file__).resolve().parent / "truth"
REAL_TRUTH = {
    # sample prefix -> truth tables (abundances, and species expected without one)
    "zymo_MSPlus": ["zymo_d6305_16s.tsv", "zymo_msplus_spikes.tsv"],
    "zymo_GutMock": ["zymo_d6331_16s.tsv"],
    "emu_zymo": ["zymo_d6305_16s.tsv"],
    "preprint_zymo": ["zymo_d6305_16s.tsv"],
    "emu_gutmock": ["emu_gutmock_ont.tsv"],
}


def load_real_truth(sample):
    """{GTDB species: expected 16S share, or nan when only presence is known}."""
    for prefix, files in REAL_TRUTH.items():
        if sample.startswith(prefix):
            t = pd.concat([pd.read_csv(TRUTH_DIR / f, sep="\t") for f in files])
            return dict(zip(t["gtdb_species"], t["abundance"]))
    raise KeyError(sample)


def compare_real(est, truth, n_reads, level="species", min_expected_reads=10):
    """Metrics of an estimate on a real mock.

    Expected species: every species in the truth tables. Recall counts those whose
    expected reads (share x n_reads) are at least min_expected_reads, plus the
    spike-ins (presence only); precision and false positives count every reported
    species outside the expected set. L1 and Bray-Curtis use the species with a known
    share, the estimate renormalised over them.
    """
    est = {k: v for k, v in est.items() if k != UNKNOWN and v > 0}
    e = collapse(est, level)
    known = {k: v for k, v in truth.items() if v == v}
    expected = collapse({k: (v if v == v else 1.0) for k, v in truth.items()}, level)
    detectable = collapse({k: (v if v == v else 1.0) for k, v in truth.items()
                           if v != v or v * n_reads >= min_expected_reads}, level)
    reported = set(e)
    fp = reported - set(expected)
    tp = len(reported & set(detectable))
    total = sum(e.values())
    t_known = collapse(known, level)
    e_known = {k: v for k, v in e.items() if k in t_known}
    s = sum(e_known.values())
    tk = sum(t_known.values())
    l1 = sum(abs(e_known.get(k, 0.0) / s - t_known[k] / tk) for k in t_known) if s else float("nan")
    precision = len(reported & set(expected)) / len(reported) if reported else float("nan")
    recall = tp / len(detectable) if detectable else float("nan")
    return {
        f"{level}_precision": precision, f"{level}_recall": recall,
        f"{level}_f1": 2 * precision * recall / (precision + recall) if precision and recall else 0.0,
        f"{level}_l1_known": l1, f"{level}_bray_curtis_known": l1 / 2 if l1 == l1 else l1,
        f"{level}_fp_count": len(fp), f"{level}_fp_abundance": sum(e[k] for k in fp) / total if total else 0.0,
        f"{level}_reported": len(reported),
    }


def real_metrics(est, sample, n_reads):
    truth = load_real_truth(sample)
    out = {}
    for level in ("species", "genus", "group"):
        out.update(compare_real(est, truth, n_reads, level))
    out["unknown_abundance"] = est.get(UNKNOWN, 0.0)
    return out
