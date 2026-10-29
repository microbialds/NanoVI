"""Shared helpers for the Phase 7a benchmark scripts.

The scripts import the pipeline's own modules from bin/ (likelihood,
variational_inference_bayesian, em_estimator, build_gtdb_db), so the benchmark
runs the code the pipeline runs. Data and results live outside the repository,
under NANOVI_BENCH (default ~/pipelines/NanoVI_benchmarks/phase7a, decision D11).
NANOVI_BIN points the scripts at another bin/, such as a checkout of commit 34d0411
for the Phase 7a arms whose options Phase 7b removed (D22, D25, D29, D31).
"""
import gzip
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = Path(os.environ.get("NANOVI_BIN", REPO / "bin")).resolve()
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))

BENCH = Path(os.environ.get("NANOVI_BENCH", Path.home() / "pipelines" / "NanoVI_benchmarks" / "phase7a"))
DOCS_DATA = REPO / "docs" / "benchmarks" / "data"
DOCS_FIG = REPO / "docs" / "benchmarks" / "figures"

GTDB = BENCH / "downloads" / "gtdb_r226"
SSU_ALL = GTDB / "ssu_all_r226.fna.gz"
SSU_REPS = GTDB / "ssu_reps_r226.fna.gz"

RANKS = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]

_RAW_SPECIES = re.compile(r"s__([^;\[]*)")


def raw_gtdb_species(header):
    """GTDB species name as in the header, suffixes kept ('Listeria monocytogenes_B')."""
    m = _RAW_SPECIES.search(header)
    return m.group(1).strip() if m else ""


def open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_taxonomy(db_dir):
    """tax_id -> lineage dict from a NanoVI taxonomy.tsv."""
    import pandas as pd
    tax = pd.read_csv(Path(db_dir) / "taxonomy.tsv", sep="\t", dtype=str, keep_default_na=False)
    tax["tax_id"] = tax["tax_id"].astype(int)
    return tax.set_index("tax_id")


def write_tsv(df, path, **kwargs):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=kwargs.pop("index", False), **kwargs)
    return path


_RAW_RANK = {rank: re.compile(prefix + r"([^;\[]*)") for prefix, rank in
             [("d__", "superkingdom"), ("p__", "phylum"), ("c__", "class"), ("o__", "order"),
              ("f__", "family"), ("g__", "genus"), ("s__", "species")]}


def raw_gtdb_lineage(header):
    """GTDB lineage as in the header, suffixes kept, rank -> name."""
    out = {}
    for rank, rx in _RAW_RANK.items():
        m = rx.search(header)
        out[rank] = m.group(1).strip() if m else ""
    return out
