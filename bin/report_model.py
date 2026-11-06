"""The model of the run report (--report; F60, D42): every number, caption and table that the PDF
(report_render.write_pdf, Typst) and the HTML site (report_render.write_site, Jinja) show, built once so the two
differ only in their short prose. The layout follows BugBuster's report model (gene2dis/BugBuster dev @ de57f34,
bin/summary_report.py): numbered figures and tables, the HTML showing up to `html_rows` rows of a table and the PDF
at most PDF_MAX_ROWS rows and PDF_MAX_COLS columns, with a pointer to the full table.

The report reads only what a run publishes, laid out as its output directory:
  sample_outputs/<sample>_rel-abundance.tsv      every sample (required)
  sample_outputs/<sample>_groups.tsv             --ambiguity_groups (D39)
  sample_outputs/<sample>_cell-abundance.tsv     --cell_abundance (D40); two '#' lines first
  sample_outputs/<sample>_read_assignments.tsv.gz  --read_assignments (D37)
  consolidated_output/otu_table_abundance.tsv, tax_table.tsv, sample_metadata.tsv
  qc/fastplong/<sample>_report.json.gz           fastplong's report
  pipeline_info/software_versions.yml
A missing optional table leaves its section out.
"""

import csv
import gzip
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from taxonomy import BELOW_DETECTION_LABEL as BELOW, UNASSIGNED_LABEL as UNASSIGNED

PLACEHOLDERS = (BELOW, UNASSIGNED)
FLAG_SHARED = 0.5      # shared_read_fraction at which a species is flagged in a sample (D42)
TOP_PER_SAMPLE = 5     # composition: the union of each sample's top 5 taxa ...
MAX_TAXA = 20          # ... capped at 20 by mean abundance; the rest is Other
HEATMAP_ABOVE = 14     # a heatmap of the top species once there are more samples than this
HEATMAP_TAXA = 30
PART_MAX = 60          # static figures with a sample axis are split into parts of at most 60 samples
PDF_MAX_ROWS = 60
PDF_MAX_COLS = 14
PORTRAIT_MM = 183      # text width of a portrait page; a wider PDF table goes on a landscape page
TOP_PROB_CLASSES = [("p >= 0.99", 0.99, 1.0 + 1e-12), ("0.5 <= p < 0.99", 0.5, 0.99), ("p < 0.5", 0.0, 0.5)]

# (id, title, HTML page)
SECTIONS = [
    ("overview", "Overview", "index.html"),
    ("composition", "Composition", "composition.html"),
    ("flagged", "Flagged species", "flagged.html"),
    ("groups", "Ambiguity groups", "groups.html"),
    ("cell", "Cell-level abundance", "cell.html"),
    ("reads", "Per-read assignments", "reads.html"),
    ("methods", "Methods", "methods.html"),
]
PAGE = {sid: page for sid, _, page in SECTIONS}


def warn(msg):
    print(f"WARNING: {msg}", file=sys.stderr)


def die(msg):
    sys.exit(f"ERROR: {msg}")


def natural_key(text):
    """V2 before V10: digits compare as numbers."""
    return [(0, int(t), "") if t.isdigit() else (1, 0, t.lower()) for t in re.split(r"(\d+)", str(text)) if t]


# ---- reading the results directory

def _read_tsv(path, comment=None):
    return pd.read_csv(path, sep="\t", comment=comment, keep_default_na=False, na_values=[""], encoding="utf-8")


def _text(col):
    """A column of names as text, blank where empty (an all-empty column reads as float NaN)."""
    return ["" if pd.isna(v) else str(v) for v in col]


def read_cell_table(path):
    """A <sample>_cell-abundance.tsv: its '#' lines (what it leaves out, D40) and the table."""
    with open(path, encoding="utf-8") as fh:
        comments = [ln[1:].strip() for ln in fh if ln.startswith("#")]
    return comments, _read_tsv(path, comment="#")


def read_fastplong(path):
    """fastplong's JSON report (gzipped as the pipeline publishes it, or plain): read counts before and after."""
    try:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as fh:
            d = json.load(fh)
    except (OSError, ValueError) as e:
        warn(f"{path}: not a fastplong JSON report ({e}); its read counts are left out")
        return None
    s, f = d.get("summary", {}), d.get("filtering_result", {})
    out = {"raw_reads": s.get("before_filtering", {}).get("total_reads"),
           "passed_reads": s.get("after_filtering", {}).get("total_reads"),
           "fastplong_version": s.get("fastplong_version")}
    for k in ("too_short_reads", "too_long_reads", "low_quality_reads", "too_many_N_reads", "chimeric_reads"):
        out[k] = f.get(k)
    return out


def read_metadata(path):
    """User sample metadata (--report_metadata): TSV or CSV with sample_id and, optionally, group."""
    sep = "," if str(path).lower().endswith(".csv") else "\t"
    m = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, encoding="utf-8")
    m.columns = [c.strip() for c in m.columns]
    if "sample_id" not in m.columns:
        die(f"{path}: the report metadata needs a sample_id column")
    if m["sample_id"].duplicated().any():
        die(f"{path}: sample_id repeated: {', '.join(sorted(set(m.loc[m['sample_id'].duplicated(), 'sample_id'])))}")
    return m


def load_results(results, metadata=None):
    """Everything the report reads from a results directory, as plain data."""
    results = Path(results)
    so = results / "sample_outputs"
    tables = sorted(so.glob("*_rel-abundance.tsv"))
    if not tables:
        die(f"no per-sample tables (*_rel-abundance.tsv) in {so}")
    run = {"results": str(results), "species": {}, "groups": {}, "cell": {}, "cell_comments": {}, "reads": {},
           "fastplong": {}}
    for t in tables:
        s = t.name[: -len("_rel-abundance.tsv")]
        df = _read_tsv(t)
        for c in ("identical_16s_species",):
            if c in df.columns:
                df[c] = _text(df[c])
        run["species"][s] = df
        g = so / f"{s}_groups.tsv"
        if g.exists():
            gd = _read_tsv(g)
            for c in ("members", "identical_16s_relatives", "lowest_common_rank", "lowest_common_taxon"):
                if c in gd.columns:
                    gd[c] = _text(gd[c])
            run["groups"][s] = gd
        c = so / f"{s}_cell-abundance.tsv"
        if c.exists():
            run["cell_comments"][s], run["cell"][s] = read_cell_table(c)
        r = so / f"{s}_read_assignments.tsv.gz"
        if r.exists():
            run["reads"][s] = r
        q = results / "qc" / "fastplong" / f"{s}_report.json.gz"
        q = q if q.exists() else results / "qc" / "fastplong" / f"{s}_report.json"
        run["fastplong"][s] = read_fastplong(q) if q.exists() else None
    co = results / "consolidated_output"
    for name in ("otu_table_abundance", "tax_table", "sample_metadata"):
        p = co / f"{name}.tsv"
        run[name] = _read_tsv(p) if p.exists() else None
    if run["otu_table_abundance"] is None or run["tax_table"] is None:
        die(f"{co} needs otu_table_abundance.tsv and tax_table.tsv")
    sv = results / "pipeline_info" / "software_versions.yml"
    run["software_versions"] = load_yaml(sv) if sv.exists() else {}
    if not sv.exists():
        warn(f"no {sv}: the methods list only the report's own software")
    samples = list(run["species"])
    groups = {s: "" for s in samples}
    run["metadata"] = None
    if metadata:
        m = read_metadata(metadata)
        missing = sorted(set(m["sample_id"]) - set(samples))
        if missing:
            warn(f"report metadata: {len(missing)} samples not in this run ({', '.join(missing[:5])})")
        m = m[m["sample_id"].isin(samples)]
        if "group" in m.columns:
            groups.update(dict(zip(m["sample_id"], m["group"])))
        run["metadata"] = m
    run["group"] = groups
    run["samples"] = order_samples(samples, groups, metadata_order=None if run["metadata"] is None else
                                   list(run["metadata"].get("group", pd.Series(dtype=str))))
    return run


def load_yaml(path):
    import yaml
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def order_samples(samples, groups, metadata_order=None):
    """Samples by group (in the order the metadata first names each group, samples without one last), then in
    natural order."""
    rank = {}
    for g in metadata_order or []:
        if g != "" and g not in rank:
            rank[g] = len(rank)
    return sorted(samples, key=lambda s: (groups.get(s, "") == "", rank.get(groups.get(s, ""), len(rank)),
                                          natural_key(s)))


def split_parts(samples, size=PART_MAX):
    """Consecutive parts of equal size (differing by at most one), none larger than `size`: 117 -> 59 + 58."""
    n = len(samples)
    k = max(1, math.ceil(n / size))
    base, extra = divmod(n, k)
    out, start = [], 0
    for i in range(k):
        m = base + (1 if i < extra else 0)
        out.append(samples[start:start + m])
        start += m
    return out


# ---- the numbers

def species_rows(df):
    return df[~df["species"].isin(PLACEHOLDERS)]


def read_fate(run):
    """One row per sample: fastplong's raw and passed reads and where the passed reads went (reported species,
    Below detection, Unassigned; from the table's estimated_counts, so empty with --keep_counts false)."""
    rows = []
    for s in run["samples"]:
        df, q = run["species"][s], run["fastplong"].get(s) or {}
        has_counts = "estimated_counts" in df.columns
        cnt = (lambda mask: int(df.loc[mask, "estimated_counts"].sum())) if has_counts else (lambda mask: None)
        r = {"sample": s, "group": run["group"].get(s, ""),
             "raw_reads": q.get("raw_reads"), "passed_reads": q.get("passed_reads"),
             "too_short": q.get("too_short_reads"), "too_long": q.get("too_long_reads"),
             "low_quality": q.get("low_quality_reads"), "chimeric": q.get("chimeric_reads"),
             "species_reads": cnt(~df["species"].isin(PLACEHOLDERS)), "below_detection_reads": cnt(df["species"] == BELOW),
             "unassigned_reads": cnt(df["species"] == UNASSIGNED),
             "n_species": int((species_rows(df)["abundance"] > 0).sum())}
        total = None
        if has_counts:
            total = r["species_reads"] + r["below_detection_reads"] + r["unassigned_reads"]
            if r["passed_reads"] is not None and total != r["passed_reads"]:
                warn(f"{s}: the table counts {total} reads, fastplong passed {r['passed_reads']}")
        for k, v in (("pct_species", "species_reads"), ("pct_below_detection", "below_detection_reads"),
                     ("pct_unassigned", "unassigned_reads")):
            r[k] = 100.0 * r[v] / total if total else (0.0 if total == 0 else None)
        rows.append(r)
    out = pd.DataFrame(rows)
    if not any(run["group"].values()):
        out = out.drop(columns="group")
    return out


def rank_abundance(run, rank):
    """taxon x sample read-level abundance at `rank` from the consolidated tables (species rows only; a genus is
    the sum of its species), samples in report order."""
    otu = run["otu_table_abundance"].set_index("species")
    otu = otu.loc[~otu.index.isin(PLACEHOLDERS)]
    otu = otu.reindex(columns=run["samples"]).fillna(0.0).astype(float)
    if rank == "species":
        return otu
    tax = run["tax_table"].set_index("species")
    key = [t if t else f"(no {rank})" for t in _text(tax[rank].reindex(otu.index))]
    return otu.groupby(key).sum()


def top_taxa(wide, per_sample=TOP_PER_SAMPLE, cap=MAX_TAXA):
    """The union of each sample's top `per_sample` taxa (abundance > 0; ties by name), capped at `cap` by mean
    abundance over the samples; ordered by mean abundance, then name."""
    union = set()
    for s in wide.columns:
        col = wide[s][wide[s] > 0]
        union.update(sorted(col.index, key=lambda t: (-col[t], str(t)))[:per_sample])
    mean = wide.mean(axis=1)
    return sorted(union, key=lambda t: (-mean[t], str(t)))[:cap]


def composition(run, rank):
    """Long table (sample, taxon, abundance) of the taxa shown and Other (the sum of every other taxon)."""
    wide = rank_abundance(run, rank)
    shown = top_taxa(wide)
    rows = []
    for s in run["samples"]:
        col = wide[s]
        for t in shown:
            rows.append({"sample": s, "taxon": t, "abundance": float(col[t])})
        rows.append({"sample": s, "taxon": "Other", "abundance": float(col.drop(index=shown).sum())})
    return pd.DataFrame(rows), shown


def heatmap_taxa(run, rank="species", n=HEATMAP_TAXA):
    wide = rank_abundance(run, rank)
    mean = wide.mean(axis=1)
    keep = sorted([t for t in wide.index if mean[t] > 0], key=lambda t: (-mean[t], str(t)))[:n]
    return wide.loc[keep]


def flagged_long(run):
    """Every reported species in every sample with its flags: shared_read_fraction >= FLAG_SHARED, and a
    non-empty identical_16s_species (D24, D27, D42)."""
    rows = []
    for s in run["samples"]:
        df = species_rows(run["species"][s])
        for r in df.itertuples(index=False):
            shared = getattr(r, "shared_read_fraction", np.nan)
            ident = getattr(r, "identical_16s_species", "") or ""
            rows.append({"sample": s, "species": r.species, "abundance": float(r.abundance),
                         "shared_read_fraction": float(shared) if shared == shared else np.nan,
                         "identical_16s_species": ident,
                         "flag_shared": bool(shared == shared and shared >= FLAG_SHARED),
                         "flag_identical": ident != ""})
    out = pd.DataFrame(rows, columns=["sample", "species", "abundance", "shared_read_fraction",
                                      "identical_16s_species", "flag_shared", "flag_identical"])
    return out


def flagged_summary(long):
    """One row per species flagged in at least one sample."""
    if long.empty:
        return long
    f = long[long["flag_shared"] | long["flag_identical"]]
    rows = []
    for sp, g in f.groupby("species", sort=False):
        allr = long[long["species"] == sp]
        idents = sorted({x for v in allr["identical_16s_species"] for x in v.split("; ") if x})
        rows.append({"species": sp, "samples_detected": int(len(allr)),
                     "samples_shared": int(allr["flag_shared"].sum()),
                     "samples_identical_16s": int(allr["flag_identical"].sum()),
                     "median_abundance": float(allr["abundance"].median()),
                     "max_abundance": float(allr["abundance"].max()),
                     "max_shared_read_fraction": float(allr["shared_read_fraction"].max()),
                     "identical_16s_species": "; ".join(idents)})
    out = pd.DataFrame(rows)
    return out.sort_values(["max_abundance", "species"], ascending=[False, True]).reset_index(drop=True)


def groups_long(run):
    """Every sample's groups of two or more members (D39)."""
    frames = []
    for s in run["samples"]:
        g = run["groups"].get(s)
        if g is None:
            continue
        g = g[~g["group"].isin(PLACEHOLDERS)]
        g = g[pd.to_numeric(g["n_members"], errors="coerce") >= 2]
        frames.append(g.assign(sample=s))
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    return out[["sample"] + [c for c in out.columns if c != "sample"]]


def group_wide(gl, samples, n=HEATMAP_TAXA):
    """member set x sample abundance of the multi-member groups, the n with the highest mean."""
    if gl.empty:
        return pd.DataFrame()
    w = gl.pivot_table(index="members", columns="sample", values="abundance", aggfunc="sum").reindex(
        columns=samples).fillna(0.0)
    mean = w.mean(axis=1)
    keep = sorted(w.index, key=lambda t: (-mean[t], str(t)))[:n]
    return w.loc[keep]


def cell_long(run):
    frames = [run["cell"][s].assign(sample=s) for s in run["samples"] if s in run["cell"]]
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    return out[["sample"] + [c for c in out.columns if c != "sample"]]


def read_assignment_summary(run, max_pairs=30):
    """Per sample: reads listed, reads by their highest probability on a species row, reads whose highest row is
    Below detection or Unassigned; and the species pairs that share the most reads (a read's two highest species
    rows, both at 0.01 or more as the file lists them, D37)."""
    per, pairs = [], {}
    for s in run["samples"]:
        p = run["reads"].get(s)
        if p is None:
            continue
        d = pd.read_csv(p, sep="\t", dtype={"read_id": str, "species": str, "probability": float},
                        keep_default_na=False)
        d = d.sort_values(["read_id", "probability", "species"], ascending=[True, False, True], kind="mergesort")
        top = d.drop_duplicates("read_id")
        row = {"sample": s, "reads": int(top.shape[0])}
        sp = top[~top["species"].isin(PLACEHOLDERS)]
        for label, lo, hi in TOP_PROB_CLASSES:
            row[label] = int(((sp["probability"] >= lo) & (sp["probability"] < hi)).sum())
        row["top: Below detection"] = int((top["species"] == BELOW).sum())
        row["top: Unassigned"] = int((top["species"] == UNASSIGNED).sum())
        two = d[~d["species"].isin(PLACEHOLDERS)].groupby("read_id", sort=False).head(2)
        multi = two.groupby("read_id", sort=False)["species"].agg(list)
        row["reads on 2+ species"] = int((multi.str.len() >= 2).sum())
        for lst in multi:
            if len(lst) >= 2:
                a, b = sorted(lst[:2])
                pairs[(s, a, b)] = pairs.get((s, a, b), 0) + 1
        per.append(row)
    pair_df = pd.DataFrame([{"sample": s, "species_a": a, "species_b": b, "reads": n}
                            for (s, a, b), n in pairs.items()], columns=["sample", "species_a", "species_b", "reads"])
    if not pair_df.empty:
        pair_df = pair_df.sort_values(["reads", "sample", "species_a", "species_b"],
                                      ascending=[False, True, True, True]).head(max_pairs).reset_index(drop=True)
    return pd.DataFrame(per), pair_df


def compute(run):
    """The data the figures and the model share."""
    data = {"samples": run["samples"], "group": run["group"], "fate": read_fate(run)}
    for rank in ("species", "genus"):
        data[f"composition_{rank}"], data[f"shown_{rank}"] = composition(run, rank)
    if len(run["samples"]) > HEATMAP_ABOVE:
        data["heatmap_species"] = heatmap_taxa(run)
    data["flagged_long"] = flagged_long(run)
    data["flagged"] = flagged_summary(data["flagged_long"])
    if run["groups"]:
        data["groups_long"] = groups_long(run)
        data["groups_wide"] = group_wide(data["groups_long"], run["samples"])
    if run["cell"]:
        data["cell_long"] = cell_long(run)
    if run["reads"]:
        data["reads"], data["read_pairs"] = read_assignment_summary(run)
    return data


# ---- formatting: tables as text for both documents (BugBuster's summary_report.py)

def fmt(v):
    """A cell as text: floats to 3 significant digits (thousands separated from 1,000), large ints separated,
    missing as blank."""
    if v is None:
        return ""
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else ""
    if isinstance(v, (float, np.floating)):
        if np.isnan(v):
            return ""
        if v == 0:
            return "0"
        if float(v).is_integer() and abs(v) >= 1000:
            return f"{v:,.0f}"
        return f"{v:.3g}" if abs(v) < 1000 else f"{v:,.0f}"
    if isinstance(v, (int, np.integer)):
        return f"{v:,}" if abs(v) >= 10000 else str(v)
    if v is pd.NA or v is pd.NaT:
        return ""
    return str(v)


def table_width_mm(header, rows, pt=6.5, wrap=36):
    width = 0.0
    for i, h in enumerate(header):
        longest = max([len(h)] + [len(r[i]) for r in rows])
        width += min(longest, wrap) * pt * 0.3528 * 0.56 + 2.5
    return width


class Numbering:
    def __init__(self):
        self.fig = 0
        self.tab = 0


def table_entry(nums, df, name, caption, html_rows=500, pdf=True, pdf_columns=None, tsv="", pdf_rows=PDF_MAX_ROWS):
    """A numbered table: header and text rows for the HTML (up to html_rows) and the PDF (capped). tsv: the
    published file that holds the full table, named in the notes."""
    d = df.reset_index(drop=True)
    nums.tab += 1
    total = len(d)
    hd = d.head(html_rows)
    out = {"no": nums.tab, "name": name, "caption": caption, "tsv": tsv,
           "header": [str(c) for c in hd.columns], "rows": [[fmt(v) for v in r] for r in hd.itertuples(index=False)],
           "n_rows": total, "html_note": "", "pdf": None}
    where = f"; the full table is {tsv}" if tsv else ""
    if total > html_rows:
        out["html_note"] = f"First {html_rows:,} of {total:,} rows{where}."
    if pdf:
        p = d[pdf_columns] if pdf_columns else d
        notes = []
        if p.shape[1] > PDF_MAX_COLS:
            notes.append(f"first {PDF_MAX_COLS} of {p.shape[1]} columns")
            p = p.iloc[:, :PDF_MAX_COLS]
        elif pdf_columns and len(pdf_columns) < d.shape[1]:
            notes.append(f"{len(pdf_columns)} of {d.shape[1]} columns")
        if total > pdf_rows:
            notes.append(f"first {pdf_rows} of {total:,} rows")
            p = p.head(pdf_rows)
        note = ("Shown: " + " and ".join(notes) + (where or "") + "; every row is in the HTML report.") if notes else ""
        rows = [[fmt(v) for v in r] for r in p.itertuples(index=False)]
        header = [str(c) for c in p.columns]
        out["pdf"] = {"header": header, "rows": rows, "note": note,
                      "landscape": table_width_mm(header, rows) > PORTRAIT_MM}
    return out


class Section:
    def __init__(self, sid, nums, records):
        self.id = sid
        self.title = dict((s, t) for s, t, _ in SECTIONS)[sid]
        self.page = PAGE[sid]
        self.nums, self.records = nums, records
        self.facts, self.figures, self.tables, self.lists = {}, {}, {}, {}

    def fig(self, key):
        r = next((x for x in self.records if x["key"] == key), None)
        if r is None:
            return None
        self.nums.fig += 1
        self.figures[key] = {"no": self.nums.fig, "key": key, "title": r["title"], "caption": r["caption"],
                             "parts": r["parts"], "interactive": r.get("interactive")}
        return self.figures[key]

    def table(self, slot, df, name, caption, **kw):
        if df is None:
            return None
        self.tables[slot] = table_entry(self.nums, df, name, caption, **kw)
        return self.tables[slot]

    def as_dict(self):
        return {"id": self.id, "title": self.title, "page": self.page, "facts": self.facts, "figures": self.figures,
                "tables": self.tables, "lists": self.lists}


# ---- references and software (BugBuster's summary_report.py)

def load_references(assets):
    """references.yml ({key: {author, title, journal, year, volume, pages, doi, type}}) and tool_references.tsv
    (software_versions.yml name -> display name, kind, citekeys)."""
    import yaml
    refs_path, map_path = Path(assets) / "references.yml", Path(assets) / "tool_references.tsv"
    refs = yaml.safe_load(refs_path.read_text(encoding="utf-8")) if refs_path.exists() else {}
    tool_refs = {}
    if map_path.exists():
        with open(map_path, encoding="utf-8") as fh:
            for row in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
                tool_refs[row["tool"]] = {"kind": row.get("kind") or "tool",
                                          "keys": [k for k in (row.get("references") or "").split(";") if k],
                                          "name": row.get("name") or row["tool"]}
    for t, r in tool_refs.items():
        for k in r["keys"]:
            if k not in (refs or {}):
                warn(f"tool_references.tsv: {t} cites {k}, which references.yml does not have")
    return refs or {}, tool_refs


def tools_rows(versions, tool_refs, extra=None):
    """One row per tool and version from software_versions.yml (process -> tool -> version; the Workflow block
    names NanoVI and Nextflow) plus the report's own (`extra`, the same form); then the 'method' entries of
    tool_references.tsv, which no process reports (Emu, NanoVI's preprint). Helpers are left out."""
    seen = {}
    blocks = dict(versions or {})
    for proc, tools in (extra or {}).items():
        blocks[proc] = tools
    for proc in sorted(blocks):
        tools = blocks[proc] or {}
        if not isinstance(tools, dict):
            continue
        for tool, ver in tools.items():
            info = tool_refs.get(tool, {"kind": "tool", "keys": [], "name": tool})
            if info["kind"] in ("helper", "method"):
                continue
            rec = seen.setdefault((info["name"], str(ver)), {"tool": info["name"], "version": str(ver),
                                                              "kind": info["kind"], "processes": set(),
                                                              "keys": info["keys"]})
            if proc != "Workflow":
                rec["processes"].add(proc)
    order = {"pipeline": 0, "tool": 1, "library": 2}
    rows = sorted(seen.values(), key=lambda r: (order.get(r["kind"], 3), r["tool"].lower(), r["version"]))
    for tool, info in sorted(tool_refs.items()):
        if info["kind"] == "method":
            rows.append({"tool": info["name"], "version": "", "kind": "method", "processes": set(),
                         "keys": info["keys"]})
    for r in rows:
        r["processes"] = ", ".join(sorted(r["processes"]))
    return rows


def format_reference(r):
    """'Authors. Title. Journal volume, pages (year). doi:...' as plain text (the HTML list)."""
    authors = r.get("author", "")
    parts = [a.strip() for a in authors.split(" and ")] if authors else []
    if len(parts) > 5 or (parts and parts[-1] == "others"):
        names = ", ".join(p for p in parts[:3] if p != "others") + " et al."
    else:
        names = ", ".join(parts[:-1]) + (" & " if len(parts) > 1 else "") + (parts[-1] if parts else "")
    venue = r.get("journal") or r.get("booktitle") or r.get("publisher") or ""
    vol = f" {r['volume']}" if r.get("volume") else ""
    pages = f", {str(r['pages']).replace('--', '-')}" if r.get("pages") else ""
    doi = f" doi:{r['doi']}" if r.get("doi") else ""
    return f"{names}. {r.get('title', '')}. {venue}{vol}{pages} ({r.get('year', '')}).{doi}"


def bibtex(refs, keys):
    """BibTeX for Typst's bibliography() of the keys cited, in that order."""
    out = []
    for k in keys:
        r = refs[k]
        fields = []
        for f in ("author", "title", "journal", "booktitle", "publisher", "year", "volume", "number", "pages", "doi", "url"):
            if r.get(f) not in (None, ""):
                v = str(r[f])
                fields.append(f"  {f} = {{{{{v}}}}}" if f == "title" else f"  {f} = {{{v}}}")
        out.append(f"@{r.get('type', 'article')}{{{k},\n" + ",\n".join(fields) + "\n}\n")
    return "\n".join(out)


PARAM_LABELS = [
    ("db", "Database (--db)"), ("taxonomy_tsv", "Taxonomy (--taxonomy_tsv)"), ("db_index", "Prebuilt index"),
    ("type", "minimap2 preset (--type)"), ("kmer_size", "minimap2 k-mer size"), ("N", "Secondary alignments (-N)"),
    ("K", "minimap2 batch (-K)"), ("min_length", "Minimum read length"), ("max_length", "Maximum read length"),
    ("min_quality", "Minimum mean read quality"), ("estimator", "Estimator"), ("min_abundance", "Detection floor"),
    ("vi_tolerance", "VB tolerance"), ("vi_max_iterations", "VB maximum iterations"), ("alpha_prior", "Dirichlet prior"),
    ("vi_elbo_threshold", "EM threshold"), ("binomial_names", "Binomial names"), ("keep_counts", "Estimated counts"),
    ("credible_intervals", "Credible intervals"), ("read_assignments", "Per-read assignments"),
    ("ambiguity_groups", "Ambiguity groups"), ("cell_abundance", "Cell-level table"),
    ("report_metadata", "Report metadata"),
]


def params_table(params):
    if not params:
        return None
    rows = [{"Parameter": label, "Option": f"--{k}", "Value": "" if params.get(k) is None else str(params[k])}
            for k, label in PARAM_LABELS if k in params]
    return pd.DataFrame(rows) if rows else None


# ---- prose shared by both documents

def caveats(data, run, params):
    out = [
        "Abundances are shares of 16S reads, not of cells: a species with more 16S copies per genome gives more "
        "reads per cell" + (" (the cell-level section corrects for this)." if run["cell"] else "."),
        "Below detection holds the reads of species pruned by the detection rules (Emu's rule: 1/R up to 1,000 "
        "reads, 10/R above, and the --min_abundance floor, D26); the thresholds, not the model, set the detection "
        "limit.",
        f"This report flags a species in a sample when its shared_read_fraction is at least {FLAG_SHARED:g} (most of "
        "its expected reads are split with another reported species) or when it shares an identical 16S with "
        "species of its genus (identical_16s_species, D24). The cutoff is a display rule, not a test.",
        "Species that share an identical 16S cannot be told apart from 16S reads: their reads go to whichever has "
        "the larger estimate (D24).",
    ]
    if any("ci_lower" in df.columns for df in run["species"].values()):
        out.append("ci_lower and ci_upper are sampling intervals of the variational posterior, not calibrated "
                   "uncertainty (D27).")
    if run["groups"]:
        out.append("Ambiguity groups join reported species of one family whose reads are shared (f >= 0.25, D39); a "
                   "group's abundance is more robust than the split between its members.")
    if run["cell"]:
        out.append("Cell-level abundance divides each species' read-level abundance by its 16S copy number "
                   "(copy_numbers.tsv of the database: GTDB ssu_count medians, filled from the genus or family "
                   "median when the species has none, D40) and renormalises; species without a copy number are "
                   "left out.")
    if run["reads"]:
        out.append("The per-read file lists each read's table rows at probability 0.01 or more (D37).")
    return out


# ---- the model

def build(run, data, records, assets, params=None, report_versions=None, date=""):
    """The report model (a JSON-ready dict)."""
    nums = Numbering()
    secs = []
    fate = data["fate"]
    samples = run["samples"]
    n = len(samples)
    has_group = any(run["group"].values())

    # overview
    ov = Section("overview", nums, records)
    sm = run.get("sample_metadata")
    tot = lambda c: int(fate[c].sum()) if c in fate and fate[c].notna().all() else None
    raw, passed = tot("raw_reads"), tot("passed_reads")
    sp_reads, bd, un = tot("species_reads"), tot("below_detection_reads"), tot("unassigned_reads")
    n_species = int((rank_abundance(run, "species").sum(axis=1) > 0).sum())
    ov.facts = {"n_samples": n, "raw_reads": fmt(raw), "passed_reads": fmt(passed),
                "pct_passed": f"{100 * passed / raw:.1f}" if raw and passed is not None else "",
                "species_reads": fmt(sp_reads),
                "pct_species": f"{100 * sp_reads / passed:.1f}" if passed and sp_reads is not None else "",
                "below_reads": fmt(bd), "unassigned_reads": fmt(un), "n_species": n_species,
                "median_species": fmt(float(fate["n_species"].median())) if n else "",
                "has_counts": sp_reads is not None, "has_group": has_group}
    cards = [(n, "samples")]
    if raw is not None:
        cards.append((raw, "raw reads"))
    if passed is not None:
        cards.append((passed, f"reads passed the length filter ({ov.facts['pct_passed']} %)" if raw else
                      "reads passed the length filter"))
    if sp_reads is not None:
        cards.append((sp_reads, f"reads on reported species ({ov.facts['pct_species']} %)" if passed else
                      "reads on reported species"))
    cards += [(n_species, "species reported (any sample)"),
              (len(data["flagged"]), "flagged species (any sample)")]
    ov.lists["cards"] = [{"value": fmt(v), "label": lab} for v, lab in cards]
    att = []
    none = [r["sample"] for _, r in fate.iterrows() if r["n_species"] == 0]
    if none:
        att.append(f"{len(none)} sample{'s' if len(none) > 1 else ''} with no reported species: {', '.join(none)}")
    if not ov.facts["has_counts"]:
        att.append("The tables have no estimated_counts (--keep_counts false): read fate shows fastplong's counts "
                   "only.")
    ov.lists["attention"] = att
    ov.fig("read_fate")
    cols = [c for c in ["sample", "group", "raw_reads", "passed_reads", "too_short", "too_long", "low_quality",
                        "chimeric", "species_reads", "below_detection_reads", "unassigned_reads", "pct_species",
                        "n_species"] if c in fate.columns]
    ov.table("fate", fate[cols], "read_fate", "Reads per sample: fastplong's raw and passed reads (too_short, "
             "too_long, low_quality: removed; chimeric: split into segments, so passed reads can exceed raw reads "
             "less the removed ones), and where the passed reads went: reported species, Below detection and "
             "Unassigned (estimated_counts of the species tables); pct_species is their share of the passed "
             "reads; n_species the species reported.", tsv="sample_outputs/<sample>_rel-abundance.tsv",
             pdf_columns=[c for c in ("sample", "raw_reads", "passed_reads", "species_reads",
                                      "below_detection_reads", "unassigned_reads", "pct_species", "n_species")
                          if c in fate.columns])
    if run["metadata"] is not None and run["metadata"].shape[1] > 1:
        ov.table("metadata", run["metadata"], "metadata", "Sample metadata (--report_metadata).",
                 tsv=str(params.get("report_metadata", "")) if params else "")
    secs.append(ov)

    # composition
    co = Section("composition", nums, records)
    co.facts = {"top": TOP_PER_SAMPLE, "cap": MAX_TAXA, "n_shown_species": len(data["shown_species"]),
                "n_shown_genus": len(data["shown_genus"]), "heatmap": "heatmap_species" in data,
                "heatmap_taxa": HEATMAP_TAXA, "heatmap_above": HEATMAP_ABOVE}
    for k in ("composition_species", "composition_genus", "heatmap_species"):
        co.fig(k)
    prev = rank_abundance(run, "species")
    prev_df = pd.DataFrame({"species": prev.index, "samples_detected": (prev > 0).sum(axis=1).values,
                            "mean_abundance": prev.mean(axis=1).values, "max_abundance": prev.max(axis=1).values})
    prev_df = prev_df[prev_df["samples_detected"] > 0].sort_values(
        ["samples_detected", "mean_abundance", "species"], ascending=[False, False, True]).head(100)
    co.table("prevalence", prev_df, "prevalence", "The 100 species detected in the most samples, with their mean "
             "and largest read-level abundance (mean over every sample, 0 where not detected).",
             tsv="consolidated_output/otu_table_abundance.tsv")
    secs.append(co)

    # flagged
    fl = Section("flagged", nums, records)
    fs_ = data["flagged"]
    fl_long = data["flagged_long"]
    fl.facts = {"cutoff": f"{FLAG_SHARED:g}", "n_flagged": len(fs_),
                "n_shared": int(fl_long.loc[fl_long["flag_shared"], "species"].nunique()) if not fl_long.empty else 0,
                "n_identical": int(fl_long.loc[fl_long["flag_identical"], "species"].nunique())
                if not fl_long.empty else 0}
    fl.fig("flagged_species")
    fl.table("flagged", fs_ if not fs_.empty else pd.DataFrame(columns=["species"]), "flagged",
             f"Species flagged in at least one sample: samples where reported, where shared_read_fraction >= "
             f"{FLAG_SHARED:g} and where identical_16s_species names a relative; median and largest abundance; the "
             "largest shared_read_fraction; the relatives.", tsv="sample_outputs/<sample>_rel-abundance.tsv",
             pdf_columns=["species", "samples_detected", "samples_shared", "samples_identical_16s",
                          "max_abundance", "max_shared_read_fraction", "identical_16s_species"]
             if not fs_.empty else None)
    secs.append(fl)

    if "groups_long" in data:
        gr = Section("groups", nums, records)
        gl = data["groups_long"]
        has_cell = "cell_abundance" in gl.columns
        gr.facts = {"n_samples": len(run["groups"]), "n_groups": len(gl),
                    "n_sets": int(gl["members"].nunique()) if not gl.empty else 0, "cell": has_cell}
        gr.fig("groups_heatmap")
        cols = [c for c in ("sample", "group", "members", "n_members", "abundance", "lowest_common_rank",
                            "lowest_common_taxon", "estimated_counts", "max_shared", "shared_read_fraction",
                            "identical_16s_relatives", "cell_abundance", "copy_number", "copy_number_level")
                if c in gl.columns] if not gl.empty else ["sample", "group", "members"]
        gr.table("groups", gl[cols] if not gl.empty else pd.DataFrame(columns=cols), "groups",
                 "Groups of two or more members in each sample (groups of one are the species rows)"
                 + ("; cell_abundance, copy_number and copy_number_level from --cell_abundance" if has_cell else "")
                 + ".", tsv="sample_outputs/<sample>_groups.tsv", html_rows=2000,
                 pdf_columns=[c for c in ("sample", "group", "members", "abundance", "lowest_common_taxon",
                                          "max_shared") if c in cols])
        secs.append(gr)

    if "cell_long" in data:
        ce = Section("cell", nums, records)
        cl = data["cell_long"]
        comments = run["cell_comments"]
        first = sorted({c[0] for c in comments.values() if c})
        ce.facts = {"n_samples": len(run["cell"]), "explanation": first[0] if len(first) == 1 else ""}
        ce.lists["comments"] = [{"sample": s, "lines": comments[s]} for s in samples if s in comments]
        ce.table("left_out", pd.DataFrame([{"sample": s, "not included": " ".join(comments[s][1:])
                                            if len(first) == 1 else " ".join(comments[s])}
                                           for s in samples if s in comments]), "cell_left_out",
                 "What each sample's cell-level table leaves out (its '#' lines).",
                 tsv="sample_outputs/<sample>_cell-abundance.tsv")
        ce.fig("cell_vs_read")
        cols = [c for c in ("sample", "species", "abundance", "read_abundance", "copy_number", "copy_number_level",
                            "genus") if c in cl.columns]
        ce.table("cell", cl[cols], "cell", "Cell-level abundance (abundance) beside the read-level abundance, the "
                 "16S copy number and the level it comes from.", tsv="sample_outputs/<sample>_cell-abundance.tsv",
                 html_rows=3000)
        secs.append(ce)

    if "reads" in data:
        rd = Section("reads", nums, records)
        rd.facts = {"n_samples": len(data["reads"])}
        rd.table("per_sample", data["reads"], "read_assignments", "Reads of each sample's per-read file by their "
                 "highest probability on a species row (p), reads whose highest row is Below detection or "
                 "Unassigned, and reads listed with two or more species.",
                 tsv="sample_outputs/<sample>_read_assignments.tsv.gz")
        rd.table("pairs", data["read_pairs"], "read_pairs", "The species pairs that share the most reads (a read's "
                 "two highest species rows).", tsv="")
        secs.append(rd)

    # methods
    me = Section("methods", nums, records)
    refs, tool_refs = load_references(assets)
    rows = tools_rows(run["software_versions"], tool_refs, report_versions)
    cited = []
    for r in rows:
        for k in r["keys"]:
            if k in refs and k not in cited:
                cited.append(k)
    me.table("params", params_table(params), "parameters", "Run parameters.", tsv="")
    tool_df = pd.DataFrame([{"Tool": r["tool"], "Version": r["version"], "Kind": r["kind"],
                             "Processes": r["processes"],
                             "Ref.": ",".join(str(cited.index(k) + 1) for k in r["keys"] if k in cited)}
                            for r in rows])
    if not tool_df.empty:
        e = me.table("tools", tool_df, "tools", "Software with its versions and the processes that used it (from "
                     "software_versions.yml and the report's own versions), and the methods NanoVI builds on.",
                     html_rows=100000, pdf_rows=100000)
        e["keys"] = [[k for k in r["keys"] if k in cited] for r in rows]
        e["pdf"]["landscape"] = False
    me.lists["caveats"] = caveats(data, run, params)
    me.lists["references"] = [{"no": i + 1, "key": k, "text": format_reference(refs[k])} for i, k in enumerate(cited)]
    me.facts["bibtex"] = bibtex(refs, cited)
    secs.append(me)

    wf = (run["software_versions"] or {}).get("Workflow", {}) or {}
    version = next((str(v) for k, v in wf.items() if "nanovi" in str(k).lower()), "")
    meta = {"title": "NanoVI run report", "generated": date, "n_samples": n, "nanovi_version": version,
            "nextflow_version": str(wf.get("Nextflow", "")), "pdf": "nanovi_report.pdf"}
    return {"run": meta, "sections": [s.as_dict() for s in secs]}
