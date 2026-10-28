#!/usr/bin/env python3
"""Clinical cohorts for the manuscript (Phase 9, decision D34; docs/benchmarks/clinical_2026.md).

Three ONT 16S cohorts, four arms each:

  nanovi_gtdb   nextflow run main.nf -profile docker, every default, GTDB r226 (default build)
  nanovi_emudb  the same pipeline on Emu's default database
  emu           Emu 3.6.2 as shipped (-N 50, k 15) on the reads NanoVI kept after filtering
  emu_same_sam  Emu 3.6.2 on nanovi_emudb's alignments (BAM -> SAM), the same alignments

Data and runs live under <bench>/../phase9 (D11); the pipeline runs from a copy of the
repository at the commit under test (src/, written by `prepare`). Steps:

    python benchmarks/clinical.py prepare  <cohort>   # subsample, samplesheet, src/
    python benchmarks/clinical.py nanovi   <cohort> <gtdb|emudb>
    python benchmarks/clinical.py emu      <cohort> [--same-sam]
    python benchmarks/clinical.py analyse             # docs/benchmarks/data/clinical_*.tsv
"""
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

import harness
from common import BENCH, DOCS_DATA, REPO, write_tsv  # first: puts bin/ on sys.path
from build_gtdb_db import normalize_gtdb_name
from comparators import EMU_DBS, EMU_IMAGE, run_emu
from names import load as load_crosswalk, to_gtdb

CLIN = BENCH.parent / "phase9"
SEED = 20261022
N_READS = 10000
SEQKIT_IMAGE = "quay.io/biocontainers/seqkit:2.14.0--hb192632_0"
SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"
PHASE7B = BENCH.parent / "phase7b"
GTDB_DB = PHASE7B / "db" / "suffixes2"  # default build with the current identical_16s.tsv (Phase 8)
GTDB_INDEX = PHASE7B / "confirm" / "suffixes" / "index_cache" / "species_taxid.a35597d73794.map-ont.k21.minimap2-2.31.mmi"
PLACEHOLDERS = {"Below detection", "Unassigned"}
EMU_SKIP = {"unassigned", "unmapped", "mapped_unclassified"}

COHORTS = {
    # key: (BioProject, subsample to N_READS?)
    "vaginal_cl": ("PRJNA1232418", False),
    "vaginal_bv": ("PRJNA723982", True),
    "crc": ("PRJNA911189", True),
}
ARMS = ["nanovi_gtdb", "nanovi_emudb", "emu", "emu_same_sam"]
PILOT = False  # --pilot: the first sample of each cohort, under <phase9>/pilot
# The four-species CRC model of the cohort paper (Aja-Macaya et al. 2025, Table 1), NCBI names
CRC_SPECIES = ["Parvimonas micra", "Fusobacterium nucleatum", "Bacteroides fragilis",
               "Agathobaculum butyriciproducens"]


def manifest(cohort):
    project = COHORTS[cohort][0]
    m = pd.read_csv(CLIN / "manifests" / f"{project}.ena.tsv", sep="\t")
    m = m[(m["instrument_platform"] == "OXFORD_NANOPORE") & (m["library_strategy"] == "AMPLICON")]
    m = m.sort_values("run_accession").reset_index(drop=True)
    return m.head(1) if PILOT else m


def labels(cohort):
    """{run: label} from the sample aliases: M1-M20; P01-P12 'cg' control / 'vg' BV; CCR / SANO."""
    m = manifest(cohort)
    out = {}
    for run, alias in zip(m["run_accession"], m["sample_alias"]):
        if cohort == "vaginal_cl":
            out[run] = "M%02d" % int(alias.rsplit("_M", 1)[1])
        elif cohort == "vaginal_bv":
            out[run] = ("BV" if "_vg_" in alias else "control", alias.split("_")[0])
        else:
            out[run] = "CRC" if alias.upper().startswith("CCR") else "healthy"
    return out


def run_dir(cohort, arm):
    return CLIN / ("pilot" if PILOT else "runs") / cohort / arm


def prepare(cohort):
    src = CLIN / "src"
    head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True,
                          check=True).stdout.strip()
    if not (src / "COMMIT").exists():
        # extract beside src/ and move it in place, so a failed archive leaves no half-filled src/
        tmp = CLIN / "src.partial"
        subprocess.run(["rm", "-rf", str(tmp)], check=True)
        tmp.mkdir(parents=True)
        subprocess.run(f"git -C {REPO} archive HEAD | tar -x -C {tmp}", shell=True, check=True)
        (tmp / "COMMIT").write_text(head + "\n")
        subprocess.run(["rm", "-rf", str(src)], check=True)
        tmp.rename(src)
    elif (src / "COMMIT").read_text().strip() != head:
        # every cohort runs the same pipeline code; remove src/ to run another commit
        print(f"note: runs use the pipeline at {(src / 'COMMIT').read_text().strip()} (src/), not HEAD {head}")
    reads = CLIN / "reads" / cohort
    reads.mkdir(parents=True, exist_ok=True)
    rows = []
    project, subsample = COHORTS[cohort]
    for run in manifest(cohort)["run_accession"]:
        raw = CLIN / "downloads" / project / f"{run}.fastq.gz"
        out = reads / f"{run}.fastq.gz"
        if subsample:
            # two steps: under pipefail, 'seqkit head' closing the pipe would fail the shuffle (SIGPIPE)
            shuf = reads / f"shuffled_{run}.fastq.gz"
            cmd = (f"seqkit shuffle -s {SEED} -j 2 {raw} -o {shuf} && seqkit head -n {N_READS} {shuf} -o {{OUT}} "
                   f"&& rm {shuf}")
            harness._run(out, harness.docker(SEQKIT_IMAGE, cmd, [raw.parent, reads]))
        elif not out.exists():
            out.symlink_to(raw)
        rows.append({"sample": run, "fastq": str(out)})
    d = run_dir(cohort, "")
    d.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(d / "samplesheet.csv", index=False)
    (d / "limit.config").write_text(
        "// Phase 9 clinical runs: at most 16 CPUs of the shared host (D20)\n"
        "executor {\n    cpus   = 16\n    memory = '96 GB'\n}\n")


def nanovi(cohort, db):
    arm = f"nanovi_{db}"
    d = run_dir(cohort, arm)
    d.mkdir(parents=True, exist_ok=True)
    base = run_dir(cohort, "")
    dbargs = (["--db", str(GTDB_DB), "--db_index", str(GTDB_INDEX)] if db == "gtdb"
              else ["--db", str(EMU_DBS["emu_default"])])
    cmd = ["nextflow", "-log", str(d / "nextflow.log"), "run", str(CLIN / "src" / "main.nf"), "-profile", "docker",
           "-c", str(base / "limit.config"), "--input", str(base / "samplesheet.csv"), *dbargs,
           "--output_dir", str(d / "out"), "--index_cache_dir", str(CLIN / "index_cache"), "-w", str(d / "work"),
           "--keep_files", "true", "--max_cpus", "8", "--max_memory", "64.GB",
           "-with-trace", str(d / "trace.tsv")]
    with open(d / "run.log", "w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True, cwd=d)


def sam_from_bam(bam, out):
    cmd = f"samtools view -h -o {{OUT}} {bam}"
    return harness._run(out, harness.docker(SAMTOOLS_IMAGE, cmd, [Path(bam).parent, Path(out).parent]))


def emu(cohort, same_sam, slots=4, threads=4):
    samples = list(manifest(cohort)["run_accession"])
    arm = "emu_same_sam" if same_sam else "emu"
    out = run_dir(cohort, arm)

    def one(run):
        if same_sam:
            # Emu takes a SAM in place of the reads and skips minimap2
            bam = run_dir(cohort, "nanovi_emudb") / "out" / "bam" / f"{run}_vi_alignments.bam"
            sam = sam_from_bam(bam, out / "sam" / f"{run}.sam")
            db, tables = EMU_DBS["emu_default"], out / "tables"
            cmd = (f"emu abundance {sam} --db {db} --threads {threads} --keep-counts "
                   f"--output-dir {tables}/tmp_{run} --output-basename {run} "
                   f"&& mv {tables}/tmp_{run}/{run}_rel-abundance.tsv {{OUT}}")
            return harness._run(tables / f"{run}_rel-abundance.tsv",
                                harness.docker(EMU_IMAGE, cmd, [sam.parent, db, tables]))
        else:
            reads = run_dir(cohort, "nanovi_gtdb") / "out" / "filtered" / f"{run}_filtered.fastq.gz"
            link = out / "reads" / f"{run}.fastq.gz"
            link.parent.mkdir(parents=True, exist_ok=True)
            if not link.exists():
                link.hardlink_to(reads)  # a symlink would point outside the container's mounts
            reads = link
        return run_emu(reads, "emu_default", out / "tables", threads=threads)

    with ThreadPoolExecutor(slots) as pool:
        list(pool.map(one, samples))


# --- analysis -----------------------------------------------------------------------------

def nanovi_table(cohort, arm, run):
    df = pd.read_csv(run_dir(cohort, arm) / "out" / "sample_outputs" / f"{run}_rel-abundance.tsv", sep="\t",
                     keep_default_na=False)
    df = df[~df["species"].isin(PLACEHOLDERS)]
    return {sp: float(a) for sp, a in zip(df["species"], df["abundance"]) if float(a) > 0}


def emu_table(path):
    t = pd.read_csv(path, sep="\t", dtype={"tax_id": str})
    if "tax_id" in t:
        t = t[~t["tax_id"].isin(EMU_SKIP)]
    t = t[t["species"].notna()]
    est = {}
    for sp, ab in zip(t["species"], t["abundance"]):
        est[sp] = est.get(sp, 0.0) + float(ab)
    return est


def normalise(est):
    total = sum(est.values())
    return {k: v / total for k, v in est.items()} if total > 0 else {}


def to_gtdb_profile(est, xwalk, gtdb_names):
    """NCBI-named profile -> GTDB species (suffixes stripped); returns (profile, unmapped share)."""
    est = normalise(est)
    out, unmapped = {}, 0.0
    for name, ab in est.items():
        g = normalize_gtdb_name(to_gtdb(name, xwalk))  # unknown names come back unchanged
        if g in gtdb_names:
            out[g] = out.get(g, 0.0) + ab
        else:
            unmapped += ab
    return normalise(out), unmapped


def strip(est):
    out = {}
    for k, v in est.items():
        n = normalize_gtdb_name(k)
        out[n] = out.get(n, 0.0) + v
    return normalise(out)


def genus_of(name, genus_map=None):
    """Genus of a species name: from genus_map (Emu's taxonomy.tsv, for NCBI names) when it
    gives one, else the first word ('Candidatus X' kept whole; GTDB names, published tables)."""
    if genus_map and genus_map.get(name):
        return genus_map[name]
    words = name.split()
    if not words:
        return name
    return " ".join(words[:2]) if words[0] == "Candidatus" and len(words) > 1 else words[0]


def emu_genera():
    t = pd.read_csv(EMU_DBS["emu_default"] / "taxonomy.tsv", sep="\t", dtype=str)
    t = t[t["species"].notna() & t["genus"].notna()]
    return dict(zip(t["species"], t["genus"]))


def genus(est, genus_map=None):
    out = {}
    for k, v in est.items():
        g = genus_of(k, genus_map)
        out[g] = out.get(g, 0.0) + v
    return out


def bray_curtis(a, b):
    keys = set(a) | set(b)
    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in keys)


OTHERS = "Others"


def fold(est, ref):
    """When a published profile lists only its main taxa plus an 'Others' bin, put
    every taxon of est that ref does not name into 'Others' as well."""
    if OTHERS not in ref:
        return est
    out = {}
    for k, v in est.items():
        key = k if k in ref else OTHERS
        out[key] = out.get(key, 0.0) + v
    return out


def published_bins(est):
    """'Not defined' (taxa the paper's genus table could not place) joins 'Others'."""
    out = {}
    for k, v in est.items():
        key = OTHERS if k in ("Not defined", "Other") else k
        out[key] = out.get(key, 0.0) + v
    return out


def dominant(est):
    return max(est, key=est.get) if est else ""


def profiles(cohort, run, xwalk, gtdb_names):
    """Every arm's profile of one sample in two naming systems: 'ncbi' (Emu's names, for
    the Emu-database arms) and 'gtdb' (every arm mapped to GTDB species, stripped)."""
    p, unmapped = {}, {}
    p["nanovi_gtdb"] = {"gtdb": strip(nanovi_table(cohort, "nanovi_gtdb", run))}
    for arm in ("nanovi_emudb", "emu", "emu_same_sam"):
        est = (nanovi_table(cohort, arm, run) if arm == "nanovi_emudb"
               else emu_table(run_dir(cohort, arm) / "tables" / f"{run}_rel-abundance.tsv"))
        g, unmapped[arm] = to_gtdb_profile(est, xwalk, gtdb_names)
        p[arm] = {"ncbi": normalise(est), "gtdb": g}
    return p, unmapped


def published(cohort, xwalk, gtdb_names):
    """{sample label: {'ncbi': profile, 'gtdb': profile, 'genus_ncbi': profile}} from the papers."""
    pub = CLIN / "downloads" / "published"
    out = {}
    if cohort == "vaginal_cl":
        # Supplementary Tables 5 (species) and 4 (genus) of the paper's workbook
        # (41598_2025_23837_MOESM2_ESM.xlsx), exported to TSV unchanged
        for name, key in (("oliva2025_table_S5_species.tsv", "ncbi"), ("oliva2025_table_S4_genus.tsv", "genus_ncbi")):
            t = pd.read_csv(pub / name, sep="\t")
            for sample, taxon, ab in t.itertuples(index=False):
                out.setdefault(sample, {}).setdefault(key, {})[taxon] = float(ab)
        for s in out:
            # map the named taxa; the 'Others' bin (and what fails to map) stays as 'Others'
            ncbi = published_bins(out[s]["ncbi"])
            rest = ncbi.pop(OTHERS, 0.0)
            named = sum(ncbi.values())
            g, unm = to_gtdb_profile(ncbi, xwalk, gtdb_names)
            out[s]["gtdb"] = {k: v * named * (1 - unm) for k, v in g.items()}
            out[s]["gtdb"][OTHERS] = rest + named * unm
            out[s]["unmapped"] = named * unm
    elif cohort == "vaginal_bv":
        res = pub / "emu_benchmark" / "results_clinicial_application"
        for i in range(1, 13):
            s = "P%02d" % i
            emu = emu_table(res / "Emu" / f"barcode{i:02d}_rel-abundance.tsv")
            br = pd.read_csv(res / "Bracken" / "genus" / f"{i}_bracken.rst", sep="\t")
            br = br[br["name"] != "Homo"]
            g, unm = to_gtdb_profile(emu, xwalk, gtdb_names)
            out[s] = {"ncbi": normalise(emu), "gtdb": g, "unmapped": unm,
                      "wgs_genus_ncbi": normalise(dict(zip(br["name"], br["new_est_reads"].astype(float)))),
                      "wgs_bacterial_reads": int(br["new_est_reads"].sum())}
    return out


def mann_whitney(x, y):
    from scipy.stats import mannwhitneyu
    if len(x) == 0 or len(y) == 0:
        return float("nan"), float("nan")
    u, p = mannwhitneyu(x, y, alternative="two-sided")
    return p, u / (len(x) * len(y))


def bh(p):
    s = pd.Series(p).dropna().sort_values()
    n = len(s)
    q = (s * n / pd.Series(range(1, n + 1), index=s.index)).iloc[::-1].cummin().iloc[::-1].clip(upper=1)
    return q.reindex(pd.Series(p).index)


def analyse():
    xwalk = load_crosswalk()
    gtdb_names = {normalize_gtdb_name(n) for n in pd.read_csv(GTDB_DB / "taxonomy.tsv", sep="\t")["species"]}
    pairs = [("nanovi_emudb", "emu_same_sam", "ncbi"), ("nanovi_emudb", "emu", "ncbi"), ("nanovi_gtdb", "emu", "gtdb"),
             ("nanovi_gtdb", "nanovi_emudb", "gtdb")]
    rows, samples = [], []
    crc_rows, paper_rows, contrib_rows = [], [], []
    gmap = emu_genera()
    db_species = database_species(gtdb_names, gmap)
    for cohort in COHORTS:
        if not (run_dir(cohort, "emu_same_sam") / "tables").exists():
            continue
        lab = labels(cohort)
        pub = published(cohort, xwalk, gtdb_names)
        crc_profiles = {arm: {} for arm in ARMS}
        contrib = {"species": {}, "genus": {}}
        for run in manifest(cohort)["run_accession"]:
            p, unmapped = profiles(cohort, run, xwalk, gtdb_names)
            label = lab[run]
            sample_label = label if cohort == "vaginal_cl" else (label[1] if cohort == "vaginal_bv" else run)
            group = label[0] if cohort == "vaginal_bv" else (label if cohort == "crc" else "")
            if sample_label in pub:
                ref = pub[sample_label]
                genus_ref = published_bins(ref.get("genus_ncbi") or genus(ref["ncbi"], gmap))
                p["published"] = {"ncbi": published_bins(ref["ncbi"]), "gtdb": published_bins(ref["gtdb"]),
                                  "genus": normalise(genus_ref)}
                unmapped["published"] = ref["unmapped"]
            base = {"cohort": cohort, "run": run, "sample": sample_label, "group": group}
            for a, b, names in pairs + ([("nanovi_gtdb", "published", "gtdb"), ("nanovi_emudb", "published", "ncbi"),
                                         ("emu", "published", "ncbi")] if "published" in p else []):
                pa, pb = p[a][names], p[b][names]
                gm = gmap if names == "ncbi" else None  # GTDB names: the genus is the first word
                ga, gb = normalise(genus(pa, gm)), normalise(genus(pb, gm))
                if b == "published":
                    if names == "ncbi":
                        gb = p[b]["genus"]  # the paper's own genus table
                    pa, ga = fold(pa, pb), fold(ga, gb)
                rows.append({**base, "arm_a": a, "arm_b": b, "names": names,
                             "species_bray_curtis": bray_curtis(pa, pb),
                             "genus_bray_curtis": bray_curtis(ga, gb),
                             "dominant_genus_a": dominant({k: v for k, v in ga.items() if k != OTHERS}),
                             "dominant_genus_b": dominant({k: v for k, v in gb.items() if k != OTHERS}),
                             "unmapped_a": unmapped.get(a, 0.0) if names == "gtdb" else 0.0,
                             "unmapped_b": unmapped.get(b, 0.0) if names == "gtdb" else 0.0})
            srow = {**base}
            for arm in p:
                prof = p[arm]["gtdb"]
                named = {k: v for k, v in genus(prof).items() if k != OTHERS}
                srow[f"{arm}_species"] = sum(1 for k in prof if k != OTHERS)
                srow[f"{arm}_lactobacillus"] = named.get("Lactobacillus", 0.0)
                srow[f"{arm}_dominant_genus"] = dominant(named)
            if cohort == "vaginal_bv":
                w = pub[sample_label]
                srow["wgs_bacterial_reads"] = w["wgs_bacterial_reads"]
                srow["wgs_lactobacillus"] = w["wgs_genus_ncbi"].get("Lactobacillus", 0.0)
                srow["wgs_dominant_genus"] = dominant(w["wgs_genus_ncbi"])
                for arm in ("nanovi_emudb", "emu"):
                    srow[f"{arm}_genus_bc_vs_wgs"] = bray_curtis(normalise(genus(p[arm]["ncbi"], gmap)),
                                                                 w["wgs_genus_ncbi"])
            samples.append(srow)
            # what separates the two databases: each taxon's share of the Bray-Curtis
            # dissimilarity between NanoVI on GTDB and on Emu's database (GTDB names)
            a, b = p["nanovi_gtdb"]["gtdb"], p["nanovi_emudb"]["gtdb"]
            for level, (x, y) in (("species", (a, b)), ("genus", (normalise(genus(a)), normalise(genus(b))))):
                for k in set(x) | set(y):
                    contrib[level][k] = contrib[level].get(k, 0.0) + 0.5 * abs(x.get(k, 0.0) - y.get(k, 0.0))
            if cohort == "crc":
                for arm in ARMS:
                    crc_profiles[arm][run] = (p[arm]["gtdb"], group)
        if cohort == "crc":
            for arm in ARMS:
                prof = crc_profiles[arm]
                species = sorted({k for v, _ in prof.values() for k in v})
                n = len(prof)
                res = []
                for sp in species:
                    present = sum(1 for v, _ in prof.values() if sp in v)
                    if present < 0.1 * n:
                        continue
                    crc = [v.get(sp, 0.0) for v, g in prof.values() if g == "CRC"]
                    healthy = [v.get(sp, 0.0) for v, g in prof.values() if g == "healthy"]
                    if not crc or not healthy:  # the pilot has one sample
                        continue
                    pval, auc = mann_whitney(crc, healthy)
                    res.append({"arm": arm, "species": sp, "prevalence": present / n, "p": pval, "auc_crc": auc,
                                "mean_crc": sum(crc) / len(crc), "mean_healthy": sum(healthy) / len(healthy)})
                if res:
                    df = pd.DataFrame(res)
                    df["q"] = bh(df["p"]).values
                    crc_rows.append(df)
                    paper_rows += paper_species(arm, prof, df, xwalk)
        n = len(manifest(cohort))
        for level, c in contrib.items():
            total = sum(c.values())
            for k in sorted(c, key=c.get, reverse=True)[:10]:
                g = genus_of(k) if level == "species" else k
                contrib_rows.append({"cohort": cohort, "level": level, "taxon": k, "mean_bc_contribution": c[k] / n,
                                     "share_of_bc": c[k] / total if total else float("nan"),
                                     "genus_species_gtdb": db_species["gtdb"].get(g, 0),
                                     "genus_species_emudb": db_species["emudb"].get(g, 0)})
    dest = CLIN / "pilot" if PILOT else DOCS_DATA
    pair_df = pd.DataFrame(rows)
    write_tsv(pair_df, dest / "clinical_pairs.tsv")
    write_tsv(pd.DataFrame(samples), dest / "clinical_samples.tsv")
    summary = (pair_df.assign(same=pair_df["dominant_genus_a"] == pair_df["dominant_genus_b"])
               .groupby(["cohort", "arm_a", "arm_b", "names"])
               .agg(samples=("run", "count"), species_bc_median=("species_bray_curtis", "median"),
                    species_bc_max=("species_bray_curtis", "max"), genus_bc_median=("genus_bray_curtis", "median"),
                    genus_bc_max=("genus_bray_curtis", "max"), dominant_genus_agree=("same", "mean"))
               .reset_index())
    write_tsv(summary, dest / "clinical_summary.tsv")
    if crc_rows:
        crc = pd.concat(crc_rows, ignore_index=True)
        write_tsv(crc, dest / "clinical_crc_species.tsv")
        write_tsv(crc_tests(crc), dest / "clinical_crc_tests.tsv")
        write_tsv(pd.DataFrame(paper_rows), dest / "clinical_crc_paper_species.tsv")
    sample_df = pd.DataFrame(samples)
    if (sample_df["cohort"] == "vaginal_bv").any():
        write_tsv(bv_lactobacillus(sample_df[sample_df["cohort"] == "vaginal_bv"]), dest / "clinical_bv_lactobacillus.tsv")
    write_tsv(pd.DataFrame(contrib_rows), dest / "clinical_db_contribution.tsv")


def database_species(gtdb_names, genus_map):
    """Species per genus in each database: GTDB r226 (default build, suffixes stripped, the
    names the comparisons use) and Emu's default database (genus from its taxonomy.tsv). The
    genus names are each database's own: a genus GTDB renamed or split is counted under
    different names in the two."""
    emu = set(pd.read_csv(EMU_DBS["emu_default"] / "taxonomy.tsv", sep="\t")["species"].dropna())
    count = lambda names, gm: pd.Series([genus_of(s, gm) for s in names]).value_counts().to_dict()
    return {"gtdb": count(gtdb_names, None), "emudb": count(emu, genus_map)}


def paper_species(arm, prof, df, xwalk):
    """The species of the CRC paper's four-species model (CRC_SPECIES) in one arm: prevalence
    and means for every one, and the test where the species passed the prevalence filter."""
    out = []
    for name in CRC_SPECIES:
        g = normalize_gtdb_name(to_gtdb(name, xwalk))
        crc = [v.get(g, 0.0) for v, grp in prof.values() if grp == "CRC"]
        healthy = [v.get(g, 0.0) for v, grp in prof.values() if grp == "healthy"]
        row = df[df["species"] == g]
        out.append({"arm": arm, "paper_species": name, "gtdb_species": g,
                    "prevalence": sum(1 for v, _ in prof.values() if g in v) / len(prof),
                    "mean_crc": sum(crc) / len(crc), "mean_healthy": sum(healthy) / len(healthy),
                    "tested": not row.empty,
                    "p": row["p"].iloc[0] if not row.empty else float("nan"),
                    "q": row["q"].iloc[0] if not row.empty else float("nan"),
                    "auc_crc": row["auc_crc"].iloc[0] if not row.empty else float("nan")})
    return out


def crc_tests(crc):
    """Per arm: species tested, species below q 0.05 and 0.10, and those species."""
    rows = []
    for arm, d in crc.groupby("arm", sort=False):
        rows.append({"arm": arm, "tested": len(d), "q_below_0.05": int((d["q"] < 0.05).sum()),
                     "q_below_0.10": int((d["q"] < 0.10).sum()),
                     "species_q_below_0.10": "; ".join(d[d["q"] < 0.10].sort_values("q")["species"])})
    return pd.DataFrame(rows)


def bv_lactobacillus(s):
    """Lactobacillus share in control and BV swabs per profile (16S arms, the published Emu
    tables and the shotgun Bracken genera), and the AUC of the share for control vs BV."""
    from scipy.stats import mannwhitneyu
    rows = []
    for arm in ARMS + ["published", "wgs"]:
        ctrl = s.loc[s["group"] == "control", f"{arm}_lactobacillus"].dropna()
        bv = s.loc[s["group"] == "BV", f"{arm}_lactobacillus"].dropna()
        if ctrl.empty or bv.empty:  # the pilot has one swab
            continue
        u, p = mannwhitneyu(ctrl, bv, alternative="two-sided")
        rows.append({"arm": arm, "control_n": len(ctrl), "control_min": ctrl.min(), "control_max": ctrl.max(),
                     "bv_n": len(bv), "bv_min": bv.min(), "bv_max": bv.max(),
                     "auc_control_higher": u / (len(ctrl) * len(bv)), "p": p})
    return pd.DataFrame(rows)


def main():
    global PILOT
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pilot", action="store_true", help="first sample of each cohort, under <phase9>/pilot")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare"); p.add_argument("cohort", choices=COHORTS)
    p = sub.add_parser("nanovi"); p.add_argument("cohort", choices=COHORTS); p.add_argument("db", choices=["gtdb", "emudb"])
    p = sub.add_parser("emu"); p.add_argument("cohort", choices=COHORTS); p.add_argument("--same-sam", action="store_true")
    sub.add_parser("analyse")
    args = ap.parse_args()
    PILOT = args.pilot
    if args.cmd == "prepare":
        prepare(args.cohort)
    elif args.cmd == "nanovi":
        nanovi(args.cohort, args.db)
    elif args.cmd == "emu":
        emu(args.cohort, args.same_sam)
    else:
        analyse()


if __name__ == "__main__":
    main()
