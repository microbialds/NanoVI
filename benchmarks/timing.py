#!/usr/bin/env python3
"""Time per step for NanoVI and Emu at matched settings (Phase 7a, item 9; W3).

For each sample and minimap2 setting (-N, k) on NanoVI's exact database, the reads
are aligned once (minimap2 2.24, the pipeline's options, SAM out) and the same SAM
goes to both tools:

  align              minimap2 -ax map-ont -N <N> -k <k> -p 0.9 -K 500M (prebuilt index)
  sam_to_bam         samtools view -b (the pipeline pipes minimap2 into it)
  nanovi_likelihood  read_likelihood.py on the BAM (Emu score, primary error profile)
  nanovi_vb          estimate_abundance.py --estimator vb
  nanovi_em          estimate_abundance.py --estimator em (Emu's EM loop in NanoVI's code)
  emu                emu abundance <sam> --db exact (Emu 3.6.2, its likelihood and EM;
                     given a SAM it skips minimap2, checked identical on 200 reads)

plus Emu as it ships ('emu_shipped': the FASTQ, its minimap2 2.30 from the FASTA,
-N 50, k 15, its default database), and the index builds, once per k. NanoVI as it
ships is the -N 3, k 21 row (cached index, so align + sam_to_bam + likelihood + vb).

Each step runs under GNU time (user and system CPU, peak RSS of the largest process,
wall). Containers get --cpuset-cpus and host steps taskset, 4 CPUs per slot, 4 slots
at most (16 cores); /proc/loadavg is read before and after each step: the host is
shared, so CPU time is the primary number and wall time is reported with the load.
Repeats are interleaved (every unit's repeat 1 is queued before any repeat 2) and the
order of the steps after the alignment rotates with the repeat. Writes
<bench>/results/timing.tsv.
"""
import argparse
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager
from pathlib import Path

import pandas as pd

import harness
from common import BENCH, BIN, write_tsv
from comparators import EMU_DBS, EMU_IMAGE

SAMPLES = ["zymo_MSPlus_sup_rep1", "emu_zymo_SRR10391201"]
SETTINGS = [(3, 21), (50, 15), (50, 21), (3, 15)]
OUT = BENCH / "timing"
GTIME = "/usr/bin/time"
FMT = "%U %S %M %e"
SINGLE_THREAD = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                 "NUMBA_NUM_THREADS": "1"}


def reads_of(sample):
    return BENCH / "real" / "runs" / sample / "filtered_default.fastq.gz"


def loadavg():
    return float(Path("/proc/loadavg").read_text().split()[0])


def timed(step, argv, timefile, cpus, container=None, mounts=()):
    """Run argv under GNU time pinned to cpus; return the measurements."""
    timefile = Path(timefile)
    if container:
        cmd = ["docker", "run", "--rm", "-u", harness.UID, "--network", "none", "--label", "nanovi-phase7a",
               "--cpuset-cpus", cpus, "-v", f"{GTIME}:/opt/gtime:ro"]
        for m in sorted({str(timefile.parent), *map(str, mounts)}):
            cmd += ["-v", f"{m}:{m}"]
        cmd += [container, "/opt/gtime", "-f", FMT, "-o", str(timefile), "bash", "-euo", "pipefail", "-c", argv]
        env = None
    else:
        cmd = ["taskset", "-c", cpus, GTIME, "-f", FMT, "-o", str(timefile)] + argv
        env = {**os.environ, **SINGLE_THREAD}
    before = loadavg()
    with open(str(timefile) + ".log", "w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True, env=env)
    after = loadavg()
    user, system, rss_kb, wall = timefile.read_text().split()[-4:]
    return {"step": step, "user_s": float(user), "system_s": float(system),
            "cpu_s": float(user) + float(system), "wall_s": float(wall), "peak_rss_mb": int(rss_kb) / 1024,
            "load1_before": before, "load1_after": after, "cpus": cpus}


def matched_unit(sample, N, k, repeat, cpus):
    run = OUT / sample / f"N{N}_k{k}" / f"r{repeat}"
    run.mkdir(parents=True, exist_ok=True)
    reads, index, db = reads_of(sample), harness.index_path("exact", k=k), EMU_DBS["nanovi_exact"]
    sam, bam = run / "aln.sam", run / "aln.bam"
    rows = [timed("align", f"minimap2 -ax map-ont -t 4 -N {N} -p 0.9 -K 500M -k {k} {index} {reads} > {sam}",
                  run / "align.time", cpus, harness.ALIGN_IMAGE["2.24"], [reads.parent, index.parent])]
    rows.append(timed("sam_to_bam", f"samtools view -b -@ 2 -o {bam} {sam}", run / "sam_to_bam.time", cpus,
                      harness.ALIGN_IMAGE["2.24"]))
    npz = run / "scores.npz"
    nanovi = [
        ("nanovi_likelihood", [harness.PYTHON, str(BIN / "read_likelihood.py"), "--alignments", str(bam),
                               "--scores", str(npz), "--counts", str(run / "counts.json")]),
        ("nanovi_vb", [harness.PYTHON, str(BIN / "estimate_abundance.py"), "--scores", str(npz),
                       "--output", str(run / "vb.json"), "--estimator", "vb"]),
        ("nanovi_em", [harness.PYTHON, str(BIN / "estimate_abundance.py"), "--scores", str(npz),
                       "--output", str(run / "em.json"), "--estimator", "em"]),
    ]
    emu = ("emu", f"emu abundance {sam} --db {db} --keep-counts --output-dir {run / 'emu'}")
    # rotate which tool goes first after the alignment
    order = ["nanovi", "emu"] if repeat % 2 else ["emu", "nanovi"]
    for tool in order:
        if tool == "emu":
            rows.append(timed(emu[0], emu[1], run / "emu.time", cpus, EMU_IMAGE, [db]))
        else:
            rows += [timed(step, argv, run / f"{step}.time", cpus) for step, argv in nanovi]
    reads_n = sum(1 for _ in pysam_names(bam))
    sam.unlink()
    bam.unlink()
    return [{"sample": sample, "N": N, "k": k, "repeat": repeat, "mode": "matched", "reads": reads_n, **r}
            for r in rows]


def pysam_names(bam):
    import pysam
    with pysam.AlignmentFile(str(bam)) as f:
        for a in f:
            if not (a.is_secondary or a.is_supplementary):
                yield a.query_name


def shipped_unit(sample, repeat, cpus):
    run = OUT / sample / "emu_shipped" / f"r{repeat}"
    if run.exists():
        shutil.rmtree(run)
    run.mkdir(parents=True)
    reads, db = reads_of(sample), EMU_DBS["emu_default"]
    row = timed("emu_shipped", f"emu abundance {reads} --db {db} --threads 4 --keep-counts --output-dir {run}",
                run / "emu_shipped.time", cpus, EMU_IMAGE, [reads.parent, db])
    return [{"sample": sample, "N": 50, "k": 15, "repeat": repeat, "mode": "emu_shipped", **row}]


def index_unit(k, cpus):
    run = OUT / "index"
    run.mkdir(parents=True, exist_ok=True)
    db = EMU_DBS["nanovi_exact"]
    out = run / f"exact.k{k}.mmi"
    row = timed("index", f"minimap2 -x map-ont -k {k} -t 4 -d {out} {db}/species_taxid.fasta",
                run / f"index_k{k}.time", cpus, harness.ALIGN_IMAGE["2.24"], [db])
    out.unlink()
    return [{"sample": "", "N": "", "k": k, "repeat": 1, "mode": "index", **row}]


def run_job(job, slots):
    cpus = slots.get()
    try:
        kind, args = job
        return {"matched": matched_unit, "shipped": shipped_unit, "index": index_unit}[kind](*args, cpus)
    finally:
        slots.put(cpus)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--cpus", default="0-15", help="CPU range split into slots of 4 (default 0-15)")
    ap.add_argument("--samples", default=",".join(SAMPLES))
    args = ap.parse_args()
    lo, hi = map(int, args.cpus.split("-"))
    slot_cpus = [f"{c}-{c + 3}" for c in range(lo, hi + 1, 4)]
    samples = args.samples.split(",")
    for k in sorted({k for _, k in SETTINGS}):
        harness.build_index(BENCH / "db" / "exact", harness.index_path("exact", k=k), threads=16, k=k)
    jobs = [("index", (k,)) for k in sorted({k for _, k in SETTINGS})]
    for repeat in range(1, args.repeats + 1):
        jobs += [("matched", (s, N, k, repeat)) for s in samples for N, k in SETTINGS]
        jobs += [("shipped", (s, repeat)) for s in samples]
    with Manager() as m:
        slots = m.Queue()
        for c in slot_cpus:
            slots.put(c)
        with ProcessPoolExecutor(len(slot_cpus)) as pool:
            rows = [r for rs in pool.map(run_job, jobs, [slots] * len(jobs)) for r in rs]
    write_tsv(pd.DataFrame(rows), BENCH / "results" / "timing.tsv")


if __name__ == "__main__":
    main()
