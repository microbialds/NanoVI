"""Steps of the abundance pipeline, run outside Nextflow so that each layer is cached
and reused across benchmark arms (Phase 7a).

The commands are those of the pipeline's modules, with the pipeline's containers:
  FILTER_READS  fastplong 0.2.2 (modules/local/filter_reads.nf, length 500-2,000,
                quality filter off unless min_quality is given)
  INDEX_DB      minimap2 2.24 -x map-ont -k 21 -d (modules/local/index_db.nf)
  ALIGN_READS   minimap2 2.24 -ax map-ont -N <N> -p 0.9 -K 500M -k 21 | samtools 1.18 view -b
  LIKELIHOOD    bin/read_likelihood.py
  RUN_VI        bin/estimate_abundance.py
The Python steps run in the benchmark environment, which has the image's package
versions (containers/requirements.txt). Every output is written under a temporary
name (tmp_ prefix) and renamed when complete, so an interrupted run never leaves a partial file behind,
and a step whose output exists is skipped. Each step writes <output>.time with its
wall time in seconds.
"""
import gzip
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from common import BENCH, BIN

PYTHON = os.environ.get("BENCH_PYTHON", sys.executable)
FASTPLONG_IMAGE = {
    "0.2.2": "quay.io/biocontainers/fastplong:0.2.2--heae3180_0",
    "0.7.0": "quay.io/biocontainers/fastplong:0.7.0--h43da1c4_0",
}
ALIGN_IMAGE = {
    # minimap2 2.24 + samtools 1.18, the pipeline's 'alignment' label
    "2.24": "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:365b17b986c1a60c1b82c6066a9345f38317b763-0",
    # minimap2 2.31 alone (D16 arm); samtools from the 2.24 image converts its SAM
    "2.31": "quay.io/biocontainers/minimap2:2.31--h118bc1c_0",
}
UID = f"{os.getuid()}:{os.getgid()}"


def docker(image, cmd, mounts, cpus=None):
    # The label lets an interrupted benchmark stop its own containers
    # (docker ps -q --filter label=nanovi-phase7a | xargs docker stop)
    args = ["docker", "run", "--rm", "-u", UID, "--network", "none", "--label", "nanovi-phase7a"]
    if cpus:
        args += ["--cpus", str(cpus)]
    for m in sorted(set(mounts)):
        args += ["-v", f"{m}:{m}"]
    return args + [image, "bash", "-euo", "pipefail", "-c", cmd]


def _run(out, argv, **kwargs):
    out = Path(out)
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name("tmp_" + out.name)  # a prefix keeps the extension, which some tools read
    start = time.time()
    with open(str(out) + ".log", "w") as log:
        subprocess.run([a.replace("{OUT}", str(tmp)) for a in argv], stdout=log, stderr=subprocess.STDOUT,
                       check=True, **kwargs)
    Path(str(out) + ".time").write_text(f"{time.time() - start:.2f}\n")
    tmp.rename(out)
    return out


def filter_reads(fastq, out, min_quality=None, version="0.2.2", threads=2):
    out = Path(out)
    quality = (f"--mean_qual {min_quality} --unqualified_percent_limit 100" if min_quality is not None
               else "--disable_quality_filtering")
    report = out.with_name(out.name.replace(".fastq.gz", ""))
    cmd = (f"fastplong -i {fastq} -o {{OUT}} --thread {threads} --html {report}.html --json {report}.json "
           f"--length_required 500 --length_limit 2000 {quality}")
    return _run(out, docker(FASTPLONG_IMAGE[version], cmd, [Path(fastq).parent, out.parent]))


def build_index(db_dir, out, threads=8, version="2.24", k=21):
    out = Path(out)
    cmd = f"minimap2 -x map-ont -k {k} -t {threads} -d {{OUT}} {db_dir}/species_taxid.fasta"
    return _run(out, docker(ALIGN_IMAGE[version], cmd, [Path(db_dir), out.parent]))


def align(reads, index, out, N, threads=4, version="2.24"):
    out = Path(out)
    if out.exists():  # the 2.31 path deletes its SAM, which _run would otherwise redo
        return out
    mm2 = (f"minimap2 -ax map-ont -t {threads} -N {N} -p 0.9 -K 500M -k 21 --split-prefix {out}.temp "
           f"{index} {reads}")
    if version == "2.24":
        cmd = f"{mm2} | samtools view -b -@ 2 -o {{OUT}} -"
        return _run(out, docker(ALIGN_IMAGE[version], cmd, [Path(reads).parent, Path(index).parent, out.parent]))
    # minimap2 2.31 writes SAM; the pipeline's samtools converts it in a second container
    sam = out.with_suffix(".sam")
    _run(sam, docker(ALIGN_IMAGE[version], f"{mm2} > {{OUT}}",
                     [Path(reads).parent, Path(index).parent, out.parent]))
    _run(out, docker(ALIGN_IMAGE["2.24"], f"samtools view -b -@ 2 -o {{OUT}} {sam}", [out.parent]))
    sam.unlink()
    return out


def subset_bam(bam, reads, out):
    """Alignments of the reads in a FASTQ, taken from a BAM of a superset of those reads.

    minimap2 maps each read independently, so this equals aligning the subset (checked:
    identical records for 250 of 500 reads on the reps index, and for 100 reads at -t 1
    against a 10k-read -t 4 alignment on the exact index). Used for filters that only
    drop reads (--min_quality).
    """
    out = Path(out)
    if out.exists():
        return out
    names = out.with_name(out.name + ".names")
    with gzip.open(reads, "rt") as fq:
        names.write_text("".join(line[1:].split()[0] + "\n" for i, line in enumerate(fq) if i % 4 == 0))
    _run(out, docker(ALIGN_IMAGE["2.24"], f"samtools view -b -N {names} -o {{OUT}} {bam}",
                     [Path(bam).parent, out.parent]))
    names.unlink()
    return out


def likelihood(bam, out, error_profile="primary", score="emu", min_identity=None, min_aligned_length=None,
               unknown_delta=None):
    # Only options that differ from the defaults are passed: Phase 7b removed the
    # experimental ones from bin/, so those arms need bin/ at commit 34d0411, while
    # the reference arm runs with the current bin/ too
    out = Path(out)
    argv = [PYTHON, str(BIN / "read_likelihood.py"), "--alignments", str(bam), "--scores", "{OUT}",
            "--counts", str(out.with_suffix(".counts.json"))]
    for flag, value in (("--error-profile", None if error_profile == "primary" else error_profile),
                        ("--score", None if score == "emu" else score),
                        ("--min-identity", min_identity), ("--min-aligned-length", min_aligned_length),
                        ("--unknown-delta", unknown_delta)):
        if value is not None:
            argv += [flag, str(value)]
    return _run(out, argv)


def estimate(npz, out, estimator="vb", min_abundance=0.001, no_pruning=False, unknown=False,
             credible_intervals=False):
    argv = [PYTHON, str(BIN / "estimate_abundance.py"), "--scores", str(npz), "--output", "{OUT}",
            "--estimator", estimator, "--min-abundance", str(min_abundance)]
    if no_pruning:
        argv.append("--no-pruning")
    if unknown:
        argv.append("--unknown-component")
    if credible_intervals:
        argv.append("--credible-intervals")
    if not Path(out).exists() and read_counts(npz)["assigned_count"] == 0:
        # a filter that keeps no reads: the pipeline stops here; record an empty estimate
        # (tmp_ and rename as in _run, with a .time file as every other step has)
        tmp = Path(out).with_name("tmp_" + Path(out).name)
        tmp.write_text(json.dumps({"abundance": {}, "below_detection_count": 0}))
        Path(str(out) + ".time").write_text("0.00\n")
        tmp.rename(out)
        return Path(out)
    return _run(out, argv, env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                                "MKL_NUM_THREADS": "1"})


def read_counts(npz):
    return json.loads(Path(npz).with_suffix(".counts.json").read_text())


def index_path(db, version="2.24", k=21):
    suffix = "" if version == "2.24" else f".minimap2-{version}"
    return BENCH / "indexes" / f"{db}.map-ont.k{k}{suffix}.mmi"
