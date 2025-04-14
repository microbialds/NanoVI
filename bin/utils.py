# utils.py

import time
import os
import subprocess
import pathlib
import pysam
from Bio import SeqIO

def timed_function(name, func, *args, **kwargs):
    """
    Print start/stop times for the wrapped function and return its result.
    """
    print(f"Starting {name}...")
    start_time = time.time()
    result = func(*args, **kwargs)
    elapsed_time = time.time() - start_time
    print(f"✅ {name} finished. Time used: {elapsed_time:.2f} seconds\n")
    return result

def filter_reads(input_fastq, output_fastq, min_length, max_length, output_dir):
    """
    Runs fastplong to filter reads, compresses the output FASTQ and reports (.html, .json).
    Returns path to the compressed FASTQ (.gz).
    """
    print("Starting read filtering with fastplong...")
    start_time = time.time()

    os.makedirs(output_dir, exist_ok=True)
    sample_name = os.path.splitext(os.path.basename(input_fastq))[0]
    fastplong_output_dir = os.path.join(output_dir, "fastplong")
    os.makedirs(fastplong_output_dir, exist_ok=True)

    cmd_fastplong = (
        f"fastplong "
        f"-i {input_fastq} "
        f"-o {output_fastq} "
        f"--html {fastplong_output_dir}/{sample_name}_report.html "
        f"--json {fastplong_output_dir}/{sample_name}_report.json "
        f"--length_required {min_length} "
        f"--length_limit {max_length}"
    )
    subprocess.run(cmd_fastplong, shell=True, check=True)

    # Compress the filtered FASTQ
    subprocess.run(f"gzip -f {output_fastq}", shell=True, check=True)
    # Compress HTML and JSON reports as well
    subprocess.run(f"gzip -f {fastplong_output_dir}/{sample_name}_report.html", shell=True, check=True)
    subprocess.run(f"gzip -f {fastplong_output_dir}/{sample_name}_report.json", shell=True, check=True)

    elapsed_time = time.time() - start_time
    print(f"fastplong finished. Time used: {elapsed_time:.2f} seconds\n")

    return f"{output_fastq}.gz"

def validate_input(path):
    """
    Checks if 'path' is a valid FASTA/FASTQ or a SAM/BAM. Raises TypeError otherwise.
    """
    sam_pysam = None
    try:
        sam_pysam = pysam.AlignmentFile(path)
    except (ValueError, OSError):
        pass
    if sam_pysam:
        return

    # If it's not SAM/BAM, try FASTA or FASTQ
    fasta_rd, fastq_rd = None, None
    try:
        fasta_rd = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
        fastq_rd = SeqIO.to_dict(SeqIO.parse(path, "fastq"))
    except (UnicodeDecodeError, ValueError):
        pass

    if not (fasta_rd or fastq_rd):
        raise TypeError("Input must be in desired format: fasta or fastq (or SAM/BAM)")

def gather_fastq_paths(input_paths):
    """
    Given a list of file or directory paths, return a sorted list of valid FASTQ paths
    (.fastq, .fastq.gz, .fq, .fq.gz). If a directory is provided, it recursively searches.
    """
    valid_exts = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}
    collected_paths = []

    for p in input_paths:
        p = os.path.abspath(p)
        if os.path.isdir(p):
            # Recursively walk through the directory
            for root, dirs, files in os.walk(p):
                for f in files:
                    fp = os.path.join(root, f)
                    if any(fp.endswith(ext) for ext in valid_exts):
                        collected_paths.append(os.path.abspath(fp))
        else:
            # It's a file, check extension
            if any(p.endswith(ext) for ext in valid_exts):
                collected_paths.append(os.path.abspath(p))

    collected_paths.sort()
    return collected_paths