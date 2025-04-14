# alignment.py

import os
import subprocess
import pathlib
import time

def generate_alignments(in_file_list, out_basename, database, args):
    """
    Uses a pre-indexed Minimap2 database if available, otherwise generates it.
    Then runs Minimap2 alignment on the input reads.

    in_file_list: List of input FASTQ/FASTA files (already filtered).
    out_basename: Base filepath for alignment outputs.
    database: Path to Emu database directory (containing species_taxid.fasta).
    args: Argument object with fields:
        - threads, kmer_size, N, K, type, split_prefix
    """
    input_file = " ".join(in_file_list)
    filetype = pathlib.PurePath(input_file).suffix

    sam_align_file = f"{out_basename}_vi_alignments.sam"
    
    db_sequence_file = os.path.join(database, 'species_taxid.fasta')
    db_index_file = os.path.join(database, 'gtdb_index.mmi')  # Pre-indexed file

    # If index doesn’t exist, build it with the correct k-mer size
    if not os.path.exists(db_index_file):
        print("🔄 Indexing database... (This happens only once)")
        index_cmd = (
            f"minimap2 -k {args.kmer_size} -d {db_index_file} {db_sequence_file}"
        )
        subprocess.run(index_cmd, shell=True, check=True)
        print(f"✅ Pre-indexed database saved as {db_index_file}")

    # Construct minimap2 command
    minimap_cmd = (
        f"minimap2 -ax {args.type} "
        f"-t {args.threads} "
        f"-N {args.N} "
        f"-p .9 "
        f"-K {args.K} "
        f"-k {args.kmer_size} "
        f"--split-prefix {args.split_prefix} "
        f"{db_index_file} {input_file} "
        f"-o {sam_align_file}"
    )

    # Print a start message for alignment
    print("\nStarting minimap2 alignment...")

    start_time = time.time()
    # Run minimap2; if you want live logs, remove .check_output() in favor of .run()
    subprocess.check_output(minimap_cmd, shell=True)
    end_time = time.time()

    elapsed = end_time - start_time
    print(f"✅ Minimap2 alignment finished. Time used: {elapsed:.2f} seconds\n")

    return sam_align_file
