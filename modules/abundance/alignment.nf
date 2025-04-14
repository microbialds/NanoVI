// modules/abundance/alignment.nf

process ALIGN_READS {

    label 'alignment'

    input:
    tuple val(sample_id), path(filtered_fastq)
    path species_fasta

    output:
    path "${sample_id}_vi_alignments.sam", emit: sam
    publishDir "${params.output_dir}", mode: 'copy'

    script:
    """
    if [ ! -f "gtdb_index.mmi" ]; then
        echo "🔄 Indexando base de datos con minimap2..."
        minimap2 -k ${params.kmer_size} -d gtdb_index.mmi ${species_fasta}
    fi

    echo "🧬 Alineando ${sample_id} con minimap2..."

    minimap2 -ax ${params.type} \\
             -t ${params.threads} \\
             -N ${params.N} \\
             -p 0.9 \\
             -K ${params.K} \\
             -k ${params.kmer_size} \\
             --split-prefix ${params.split_prefix} \\
             gtdb_index.mmi \\
             ${filtered_fastq} \\
             -o ${sample_id}_vi_alignments.sam
    """
}

