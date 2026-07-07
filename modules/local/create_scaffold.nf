process CREATE_SCAFFOLD {
    tag "$meta.id"
    label 'process_medium'
    container "quay.io/michellejlin/tpallidum_wgs:latest"

    // NOTE: bam is assumed to be sorted
    input:
    tuple val(meta), path(bam), path(ref)
    val(min_contig_len)

    output:
    tuple val(meta), path("*initial_consensus*.fasta"), emit: scaffold
    tuple val(meta), path("*_contig_stats.tsv"), emit: contig_stats

    script:
    def prefix = "${meta.id}"

    """
    tp_make_seq.R $prefix $bam $ref $min_contig_len
    """
    }
