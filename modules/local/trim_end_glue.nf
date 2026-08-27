process TRIM_END_GLUE {
    label "process_single"
    tag "${meta.id}"
    container "quay.io/epil02/mafft_biopython:3.10"

    input:
    tuple val(meta), path(query), path(chosen_ref)
    val(n_remove)

    output:
    tuple val(meta), path("${meta.id}_end_trimmed.fasta"), emit: fasta

    script:

    def out = "${meta.id}_end_trimmed.fasta"

    """
    cat ${chosen_ref} ${query} > ${meta.id}_mafft_in
    mafft --op 2 --ep 0.25 ${meta.id}_mafft_in > ${meta.id}_aln.fa
    trim_query_overhangs.py ${meta.id}_aln.fa > ${out}
    """
    }
