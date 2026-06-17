process TRIM_END_GLUE {
    label "process_single"
    tag "${meta.id}"
    container "quay.io/biocontainers/sed:4.9"

    input:
    tuple val(meta), path(fasta)
    val(n_remove)

    output:
    tuple val(meta), path("${meta.id}_unglued.fasta"), emit: fasta

    script:

    def out = "${meta.id}_unglued.fasta"

    """
    head -n 1 ${fasta} > ${out}
    tail -n +2 ${fasta} | sed -E 's/^.{$n_remove}|.{$n_remove}\$//g' >> ${out}
    """
    }
