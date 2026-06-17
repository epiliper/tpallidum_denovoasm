process MERGE_FASTQ {
    tag "$meta.id"
    label "process_medium"
    container "quay.io/biocontainers/pigz:2.8"

    input:
    tuple val(meta), path(fq1), path(fq2), path(fq3)

    output:
    tuple val(meta), path("*MERGED*.fastq.gz"), emit: fastq

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end)
        """
        cat ${fq1} ${fq2} ${fq3} > ${prefix}_MERGED.fastq.gz
        """
    else
        """
        cat ${fq1[0]} ${fq2[0]} ${fq3[0]} > ${prefix}_MERGED_1.fastq.gz
        cat ${fq1[1]} ${fq2[1]} ${fq3[1]} > ${prefix}_MERGED_2.fastq.gz
        """
}
