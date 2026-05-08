process MERGE_FASTQ {
       tag "$meta.id" 
       label "process_medium"
       container "quay.io/biocontainers/pigz:2.8"

       input:
       tuple val(meta), path(reads1), path(reads2)

       output:
       tuple val(meta), path("*MERGED.fastq.gz"), emit: fastq

       script:
       def prefix = task.ext.prefx ?: "${meta.id}"


       """
       cat $reads1 $reads2 >> ${prefix}_MERGED.fastq.gz
       """

    }
