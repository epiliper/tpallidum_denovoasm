process MERGE_FASTQ {
       tag "$meta.id" 
       label "process_medium"
       container "quay.io/biocontainers/pigz:2.8"

       input:
       tuple val(meta), val(reads)

       output:
       tuple val(meta), path("*MERGED*.fastq.gz"), emit: fastq

       script:
       def prefix = task.ext.prefx ?: "${meta.id}"
       def in1 = meta.single_end ? reads.se : reads.r1.join(' ')
       def in2 = meta.single_end ? "" : reads.r2.join(' ')

       def merge_cmd = meta.single_end ? "cat $in1 > ${prefix}_MERGED.fastq.gz" : "cat ${in1} > ${prefix}_MERGED_1.fastq.gz && cat ${in2} > ${prefix}_MERGED_2.fastq.gz"

       """
       ${merge_cmd}
       """

    }
