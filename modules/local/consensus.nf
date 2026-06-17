process CONSENSUS {
        tag "${meta.id}"
        label "process_medium"
        container "quay.io/michellejlin/tpallidum_wgs"

        input:
        tuple val(meta), path(bam), path(bai), path(ref)
        val(suffix)
        val(min_depth)

        output:
        tuple val(meta), path("*consensus.fasta"), emit: fasta

        script:
        def prefix = task.ext.prefix ?: meta.id
        def sampleprefix = "${prefix}_${suffix}"

        """
        tp_make_consensus.R ${sampleprefix} ${bam} ${ref} ${min_depth}
        """
    }
