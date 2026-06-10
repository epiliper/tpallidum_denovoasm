process MINIMAP2_IDX_ALN_ASM {
    tag "${meta.id}"
    label "process_high"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data' :
            'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce' }"

    input:
    tuple val(meta), path(reads), path(ref)
    val suffix

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai

    script:

    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    minimap2 \\
        ${args} \\
        -ax asm5 -t ${task.cpus} \\
        ${ref} \\
        ${reads} | \\
        samtools view -h -F 2308 | \\
        samtools sort -@ ${task.cpus} -o ${prefix}_${suffix}.bam -

    samtools index ${prefix}_${suffix}.bam -@ ${task.cpus}
    """
    }
