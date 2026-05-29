process DENOVO_ASSEMBLE {
    // container "quay.io/biocontainers/unicycler:0.4.4--py37h13b99d1_3"
    container "quay.io/biocontainers/unicycler:0.5.1--py312hdcc493e_5 "
    label 'process_high'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(reads)

    output:
        tuple val(meta),file("${meta.id}_assembly.gfa"),file("${meta.id}_assembly.fasta")// into Unicycler_ch
        tuple val(meta), file("*_assembly.gfa"), emit: gfa
        tuple val(meta), file("*assembly.fasta"), emit: fasta
        file("*")// into Unicycler_dump_ch

    script:

    def run = meta.single_end ? "unicycler -s $reads" : "unicycler -1 ${reads[0]} -2 ${reads[1]}"
    """
    #!/bin/bash

    ${run} -o ./ -t ${task.cpus}

    cp assembly.gfa ${meta.id}_assembly.gfa
    cp assembly.fasta ${meta.id}_assembly.fasta

    """
}
