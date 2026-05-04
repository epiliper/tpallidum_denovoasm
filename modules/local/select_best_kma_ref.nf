process SELECT_BEST_KMA_REF {
    tag "$meta.id"
    label 'process_single'

    // container

    input:
    tuple val(meta), path(kma_res), path(ref_fasta)

    output:
    tuple val(meta), path("*_best_ref.fa"), emit: chosen_ref

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    select_best_ref_from_kma.py -i $kma_res --prefix ${prefix} -r ${ref_fasta}
    """

    }
