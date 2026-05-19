process POLISH {
    tag "${meta.id}"
    label 'process_single'
    container "quay.io/epil02/mafft_python310:0.02"

    input: 
    tuple val(meta), path(denovo_fasta), path(chosen_ref)

    output:
    tuple val(meta), path("*polished.fa")

    script:
    def prefix = task.ext.prefix ?: meta.id

    """
    polish.py --query ${denovo_fasta} --reference ${chosen_ref} --header_name ${prefix} --output ${prefix}_polished.fa
    """
    }
