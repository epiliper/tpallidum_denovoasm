process SPLIT_CONTIGS {
        tag "${meta.id}"
        label "process_single"
        container "quay.io/epil02/biopython_pandas_numpy_regex:0.02"

        input:
        tuple val(meta), path(contigs), path(ref)

        output:
        tuple val(meta), path("*split*.fa*"), emit: contigs

        script:
        """
        handle_contig_termini.py --contigs ${contigs} --ref ${ref} \\
            --output_prefix ${meta.id}
        """
    }
