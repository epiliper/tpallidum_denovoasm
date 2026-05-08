process MASK_TP_FASTA {
        tag "${meta.id}"
        label "process_single"
        container "quay.io/epil02/biopython_pandas_numpy_regex:0.02"

        input: 
        tuple val(meta), path(fasta)

        output:
        tuple val(meta), path("*_masked_NCBI.fasta"), emit: fasta

        script:
        """
        mask_NCBI.py -f $fasta
        """
    }
