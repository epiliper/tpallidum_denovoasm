process CONTIGUITY_GENE_MASK {
        tag "${meta.id}"
        label "process_single"
        container "quay.io/gdsc/biopython-pysam:3.12.1"

        input:
        tuple val(meta), path(fasta), path(bam), path(bai)
        path(instr_sheet)

        output:
        tuple val(meta), path("*genemasked*.fa*"), emit: fasta
        tuple val(meta), path("*genemask.tsv"), emit: masklog

        script:
        """
        check_contig_gene_coverage.py --instructions ${instr_sheet} \\
            --bam ${bam} \\
            --fasta ${fasta} \\
            --output_prefix ${meta.id}
        """

    }
