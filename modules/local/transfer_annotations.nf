process TRANSFER_ANNOTATIONS {
    tag "${meta.id}"
    label "process_medium"
    container "quay.io/epil02/mafft_r:0.01"

    input:
    tuple val(meta), path(input_fasta), path(ref_fasta)
    path main_coord_xlsx
    path coord_lookups

    output:
    tuple val(meta), path("*extracted_ORFs.csv")
    tuple val(meta), path("*annot_aln.fa")

    script:
    // annotation_transfer.R ${meta.id}_aln ${coord_xlsx} ${ref_fasta} ${meta.id}
    """
    cat ${ref_fasta} ${input_fasta} >> ${meta.id}_mafft_in
    mafft --op 2 --ep 0.25 ${meta.id}_mafft_in > ${meta.id}_annot_aln.fa
    annotation_transfer_ANY_REF.R ${meta.id}_annot_aln.fa ${main_coord_xlsx} ${input_fasta} ${meta.id}
    """
}
