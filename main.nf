#!/usr/bin/env nextflow
include { TPALLIDUM_DENOVOASM } from './workflows/tpallidum_denovoasm'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_tpallidum_denovoasm_pipeline'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_tpallidum_denovoasm_pipeline'

workflow {

    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    //
    // WORKFLOW: Run pipeline
    //
    TPALLIDUM_DENOVOASM(
        PIPELINE_INITIALISATION.out.samplesheet,
        params.refs,
        params.bbduk_filter,
        params.min_contig_len_bp,
        params.min_consensus_depth,
        params.end_glue_bp,
        params.mask_sheet,
        params.annot_master_coords,
        params.annot_coord_lookups,
        params.outdir,
    )
}
