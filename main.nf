#!/usr/bin/env nextflow
include { TPALLIDUM_DENOVOASM } from './workflows/tpallidum_denovoasm'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_tpallidum_denovoasm_pipeline'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_tpallidum_denovoasm_pipeline'

workflow {
    main:

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
        params.n_end_glue_bases,
        params.outdir
    )

    // emit:
    // multiqc_report = TPALLIDUM_DENOVOASM.out.multiqc_report // channel: /path/to/multiqc_report.html
}
