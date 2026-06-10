/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_tpallidum_denovoasm_pipeline'

include { BBMAP_BBDUK as BBDUK_REMOVE } from '../modules/nf-core/bbmap/bbduk/main'

include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_COMBINE_GNA_TRNA } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FASTQ } from '../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_RNA } from '../modules/nf-core/samtools/fastq/main'

include { PICARD_MARKDUPLICATES } from '../modules/nf-core/picard/markduplicates/main'
include { PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_RRNA } from '../modules/nf-core/picard/markduplicates/main'

include { MEGAHIT } from '../modules/nf-core/megahit/main'
include { DENOVO_ASSEMBLE } from '../modules/local/unicycler'

include { KMA_KMA } from '../modules/nf-core/kma/kma/main'
include { KMA_INDEX } from '../modules/nf-core/kma/index/main'
include { SELECT_BEST_KMA_REF } from '../modules/local/select_best_kma_ref'

include { MINIMAP2_IDX_ALN_ASM } from '../modules/local/minimap2'

include { BOWTIE2_BUILD as BOWTIE2_INDEX_DB_REF } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_INDEX_CHOSEN_REF } from '../modules/nf-core/bowtie2/build/main'

include { BOWTIE2_BUILD as BOWTIE2_INDEX_ROUND1 } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_INDEX_ROUND2 } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_INDEX_ROUND3 } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_INDEX_ROUND4 } from '../modules/nf-core/bowtie2/build/main'

include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_GNA} from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_TRNA } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_RRNA } from '../modules/nf-core/bowtie2/align/main'

include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_TO_DB_REF } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_ROUND1 } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_ROUND2 } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_ROUND3 } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_ROUND4 } from '../modules/nf-core/bowtie2/align/main'

include { CREATE_SCAFFOLD } from '../modules/local/create_scaffold'
include { CREATE_SCAFFOLD  as CREATE_SCAFFOLD2 } from '../modules/local/create_scaffold'

include { MERGE_FASTQ as MERGE_ALL_NA } from '../modules/local/merge_fastq'

include { CONSENSUS as CONSENSUS_ROUND1 } from '../modules/local/consensus'
include { CONSENSUS as CONSENSUS_ROUND2 } from '../modules/local/consensus'
include { CONSENSUS as CONSENSUS_ROUND3 } from '../modules/local/consensus'

include { PILON } from '../modules/nf-core/pilon/main'

include { MASK_TP_FASTA } from '../modules/local/mask_tp_fasta'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TPALLIDUM_DENOVOASM {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ref_fasta
    bbduk_remove
    min_contig_len_bp
    // multiqc_config
    // multiqc_logo
    // multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    // split up into our nucleic acids of choice. for details on how we parse out the different NA types from the samplesheet, see PIPELINE_INITIALISATION
    ch_samplesheet
        .multiMap { meta, na ->
            gna: [meta, na[0]]
            trna: [meta, na[1]]
            rrna: [meta, na[2]]
        }
        .set { na_channels }

    ///////////////////////////////
    /// BEGIN: CORE WORKFLOW EP ///

    /////////////////////////////////////////////////////////////////////
    // PART 1: get TP-mapping reads, align gNA and tRNA, dedup and filter
    /////////////////////////////////////////////////////////////////////

    channel.value([[id: file(ref_fasta).baseName], file(ref_fasta)]).set { ref_ch }

    SAMTOOLS_FAIDX(ref_ch, false)
    KMA_INDEX(ref_ch)
    KMA_INDEX.out.index.map { _meta, index -> index }.set { db_idx_kma }

    BOWTIE2_INDEX_DB_REF(ref_ch)
    BOWTIE2_INDEX_DB_REF.out.index.map { _meta, index -> index }.set { db_idx_bwa }

    // don't align rrna
    BOWTIE2_ALIGN_GNA(na_channels.gna.combine(db_idx_bwa), "GNA", false, true)
    BOWTIE2_ALIGN_TRNA(na_channels.trna.combine(db_idx_bwa), "TRNA", false, true)

    BOWTIE2_ALIGN_GNA.out.bam
        .join(BOWTIE2_ALIGN_TRNA.out.bam)
        .join(BOWTIE2_ALIGN_GNA.out.bai)
        .join(BOWTIE2_ALIGN_TRNA.out.bai)
        .map { meta, bam1, bam2, bai1, bai2 -> [ meta, [ bam1, bam2 ], [ bai1, bai2 ]]}
        .set { ch_merge_in }

    // join by meta here
    SAMTOOLS_MERGE_COMBINE_GNA_TRNA(ch_merge_in)

    PICARD_MARKDUPLICATES(
        SAMTOOLS_MERGE_COMBINE_GNA_TRNA.out.bam
        )

    // convert back to fastq for subsequent steps
    SAMTOOLS_FASTQ(
        PICARD_MARKDUPLICATES.out.bam, false)
    SAMTOOLS_FASTQ.out.fastq.set { fastq_raw_ch }

    BBDUK_REMOVE(
        fastq_raw_ch, bbduk_remove)
    BBDUK_REMOVE.out.reads.set { fastq_ch }

    ///////////////////////////////////////////////
    // PART 2: DENOVO ASSEMBLY, REFERENCE SELECTION
    ///////////////////////////////////////////////

    // MEGAHIT(fastq_ch)
    DENOVO_ASSEMBLE(fastq_ch)

    // MEGAHIT.out.contigs.set { contigs_ch }
    DENOVO_ASSEMBLE.out.fasta.set { contigs_ch }

    KMA_KMA(contigs_ch.combine(db_idx_kma))

    SELECT_BEST_KMA_REF(KMA_KMA.out.res, ref_fasta).set { chosen_ref_ch }

    /////////////////////////////////////
    // PART 3: GENERATE INITIAL CONSENSUS
    ////////////////////////////////////

    BOWTIE2_INDEX_CHOSEN_REF(chosen_ref_ch)

    MINIMAP2_IDX_ALN_ASM(contigs_ch.join(chosen_ref_ch), "contig")

    CREATE_SCAFFOLD(MINIMAP2_IDX_ALN_ASM.out.bam.join(chosen_ref_ch),
    min_contig_len_bp)

    // BWA_INDEX_CHOSEN_REF(chosen_ref_ch)
    // BWA_MEM_ALIGN_TO_CHOSEN_REF(contigs_ch.join(BWA_INDEX_CHOSEN_REF.out.index), "contig", true)

    // CREATE_SCAFFOLD(BWA_MEM_ALIGN_TO_CHOSEN_REF.out.bam.join(chosen_ref_ch), min_contig_len_bp)

    ///////////////////////////////////////////////////////
    // PART 3.5: ALIGN AND DEDUP rRNA, MERGE WITH OTHER NAS
    ///////////////////////////////////////////////////////

    BOWTIE2_ALIGN_RRNA(na_channels.rrna.join(BOWTIE2_INDEX_CHOSEN_REF.out.index), "rrna", false, true)
    PICARD_MARKDUPLICATES_RRNA(BOWTIE2_ALIGN_RRNA.out.bam)
    SAMTOOLS_FASTQ_RNA(PICARD_MARKDUPLICATES_RRNA.out.bam, false)

    // combine GA, TRNA, RRNA, and bbduk-filtered reads
    MERGE_ALL_NA(
        fastq_ch
            .join(SAMTOOLS_FASTQ_RNA.out.fastq)
            .join(BBDUK_REMOVE.out.contam)
            .map { meta, fq1, fq2, fq3 -> meta.single_end
                ? [ meta, [ se: [ fq1, fq2, fq3 ] ] ]
                : [ meta, [ r1: [ fq1[0], fq2[0], fq3[0] ], r2: [ fq1[1], fq2[1], fq3[1] ] ] ] }
            .dump(tag: "merge all input")
    )

    MERGE_ALL_NA.out.fastq.set { all_na_fastq }


    //////////////////////////////
    // PART 4: ITERATIVE ALIGNMENT
    /////////////////////////////
    BOWTIE2_INDEX_ROUND1(CREATE_SCAFFOLD.out.scaffold)
    BOWTIE2_ALIGN_ROUND1(all_na_fastq.join(BOWTIE2_INDEX_ROUND1.out.index), "to_scaffold", false, true)
    CONSENSUS_ROUND1(
        BOWTIE2_ALIGN_ROUND1.out.bam
        .join(BOWTIE2_ALIGN_ROUND1.out.bai)
        .join(CREATE_SCAFFOLD.out.scaffold),
        "intermediate"
    )

    BOWTIE2_INDEX_ROUND2(CONSENSUS_ROUND1.out.fasta)
    BOWTIE2_ALIGN_ROUND2(all_na_fastq.join(BOWTIE2_INDEX_ROUND2.out.index), "to_intermediate", false, true)
    CONSENSUS_ROUND2(
        BOWTIE2_ALIGN_ROUND2.out.bam
        .join(BOWTIE2_ALIGN_ROUND2.out.bai)
        .join(CONSENSUS_ROUND1.out.fasta),
        "final"
    )

    // BOWTIE2_INDEX_ROUND3(CONSENSUS_ROUND2.out.fasta)
    // BOWTIE2_ALIGN_ROUND3(all_na_fastq.join(BOWTIE2_INDEX_ROUND3.out.index), "pre_pileon", false, true)

    // PILON(
    //     CONSENSUS_ROUND2.out.fasta
    //     .join(BOWTIE2_ALIGN_ROUND3.out.bam)
    //     .join(BOWTIE2_ALIGN_ROUND3.out.bai),
    //     "frags"
    // )

    // BOWTIE2_INDEX_ROUND4(PILON.out.improved_assembly)
    // BOWTIE2_ALIGN_ROUND4(all_na_fastq.join(BOWTIE2_INDEX_ROUND4.out.index), "final", false, true)
    // CONSENSUS_ROUND3(
    //     BOWTIE2_ALIGN_ROUND4.out.bam
    //     .join(BOWTIE2_ALIGN_ROUND4.out.bai)
    //     .join(CONSENSUS_ROUND1.out.fasta),
    //     "complete"
    // )

    MASK_TP_FASTA(CONSENSUS_ROUND2.out.fasta)


    /// END: CORE WORKFLOW EP ///
    ////////////////////////////

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'tpallidum_denovoasm_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    // //
    // // MODULE: MultiQC
    // //
    // ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    // def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    // def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // def ch_multiqc_custom_methods_description = multiqc_methods_description
    //     ? file(multiqc_methods_description, checkIfExists: true)
    //     : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    // def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    // MULTIQC(
    //     ch_multiqc_files.flatten().collect().map { files ->
    //         [
    //             [id: 'tpallidum_denovoasm'],
    //             files,
    //             multiqc_config
    //                 ? file(multiqc_config, checkIfExists: true)
    //                 : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
    //             multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
    //             [],
    //             [],
    //         ]
    //     }
    // )

    emit:
    // multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions = ch_versions // channel: [ path(versions.yml) ]
}
