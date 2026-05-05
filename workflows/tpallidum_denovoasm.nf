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

include { BOWTIE2_BUILD } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_GNA } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_TRNA } from '../modules/nf-core/bowtie2/align/main'

include { BBMAP_BBDUK as BBDUK_REMOVE } from '../modules/nf-core/bbmap/bbduk/main'

include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_COMBINE_GNA_TRNA } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FASTQ } from '../modules/nf-core/samtools/fastq/main'

include { PICARD_MARKDUPLICATES } from '../modules/nf-core/picard/markduplicates/main'

include { MEGAHIT } from '../modules/nf-core/megahit/main'

include { KMA_KMA } from '../modules/nf-core/kma/kma/main'
include { KMA_INDEX } from '../modules/nf-core/kma/index/main'
include { SELECT_BEST_KMA_REF } from '../modules/local/select_best_kma_ref'

include { BWA_INDEX as BWA_INDEX_DB_REF } from '../modules/nf-core/bwa/index/main'
include { BWA_INDEX as BWA_INDEX_ROUND1 } from '../modules/nf-core/bwa/index/main'
include { BWA_INDEX as BWA_INDEX_ROUND2 } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_MEM_ALIGN_TO_DB_REF } from '../modules/nf-core/bwa/mem/main'
include { BWA_MEM as BWA_MEM_ALIGN_ROUND1 } from '../modules/nf-core/bwa/mem/main'
include { BWA_MEM as BWA_MEM_ALIGN_ROUND2 } from '../modules/nf-core/bwa/mem/main'

include { CREATE_SCAFFOLD } from '../modules/local/create_scaffold'

include { IVAR_CONSENSUS as IVAR_CONSENSUS_ROUND1 } from '../modules/nf-core/ivar/consensus/main'
include { IVAR_CONSENSUS as IVAR_CONSENSUS_ROUND2 } from '../modules/nf-core/ivar/consensus/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def bwa_mem_aln_channel(reads_ch, index_ch, fasta_ch) {
    return reads_ch
        .join(index_ch)
        .join(fasta_ch)
        .multiMap { meta, reads, index, fasta ->
            reads: [meta, reads]
            index: [meta, index]
            fasta: [meta, fasta]
        }
}

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

    // na_channels.gna.view{ it -> "GNA: ${it}" }
    // na_channels.trna.view{ it -> "TRNA ${it}" }
    // na_channels.rrna.view{ it -> "RRNA: ${it}" }


    // FASTQC(ch_samplesheet)
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, file -> file })


    ///////////////////////////////
    /// BEGIN: CORE WORKFLOW EP ///

    /////////////////////////////////////////////////////////////////////
    // PART 1: get TP-mapping reads, align gNA and tRNA, dedup and filter
    /////////////////////////////////////////////////////////////////////

    channel.value(file(ref_fasta)).map { fasta_file -> [[id: fasta_file.baseName], fasta_file ] }.set { ref_ch }
    SAMTOOLS_FAIDX(ref_ch, false)
    ref_ch.join(SAMTOOLS_FAIDX.out.fai).set { ref_ch_2 }

    BOWTIE2_BUILD(ref_ch_2)

    // don't align rrna
    BOWTIE2_ALIGN_GNA(na_channels.gna, BOWTIE2_BUILD.out.index, ref_ch_2, "GNA", false, true)
    BOWTIE2_ALIGN_TRNA(na_channels.trna, BOWTIE2_BUILD.out.index, ref_ch_2, "TRNA", false, true)

    // join by meta here
    SAMTOOLS_MERGE_COMBINE_GNA_TRNA(BOWTIE2_ALIGN_GNA.out.bam.join(BOWTIE2_ALIGN_TRNA.out.bam))

    PICARD_MARKDUPLICATES(SAMTOOLS_MERGE_COMBINE_GNA_TRNA.out.bam, ref_ch_2)

    // convert back to fastq for subsequent steps
    SAMTOOLS_FASTQ(PICARD_MARKDUPLICATES.out.bam, false)
    SAMTOOLS_FASTQ.out.fastq.set { fastq_ch }

    BBDUK_REMOVE(fastq_ch, bbduk_remove)
    BBDUK_REMOVE.out.reads.set { fastq_ch }

    ///////////////////////////////////////////////
    // PART 2: DENOVO ASSEMBLY, REFERENCE SELECTION
    ///////////////////////////////////////////////

    MEGAHIT(fastq_ch)
    MEGAHIT.out.contigs.set { contigs_ch }

    KMA_INDEX(ref_ch_2.map { meta, fasta, _fai -> [meta, fasta]})
    KMA_KMA(contigs_ch, KMA_INDEX.out.index)

    SELECT_BEST_KMA_REF(KMA_KMA.out.res, ref_ch_2.map { _meta, fasta, _fai -> fasta}).set { chosen_ref_ch }

    /////////////////////////////////////
    // PART 3: GENERATE INITIAL CONSENSUS
    ////////////////////////////////////

    BWA_INDEX_DB_REF(chosen_ref_ch)

    bwa_mem_aln_channel(
        contigs_ch,
        BWA_INDEX_DB_REF.out.index,
        chosen_ref_ch,
    ).set { aln_channels }

    // contigs_ch
    //     .join(BWA_INDEX_DB_REF.out.index)
    //     .join(chosen_ref_ch)
    //     .multiMap { meta, contigs, index, fasta -> 
    //         reads: [ meta, contigs ]
    //         index: [ meta, index   ] 
    //         fasta: [ meta, fasta   ] 
    //         }
    //     .set { aln_channels }

    BWA_MEM_ALIGN_TO_DB_REF(aln_channels.reads, aln_channels.index, aln_channels.fasta, true)

    CREATE_SCAFFOLD(BWA_MEM_ALIGN_TO_DB_REF.out.bam.join(chosen_ref_ch), min_contig_len_bp)

    //////////////////////////////
    // PART 4: ITERATIVE ALIGNMENT
    /////////////////////////////

    BWA_INDEX_ROUND1(CREATE_SCAFFOLD.out.scaffold)
    bwa_mem_aln_channel(fastq_ch, BWA_INDEX_ROUND1.out.index, chosen_ref_ch).set { aln_channels }
    BWA_MEM_ALIGN_ROUND1(aln_channels.reads, aln_channels.index, aln_channels.fasta, true)

    BWA_MEM_ALIGN_ROUND1.out.bam.join(chosen_ref_ch).multiMap{ meta, bam, chosen_ref -> 
        bam_ch: [ meta, bam ]
        ref_ch: [ meta, chosen_ref ]
    }.set { ivar_input }

    IVAR_CONSENSUS_ROUND1(ivar_input.bam_ch, ivar_input.ref_ch, false)

    BWA_INDEX_ROUND2(IVAR_CONSENSUS_ROUND1.out.fasta)
    bwa_mem_aln_channel(fastq_ch, BWA_INDEX_ROUND2.out.index, chosen_ref_ch).set { aln_channels }
    BWA_MEM_ALIGN_ROUND2(aln_channels.reads, aln_channels.index, aln_channels.fasta, true)

    BWA_MEM_ALIGN_ROUND2.out.bam.join(IVAR_CONSENSUS_ROUND1.out.fasta).multiMap{ meta, bam, chosen_ref -> 
        bam_ch: [ meta, bam ]
        ref_ch: [ meta, chosen_ref ]
    }.set { ivar_input }


    IVAR_CONSENSUS_ROUND2(ivar_input.bam_ch, ivar_input.ref_ch, false)

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
