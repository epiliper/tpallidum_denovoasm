/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                                            } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                            } from '../subworkflows/local/utils_nfcore_tpallidum_denovoasm_pipeline'

include { BOWTIE2_BUILD                                     } from '../modules/nf-core/bowtie2/build/main'       
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_GNA                } from '../modules/nf-core/bowtie2/align/main'       
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_TRNA               } from '../modules/nf-core/bowtie2/align/main'       

include { BBMAP_BBDUK as BBDUK_REMOVE                       } from '../modules/nf-core/bbmap/bbduk/main'

include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_COMBINE_GNA_TRNA } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_FAIDX                                    } from '../modules/nf-core/samtools/faidx/main' 
include { SAMTOOLS_FASTQ } from '../modules/nf-core/samtools/fastq/main'

include { PICARD_MARKDUPLICATES                             } from '../modules/nf-core/picard/markduplicates/main'

include { MEGAHIT                                           } from '../modules/nf-core/megahit/main'

include { KMA_KMA                                           } from '../modules/nf-core/kma/kma/main'
include { KMA_INDEX                                         } from '../modules/nf-core/kma/index/main'
include { SELECT_BEST_KMA_REF                               } from '../modules/local/select_best_kma_ref'

 include { BWA_INDEX                                        } from '../modules/nf-core/bwa/index/main'
 include { BWA_MEM as BWA_MEM_ALIGN_TO_DB_REF               } from '../modules/nf-core/bwa/mem/main'
 include { BWA_MEM as BWA_MEM_ALIGN_ROUND1                  } from '../modules/nf-core/bwa/mem/main'
 include { BWA_MEM as BWA_MEM_ALIGN_ROUND2                  } from '../modules/nf-core/bwa/mem/main'
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
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    // split up into our NAs of choice, for details on how we parse out the different NA types from the samplesheet, see PIPELINE_INITIALISATION
    ch_samplesheet.multiMap{ meta, gna, trna, rrna -> 
        gna: [meta, gna]
        trna: [meta, trna]
        rrna: [meta, rrna]
        }.set { na_channels }

    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })


    ///////////////////////////////
    /// BEGIN: CORE WORKFLOW EP ///

    // PART 1: get TP-mapping reads, align gNA and tRNA, dedup and filter
    /////////////////////////////////////////////////////////////////////

    ref_fasta.map{ fasta_file -> [null, fasta_file, null] }.set { ref_ch }
    
    SAMTOOLS_FAIDX(ref_ch, null)
    SAMTOOLS_FAIDX.out.fa.join(SAMTOOLS_FAIDX.out.fai).set { ref_ch }

    BOWTIE2_BUILD(ref_ch)

    // don't align rrna
    BOWTIE2_ALIGN_GNA(na_channels.gna, BOWTIE2_BUILD.out.index, null, false, true)
    BOWTIE2_ALIGN_TRNA(na_channels.trna, BOWTIE2_BUILD.out.index, null, false, true)

    // join by meta here
    SAMTOOLS_MERGE_COMBINE_GNA_TRNA(BOWTIE2_ALIGN_GNA.out.bam.join(BOWTIE2_ALIGN_TRNA.out.bam), null)

    PICARD_MARKDUPLICATES(SAMTOOLS_MERGE_COMBINE_GNA_TRNA.out.bam, ref_ch)

    // convert back to fastq for subsequent steps
    SAMTOOLS_FASTQ(PICARD_MARKDUPLICATES.out.bam, false).out.set { fastq_ch }

    BBDUK_REMOVE(fastq_ch, bbduk_remove).out.reads.set { fastq_ch }

    // PART 2: DENOVO ASSEMBLY, REFERENCE SELECTION
    //////////////////////////

    MEGAHIT(fastq_ch).out.contigs.set { contigs_ch }

    KMA_INDEX(ref_fasta)

    KMA_KMA(contigs_ch, KMA_INDEX.out.index)

    SELECT_BEST_KMA_REF(KMA_KMA.out.res)

    // PART 3: GENERATE INITIAL CONSENSUS
    ////////////////////////////

    BWA_INDEX(SELECT_BEST_KMA_REF.out.chosen_ref)

    contigs_ch
        .join(BWA_INDEX.out.index)
        .join(SELECT_BEST_KMA_REF.out.chosen_ref)
        .multiMap { meta, contigs, index, fasta -> 
            reads: [ meta, contigs ]
            index: [ meta, index   ] 
            fasta: [ meta, fasta   ] 
            }
        .set { aln_channels }

    BWA_MEM_ALIGN_TO_DB_REF(aln_channels.contigs, aln_channels.index, aln_channels.fasta, true)

    // PART 4: ITERATIVE ALIGNMENT

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
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'tpallidum_denovoasm_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'tpallidum_denovoasm'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
