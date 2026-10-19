nextflow.enable.dsl = 2

params.outdir   = "results"
params.adapter  = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.genome   = "reference/GRCh38_clean.fa"
params.mature   = "reference/mature_hsa_clean.fa"
params.hairpin  = "reference/hairpin_hsa_clean.fa"
params.metadata = "metadata.csv"
params.bowtie_index = "${projectDir}/results/bowtie_index/hg38_bowtie_index"
//---------------- PREPROCESSING ----------------//

include { QUALITY_CHECK } from './modules/preprocessing/01_quality_check.nf'
include { TRIM } from './modules/preprocessing/02_trimming.nf'
include { POST_QC } from './modules/preprocessing/03_post_qc.nf'
include { BOWTIE_ALIGN } from './modules/preprocessing/bowtie_align.nf'

//---------------- miRDeep2 ----------------//

include { MIRDEEP_ALIGN } from './modules/mirdeep2/04_mirdeep_align.nf'
include { MIRDEEP_PREDICT } from './modules/mirdeep2/05_mirdeep_predict.nf'
include { MIRDEEP_QUANTIFICATION } from './modules/mirdeep2/06_mirdeep_quantify.nf'

//---------------- Known miRNA ----------------//

include { PREPARE_KNOWN_COUNT_MATRIX } from './modules/known_mirna/07_known_count_matrix.nf'
include { KNOWN_DESEQ2 } from './modules/known_mirna/08_known_DEA.nf'
include { KNOWN_TARGET_PREDICT } from './modules/known_mirna/08_known_target_prediction.nf'
include { KNOWN_NETWORK_ANALYSIS } from './modules/known_mirna/09_known_network_analysis.nf'
include { KNOWN_ENRICHMENT } from './modules/known_mirna/10_known_enrichment.nf'

//---------------- MultiQC ----------------//

include { MULTIQC } from './modules/multiqc.nf'

workflow {

    //---------------- Input FASTQ ----------------//

    read_ch = Channel
        .fromPath("data/*.fastq")
        .map { file ->
            tuple(file.baseName, file)
        }

    //---------------- QC ----------------//

    QUALITY_CHECK(read_ch)

    //---------------- Trimming ----------------//

    trimmed_ch = TRIM(read_ch)

    //---------------- Post-trim QC ----------------//

    POST_QC(trimmed_ch.trimmed)

    //---------------- Bowtie Alignment ----------------//

    bowtie_result = BOWTIE_ALIGN(
        trimmed_ch.trimmed
    )

    bowtie_logs = bowtie_result.map { sample, sam, log ->
    log
}

    //---------------- miRDeep2 Alignment ----------------//

    aligned = MIRDEEP_ALIGN(
        trimmed_ch.trimmed
    )

    genome_ref_ch  = Channel.value(file(params.genome))
    mature_ref_ch  = Channel.value(file(params.mature))
    hairpin_ref_ch = Channel.value(file(params.hairpin))

    //---------------- Novel miRNA Prediction ----------------//

    MIRDEEP_PREDICT(
        aligned,
        genome_ref_ch,
        mature_ref_ch,
        hairpin_ref_ch
    )

    //---------------- Quantification ----------------//

    collapsed_ch = aligned.map { sample, collapsed, arf ->
        tuple(sample, collapsed)
    }

    quantification_result = MIRDEEP_QUANTIFICATION(
        collapsed_ch,
        mature_ref_ch,
        hairpin_ref_ch
    )

    expression_files = quantification_result.counts
        .map { sample, file -> file }
        .collect()

    //---------------- Count Matrix ----------------//

    count_matrix_result = PREPARE_KNOWN_COUNT_MATRIX(
        expression_files
    )

    metadata_ch = Channel.value(
        file(params.metadata)
    )

    //---------------- Differential Expression ----------------//

    deseq_script = Channel.value(
    file("${projectDir}/r_scripts/01_diff_exp_analysis.R")
)

known_deseq_result = KNOWN_DESEQ2(
    count_matrix_result.filtered_matrix,
    metadata_ch,
    deseq_script
)
    //---------------- Target Prediction ----------------//

    target_script = Channel.value(
        file("${projectDir}/r_scripts/02_target_prediction.R")
    )

    target_result = KNOWN_TARGET_PREDICT(
        known_deseq_result.early_sig,
        known_deseq_result.late_sig,
        known_deseq_result.late_early_sig,
        target_script
    )

    //---------------- Target Files ----------------//

    target_prediction_files = target_result.csv
        .flatten()
        .filter { f ->
            !f.name.endsWith("_summary.csv") &&
            !f.name.endsWith("_unique_target_genes.csv")
        }

    //---------------- Network Analysis ----------------//

    network_script = Channel.value(
        file("${projectDir}/r_scripts/03_network.R")
    )

    KNOWN_NETWORK_ANALYSIS(
        target_prediction_files.collect(),
        network_script
    )

    //---------------- Enrichment ----------------//

    KNOWN_ENRICHMENT(
        target_prediction_files.collect()
    )

    //---------------- MultiQC ----------------//

    multiqc_input = Channel.empty()

    multiqc_input = multiqc_input.mix(QUALITY_CHECK.out.zip)
    multiqc_input = multiqc_input.mix(QUALITY_CHECK.out.html)

    multiqc_input = multiqc_input.mix(TRIM.out.report)

    multiqc_input = multiqc_input.mix(POST_QC.out.zip)
    multiqc_input = multiqc_input.mix(POST_QC.out.html)

    multiqc_input = multiqc_input.mix(bowtie_logs)

    MULTIQC(
    multiqc_input.collect()
)

}
