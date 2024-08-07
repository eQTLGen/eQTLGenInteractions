#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage() {
  log.info"""
  =======================================================
    Interaction analysis v${workflow.manifest.version}
  =======================================================
    
  """.stripIndent()
}

/*
 * Parameters
 */

params.bfile = ''
params.vcf_dir = ''
params.bgen_dir = ''

params.genes_to_test = ''
params.qtls_to_test = ''
params.lab_cell_perc = ''

params.signature_matrix_name = "LM22"
params.deconvolution_method = "dtangle"
params.num_perm = 0

params.run_stratified = false
params.preadjust = false
params.cell_perc_interactions = false
params.expr_pcs = ''
params.num_expr_PCs = 20

/*
 * Channel declarations
 */
raw_expr_ch = Channel.fromPath(params.raw_expfile)
filt_exp_ch = Channel.fromPath(params.norm_expfile)
outdir_ch = Channel.fromPath(params.outdir, type: 'dir')
covars_ch = Channel.fromPath(params.covariates)
gte_ch = Channel.fromPath(params.gte)

annotation_ch = Channel.fromPath("$projectDir/data/LimixAnnotationFile.GRCh38.110.txt.gz")
gene_lengths_ch = Channel.fromPath("$projectDir/data/GeneLengths_GRCh38.110_ensg.txt.gz")

// chunk channel for the interaction analysis: chrom -> chunk
Channel
    .fromPath(params.chunk_file)
    .splitCsv( header: false )
    .map { row -> tuple(row[0].split(':')[0], row[0]) }
    .set { chunk_ch }

// genotype channel declaration. Interaction analysis doesn't run with bgen files yet
if (params.bgen_dir != '') {
  Channel.from(1..22)
  .map { chr -> tuple("$chr", file("${params.bgen_dir}/*chr${chr}.bgen"), file("${params.bgen_dir}/*chr${chr}.sample")) }
  .ifEmpty { exit 1, "Input .bgen files not found!" }
  .set { chr_bgen_pairs }
} else if (params.vcf_dir != '') {
  Channel.from(1..22)
  .map { chr -> tuple("$chr", file("${params.vcf_dir}/*chr${chr}.filtered.vcf.gz")) }
  .ifEmpty { exit 1, "Input .vcf.gz files not found!" }
  .set { chr_vcf_pairs }
  Channel.empty()
    .set { chr_bgen_pairs }
} else {
  Channel
    .from(params.bfile)
    .ifEmpty { exit 1, "Input plink prefix not found!" }
    .map { genotypes -> [file("${genotypes}.bed"), file("${genotypes}.bim"), file("${genotypes}.fam")]}
    .set { bfile_ch }
} 


include { PREPARE_COVARIATES; NormalizeExpression; ConvertVcfToBgen; ConvertVcfToPlink; MergePlinkPerChr } from './modules/prepare_data.nf'
include { RUN_INTERACTION_QTL_MAPPING; IeQTLmapping; IeQTLmapping_InteractionCovariates; SplitCovariates; PreadjustExpression } from './modules/interaction_analysis2.nf'
include { RUN_STRATIFIED_ANALYSIS; RunEqtlMappingPerGenePlink } from './modules/stratified_analysis.nf'

/* 
 * Analysis
 */

workflow {
    params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
    
    /*
     * Prepare normalized expression and covariate data
     */
    NormalizeExpression(raw_expr_ch, filt_exp_ch, params.exp_platform, gte_ch )
    norm_exp_ch = NormalizeExpression.out.norm_expression_table
    covariates_ch = PREPARE_COVARIATES(params.exp_platform, raw_expr_ch, norm_exp_ch, params.signature_matrix_name, params.deconvolution_method,covars_ch, gene_lengths_ch, annotation_ch, Channel.fromPath(params.genotype_pcs), Channel.fromPath(params.gte))
    
    /*
     * Prepare genotype data
     */
    Channel.empty()
      .set { bfile_ch }
    // if no plink genotypes provided convert VCF to plink and combine chromosomes together
    if (params.bfile == '') {
      ConvertVcfToPlink(chr_vcf_pairs)
      MergePlinkPerChr(ConvertVcfToPlink.out.bfile_per_chr_ch.collect()).set {bfile_ch}
    } else {
      Channel
        .from(params.bfile)
        .ifEmpty { exit 1, "Input plink prefix not found!" }
        .map { genotypes -> [file("${genotypes}.bed"), file("${genotypes}.bim"), file("${genotypes}.fam")]}
        .set { bfile_ch }
    }

    /*
     * Run interaction analysis
     */
    // if SNP-gene pairs (qtls_to_test) are provided use them in the interaction analysis, otherwise test all SNPs around the gene
    if (params.qtls_to_test == ''){
      features_to_test_ch = Channel.fromPath(params.genes_to_test)
    } else {
      features_to_test_ch = Channel.fromPath(params.qtls_to_test)
    }
    RUN_INTERACTION_QTL_MAPPING(norm_exp_ch, bfile_ch, covariates_ch, annotation_ch, Channel.of(params.covariate_to_test), chunk_ch.map { it[1] }, features_to_test_ch)
          
    /*
     * Stratified analysis if required
     */
    if (params.run_stratified){
      RUN_STRATIFIED_ANALYSIS(norm_exp_ch, bfile_ch, covariates_ch, annotation_ch, gte_ch, chunk_ch)
    }  
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
