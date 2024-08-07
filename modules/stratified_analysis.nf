#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Split the combined covariate table into 2 according to the covariate of interest (E.g. sex or age). 
 * In case of a binary covariate the 2 values will be used for the split (E.g. males and females separately), 
 * in case of a quantitative covariate, samples falling in the first and last quartile will be written to the output files.
 */
process StratifyData {
    label "medium2"
    publishDir params.outdir, mode: 'copy'
    input:
        path (combined_covariates)
        val (covariate_name)
        path (norm_expression)
    output:
        path("*covariates*.txt"), emit: strat_covariates_ch
        path("*expression_summary*txt"), emit: strat_expr_summary_ch
    script:
    """
        Rscript $projectDir/bin/split_covariate_into_bins.R $combined_covariates $covariate_name $norm_expression ./
    """
}

/*
 * Run the eQTL mapping for the selected covariate level, testing all SNPs within a window from the genes to test
 */
process RunEqtlMappingPerGenePlink{
    label "long"
    tag "Chunk: $chunk"
    echo true

    input:
    tuple path(norm_expression), path(bed), path(bim), path(fam), path (covariates), path(limix_annotation), path(gte), val(covariate_to_test), val(chunk)


    output:
    //path "limix_out*/*"

    shell:
    '''
    geno=!{bed}
    plink_base=${geno%.bed}
    covar=!{covariates}
    tmp=${covar##*covariates.}
    outdir=${PWD}/limix_out_${tmp%.txt}/
    echo $outdir
    mkdir $outdir

    awk 'BEGIN {OFS="\\t"}; {print $2, $2}' !{fam} > gte.txt

    genes_to_test = !{params.genes_to_test}
    eval "$(conda shell.bash hook)"
    conda activate py39

    python /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_limix/limix_qtl/Limix_QTL/run_QTL_analysis.py \
     --plink ${plink_base} \
      -af !{limix_annotation} \
      -cf !{covariates} \
      -pf !{norm_expression} \
      -ff !{genes_to_test} \
      -od ${outdir} \
      -np 0 \
      -maf 0.05 \
      -c -gm gaussnorm \
      -w 1000000 \
      -hwe 0.0001 \
      -gr !{chunk} \
      -smf gte.txt
      
      #ls ${outdir}

      if [ ! -d !{params.outdir}/limix_out_${tmp%.txt}/ ]
      then
        mkdir !{params.outdir}/limix_out_${tmp%.txt}/
      fi
      if [ ! -z "$(ls -A ${outdir}/)" ]
      then
        cp ${outdir}/* !{params.outdir}/limix_out_${tmp%.txt}/
      else
        echo "No limix output to copy"
      fi
      #ls -la ${outdir}/

    '''
}

workflow RUN_STRATIFIED_ANALYSIS {
    take:
        norm_expression
        plink_geno
        covariates
        limix_annotation
        gte
        chunk

    main:
        covar_to_test_ch = Channel.of(params.covariate_to_test)
        StratifyData(covariates, covar_to_test_ch, norm_expression)
        RunEqtlMappingPerGenePlink(norm_expression.combine(plink_geno).combine(StratifyData.out.strat_covariates_ch.flatten()).combine(limix_annotation).combine(gte).combine(covar_to_test_ch).combine(chunk.map { it[1] }))
}

