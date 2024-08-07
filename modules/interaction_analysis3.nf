#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * run QTL mapping per SNP-Gene pair or for all SNPs around the gene depending on command line parameters.
 * Covariates will be included in the model as linear terms, and the $covariate_to_test will also be included as an interaction with genotype
 */
process IeQTLmapping {
    tag "Chunk: $chunk"

    //publishDir "${params.outdir}", mode: 'copy', overwrite: true, failOnError: true


    input:
        tuple path(tmm_expression), path(bed), path(bim), path(fam), path(covariates), path(limix_annotation), val(covariate_to_test), val(chunk), path(qtl_ch)
    
    // make the output optional for the case when there are no eQTLs to test and the output is empty. If it's not optional then .collect() in the workflow description will not work
    output:
        path "limix_out/*", optional: true

    shell:
    '''

    geno=!{bed}
    plink_base=${geno%.bed}
    outdir=${PWD}/limix_out/
    mkdir -p $outdir

    # make a fake gte because limix doesn't work without it
    awk 'BEGIN {OFS="\\t"}; {print $2, $2}' !{fam} > gte.txt

    # determine whether to test SNP-gene pairs or all SNPs within 1Mb
    qtls=!{params.qtls_to_test}
    genes=!{params.genes_to_test}
    if [ "${#qtls}" -gt 1 ]
    then 
        arg_line="-fvf !{qtl_ch}"
    else
        arg_line="-ff !{qtl_ch} -w 1000000 "
    fi

    # --interaction_term is a colum from -cf
    # all variables from -cf are used in model without interaction
    # -gm gaussnorm does an inverse normal transformation of the expression data
    
    python /limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \
     --plink ${plink_base} \
      -af !{limix_annotation} \
      -cf !{covariates} \
      -pf !{tmm_expression} \
      -smf gte.txt \
      ${arg_line} \
      -od ${outdir} \
      --interaction_term !{covariate_to_test} \
      -gr !{chunk} \
      -np !{params.num_perm} \
      -maf 0.01 \
      -c \
      -gm gaussnorm \
      -hwe 0.0001 \
      --write_permutations \
      --write_zscore


      
    '''
}



/*
 * Plot the interaction of rs1981760-NOD2 cis-eQTL with STX3 gene expression (a proxy for neutrophil percentage) - just a sanity check
 */
process PlotSTX3NOD2 {
    label "short"
    publishDir params.outdir, mode: 'copy', overwrite: true, failOnError: true

    input:
        tuple path (expression_path), path (covariate_path), path (bed), path (bim), path (fam), path (exp_PCs)

    output:
        path ("*.pdf")

    shell:
    if (params.expr_pcs == ''){   
        '''
        geno=!{bed}
        plink_base=${geno%.bed}
        echo !{bed}
            !{projectDir}/tools/plink --bfile $plink_base --snp rs1981760 --recode 12 --out snp_geno
        Rscript !{projectDir}/bin/plot_STX3_NOD2.R -e !{expression_path} -c !{covariate_path} -b snp_geno.ped
        '''
    } else {
        '''
        geno=!{bed}
        plink_base=${geno%.bed}
        !{projectDir}/tools/plink --bfile $plink_base --snp rs1981760 --recode 12 --out snp_geno
        Rscript !{projectDir}/bin/plot_STX3_NOD2.R -e !{expression_path} -c !{covariate_path} -b snp_geno.ped -p !{exp_PCs}
        '''
    }
}

/*
 * Convert the interaction analysis output folder to a text file
 */
process ConvertIeQTLsToText {
    echo true

     publishDir params.outdir, mode: 'copy', overwrite: true, failOnError: true

    input:
    path limix_out_files
    
    output:
        path "${params.covariate_to_test}.iqts.txt.gz*"
        //path "limix_out_text/*txt.gz"
        //path "limix_out_text/limix_out_text.md5"


    script:
    """

    python /limix_qtl/Limix_QTL/post_processing/minimal_interaction_postprocess.py \
      -id ./ \
      -od  ./ \
      -sfo \
      --write_compressed

    mv iqtl_results_all.txt.gz ${params.covariate_to_test}.iqts.txt.gz
    md5sum  ${params.covariate_to_test}.iqts.txt.gz > ${params.covariate_to_test}.iqts.txt.gz.md5


    """
}



workflow RUN_INTERACTION_QTL_MAPPING {   
    take:
        tmm_expression
	    plink_geno
        covariates_ch
	    limix_annotation
        covariate_to_test
	    chunk
        qtl_ch

    main:

    expr_pcs_ch = params.expr_pcs
            ? Channel.fromPath(params.expr_pcs, checkIfExists:true)
            : Channel.fromPath('EMPTY')
    // if run interaction analysis with covariate * genotype interaction terms, preadjust gene expression for other, linear covariates

        interaction_ch = tmm_expression.combine(plink_geno).combine(covariates_ch).combine(limix_annotation).combine(covariate_to_test).combine(chunk).combine(qtl_ch)
        ConvertIeQTLsToText(IeQTLmapping(interaction_ch).collect())

    // Plot the interaction of NOD2 cis-eQTL with STX3 (neutrophil proxy) as a sanity check
    PlotSTX3NOD2(tmm_expression.combine(covariates_ch).combine(plink_geno).combine(expr_pcs_ch))
}
