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
 * run QTL mapping per SNP-Gene pair or for all SNPs around the gene depending on command line parameters.
 * covariates will be included both as linear terms and as interactions with genotypes. 
 */
process IeQTLmapping_InteractionCovariates {
    tag "Chunk: $chunk"

    //publishDir "${params.outdir}", mode: 'copy', overwrite: true, failOnError: true

    input:
        tuple path(preadjusted_expression), path(bed), path(bim), path(fam), path(interaction_covariates), path(limix_annotation), val(covariate_to_test), val(chunk), path(qtl_ch)
    
    output:
        path "limix_out/*", optional: true

    shell:
    '''
    geno=!{bed}
    plink_base=${geno%.bed}
    outdir=${PWD}/limix_out/
    mkdir $outdir

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
      -cf !{interaction_covariates} \
      -pf !{preadjusted_expression} \
      -smf gte.txt \
      $arg_line \
      -od ${outdir} \
      --interaction_term !{covariate_to_test} \
      -gr !{chunk} \
      -np !{params.num_perm} \
      -maf 0.01 \
      -c \
      -gm gaussnorm \
      -hwe 0.0001 \
      --interaction_covariates \
      --write_permutations \
      --write_zscore
      

    '''
}


/*
 * This process is used when the user asks to add covariate * genotype interaction terms to the Limix interaction analysis model
 * The process splits the covariate table into 2 tables: linear covariates that should be regressed out before the interaction analysis and interaction covariates that will be included in the interaction model together with their interaction with genotypes
 */
process SplitCovariates {
    tag "Split covariates"
    label "short"
    //publishDir params.outdir, mode: 'copy'

    input:
        path covariates_path

    output:
        path ("covariates_preadjust.txt"), emit: linear_covariates_ch
        path ("covariates_interaction.txt"), emit: interaction_covariates_ch

    shell:
    '''
    # The covariate of interest (e.g. sex or age) should always be added as an interaction term with genotype. If the data is RNA-seq, then RNA quality should also be included as interaction term 
    if [ !{params.exp_platform} = "RNAseq" -o !{params.exp_platform} = "RNAseq_HGNC" ]; then
          interaction_covariates=`echo "AvgExprCorrelation,!{params.covariate_to_test}"`
    else 
          interaction_covariates = `echo "!{params.covariate_to_test}"`
    fi
    
    # if the user specified to add interactions of cell proportions with genotype to the model, then get the column names of cell counts file and include them in the list of interaction covariates
    # then split the combined covariate table into interaction and linear covariates
    if [ !{params.cell_perc_interactions} = false ]; then
        python !{projectDir}/bin/extract_columns_from_file.py !{covariates_path} ${interaction_covariates} covariates_interaction.txt
        python !{projectDir}/bin/extract_columns_from_file.py !{covariates_path} ${interaction_covariates} covariates_preadjust.txt -v
    else 
        cell_types=`cut -f2- !{params.outdir}/cell_counts.txt | awk 'BEGIN {FS="\t"; OFS=","}; {if (NR == 1) {$1=$1; print}}'`
        interaction_covariates2=`echo "$interaction_covariates,$cell_types"`

        python !{projectDir}/bin/extract_columns_from_file.py !{covariates_path} ${interaction_covariates2} covariates_interaction.txt
        python !{projectDir}/bin/extract_columns_from_file.py !{covariates_path} ${interaction_covariates2} covariates_preadjust.txt -v
    fi
    '''
}

/*
 * Regress out linear covariates from the expression table prior to the interaction analysis. Also regress out expression PCs if provided in the parameters.
 */
process PreadjustExpression {
    tag "Preadjust expression"
    label "short"
    //publishDir params.outdir, mode: 'copy'

    input:
    path expression_path
    path covariates_path
    path pcs_path

    output:
    path ("expression_corrected.noINT.txt")

    shell:
    if (params.expr_pcs == ''){
        '''
        Rscript !{projectDir}/bin/regress_linear_covariates.R !{expression_path} !{covariates_path} ./
        '''
    } else {
    '''
      n=!{params.num_expr_PCs}
      n2=$((n + 1))
      cut -f1-$n2 !{pcs_path} > pcs.txt
      Rscript !{projectDir}/bin/regress_linear_covariates.R !{expression_path} !{covariates_path} ./ pcs.txt
    
    '''
    }
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
    if (params.preadjust){   


        SplitCovariates(covariates_ch)
        PreadjustExpression(tmm_expression, SplitCovariates.out.linear_covariates_ch, expr_pcs_ch)

        interaction_ch = PreadjustExpression.out.combine(plink_geno).combine(SplitCovariates.out.interaction_covariates_ch).combine(limix_annotation).combine(covariate_to_test).combine(chunk).combine(qtl_ch)
        ConvertIeQTLsToText(IeQTLmapping_InteractionCovariates(interaction_ch).collect())
    // Run interaction analysis with covariates added to the model as linear term (except the covariate of interest, which is always included as an interaction term)
    } else {
        interaction_ch = tmm_expression.combine(plink_geno).combine(covariates_ch).combine(limix_annotation).combine(covariate_to_test).combine(chunk).combine(qtl_ch)
        ConvertIeQTLsToText(IeQTLmapping(interaction_ch).collect())
    }
    // Plot the interaction of NOD2 cis-eQTL with STX3 (neutrophil proxy) as a sanity check
    PlotSTX3NOD2(tmm_expression.combine(covariates_ch).combine(plink_geno).combine(expr_pcs_ch))
}
