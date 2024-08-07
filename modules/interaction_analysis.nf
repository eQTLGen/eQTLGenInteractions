#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * run QTL mapping per SNP-Gene pair
 */
process IeQTLmappingPerSNPGene {
    tag "Chunk: $chunk"

    input:
    tuple path(tmm_expression), path(bed), path(bim), path(fam), path(covariates), path(limix_annotation), path(gte), path(qtls_to_test), val(covariate_to_test), val(chunk)
    

    output:
    //path "limix_out/*"

    shell:
    '''
     geno=!{bed}
     plink_base=${geno%.bed}
     outdir=${PWD}/limix_out/
     mkdir $outdir
     ls $PWD

     python /limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \
     --plink ${plink_base} \
      -af !{limix_annotation} \
      -cf !{covariates} \
      -pf !{tmm_expression} \
      -smf !{gte} \
      -fvf !{qtls_to_test} \
      -od ${outdir} \
      --interaction_term !{covariate_to_test} \
      -gr !{chunk} \
      -np 0 \
      -maf 0.05 \
      -c -gm gaussnorm \
      -w 1000000 \
      -hwe 0.0001
      
      echo !{params.outdir}
      if [ ! -d !{params.outdir}/limix_output/ ]
      then
        echo "test folder create"

      	mkdir !{params.outdir}/limix_output/
      fi
      if [ ! -z "$(ls -A ${outdir}/)" ]
      then
      	cp ${outdir}/* !{params.outdir}/limix_output/
      else
	echo "No limix output to copy"      
      fi
      ls -la ${outdir}/  
      
    '''
}


process IeQTLmappingPerGeneTMP {
    label "long"

    tag "Chunk: $chunk"
    echo true

    input:
    tuple path(tmm_expression), path(bed), path(bim), path(fam), path(covariates), path(limix_annotation), path(gte), path(genes_to_test), val(covariate_to_test), val(chunk)


    output:

    shell:
    '''
    echo "Running chunk !{chunk}"
    '''
}

/*
 * run QTL mapping for all SNPs arounnd the specified Gene for a given chunk
 */
process IeQTLmappingPerGene {
    label "verylong"
    errorStrategy 'ignore'
    tag "Chunk: $chunk"
    //echo true

    input:
    tuple path(tmm_expression), path(bed), path(bim), path(fam), path(covariates), path(limix_annotation), path(gte), path(genes_to_test), val(covariate_to_test), val(chunk)


    output:
    //path "limix_out/*"

    shell:
    '''
    geno=!{bed}
    plink_base=${geno%.bed}
    outdir=${PWD}/limix_out/
    mkdir $outdir

    awk 'BEGIN {OFS="\\t"}; {print $2, $2}' !{fam} > gte.txt

    eval "$(conda shell.bash hook)"
    conda activate py39
    #python /limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py

    python /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_limix/limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \
    --plink ${plink_base} \
    -af !{limix_annotation} \
    -cf !{covariates} \
    -pf !{tmm_expression} \
    -ff !{genes_to_test} \
    -od ${outdir} \
    --interaction_term !{covariate_to_test} \
    -np 20 \
    -maf 0.05 \
    -c -gm gaussnorm \
    -w 1000000 \
    -hwe 0.0001 \
    -gr !{chunk} \
    -smf gte.txt \
    --write_permutations --write_zscore
    
    
    if [ ! -d !{params.outdir}/limix_output/ ]
    then
            echo "test folder create2"
      mkdir !{params.outdir}/limix_output/
    fi
    if [ ! -z "$(ls -A ${outdir}/)" ]
    then
      cp ${outdir}/* !{params.outdir}/limix_output/
    else
      echo "No limix output to copy"
    fi
    ls -la ${outdir}/

    '''
}

process IeQTLmappingPerGeneNoChunks {
    label "long"
    //echo true
    input:
    tuple path(tmm_expression), path(bed), path(bim), path(fam), path(covariates), path(limix_annotation), path(gte), path(genes_to_test), val(covariate_to_test)
    

    output:
    //path "limix_out/*"
     
    shell:
    '''
     geno=!{bed}
     plink_base=${geno%.bed}
     outdir=${PWD}/limix_out/
     mkdir $outdir
     
     python /limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \
     --plink ${plink_base} \
      -af !{limix_annotation} \
      -cf !{covariates} \
      -pf !{tmm_expression} \
      -ff !{genes_to_test} \
      -od ${outdir} \
      --interaction_term !{covariate_to_test} \
      -np 0 \
      -maf 0.05 \
      -c -gm gaussnorm \
      -w 1000000 \
      -hwe 0.0001 \
      -gr 2:1-400000000
      
      echo !{params.outdir}
      if [ ! -d !{params.outdir}/limix_output/ ]
      then
        mkdir !{params.outdir}/limix_output/
      fi
      if [ ! -z "$(ls -A ${outdir}/)" ]
      then
              echo "test folder create3"
        cp ${outdir}/* !{params.outdir}/limix_output/
      else
        echo "No limix output to copy"
      fi
      ls -la ${outdir}/

    '''
}


/*
 * run QTL mapping for all SNPs arounnd the specified Gene for a given chunk
 */
process IeQTLmappingPerGeneBgen {
    label "long"
    tag "Chunk: $chunk"
    //echo true

    input:
    tuple path(tmm_expression), path(covariates), path(limix_annotation), path(gte), path(genes_to_test), val(covariate_to_test), val(chr), val(chunk), path(bgen), path(bgen_sample)


    output:
    //path "limix_out/*"

    shell:
    '''
     geno=!{bgen}
     bgen_base=${geno%.bgen}
     outdir=${PWD}/limix_out/
     mkdir $outdir

     python /limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py \
     --bgen ${bgen_base} \
      -af !{limix_annotation} \
      -cf !{covariates} \
      -pf !{tmm_expression} \
      -ff !{genes_to_test} \
      -od ${outdir} \
      --interaction_term !{covariate_to_test} \
      -np 0 \
      -maf 0.05 \
      -c -gm gaussnorm \
      -w 1000000 \
      -hwe 0.0001 \
      -gr !{chunk} \
      -smf /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/output/tmp_gte.txt
      
      echo !{params.outdir}
      if [ ! -d !{params.outdir}/limix_output/ ]
      then
              echo "test folder create4"
        mkdir !{params.outdir}/limix_output/
      fi
      if [ ! -z "$(ls -A ${outdir}/)" ]
      then
        cp ${outdir}/* !{params.outdir}/limix_output/
      else
        echo "No limix output to copy"
      fi
      ls -la ${outdir}/

    '''
}

process FilterGenesToTest {
    label 'small'
    echo true
    publishDir params.outdir, mode: 'copy'

    input:
        tuple path (genes), path (normalized_expression_data)
    
    output:
        path "genes_to_test.txt"
    
    """
#!/usr/bin/env python
import pandas as pd

genes_to_test = pd.read_csv(\"$genes\", sep = "\t")
genes_expressed = pd.read_csv(\"$normalized_expression_data\", sep = "\t", usecols=[0])

gene_overlap = set(genes_to_test.iloc[:,0]).intersection(set(genes_expressed.iloc[:,0]))

out_f = open('genes_to_test.txt', 'w')
for g in gene_overlap:
  _ = out_f.write(g + '\\n')

"""

}
