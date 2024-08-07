#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name="InteractionAnalysis"

# These are needed modules in UT HPC to get singularity and Nextflow running. Replace with appropriate ones for your HPC.
ml nextflow

# We set the following variables for nextflow to prevent writing to your home directory (and potentially filling it completely)
# Feel free to change these as you wish.
export SINGULARITY_CACHEDIR=../singularitycache
export NXF_HOME=../nextflowcache

# Disable pathname expansion. Nextflow handles pathname expansion by itself.
set -f

# Cohort name
c=NTR_GoNL
# Genotype data
#[full path to the folder with imputed filtered vcf files produced by eQTLGen pipeline 2_Imputation step (postimpute folder)]
vcf_dir_path=/groups/umcg-bios/tmp01/projects/BIOS_for_eQTLGenII/pipeline/20220426/2_Imputation/out/${c}/postimpute/

# raw expression data (same as input to DataQC step)
raw_exp_path=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/run1/data/${c}/${c}_raw_expression.txt.gz
# normalized expression data (output of the DataQC step)
norm_exp_path=/groups/umcg-bios/tmp01/projects/BIOS_for_eQTLGenII/pipeline/20220426/1_DataQC/out/${c}/outputfolder_exp/exp_data_QCd/exp_data_preprocessed.txt
# File that contains cohort covariates: E.g. sex and age. Sample ids should be the same as in the genotype data
covariate_path=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/run1/BIOS_covariates.txt
# genotype to expression coupling file
gte_path=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/run1/data/${c}/${c}.gte

# covariate to test (name of the sex column in the covariate file)
covariate_to_test=gender_F1M2

# Path to genotype PCs (output of dataQC step)
genotype_pcs_path=/groups/umcg-bios/tmp01/projects/BIOS_for_eQTLGenII/pipeline/20220426/1_DataQC/out/${c}/outputfolder_gen/gen_PCs/GenotypePCs.txt
# Path to expression PCs (output of dataQC step)
expression_pcs_path=/groups/umcg-bios/tmp01/projects/BIOS_for_eQTLGenII/pipeline/20220426/1_DataQC/out/${c}/outputfolder_exp/exp_PCs/exp_PCs.txt

# output folder (needs to exist)
output_path=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/run2/results/${c}_interactions/


# Path to the nextflow interaction analysis folder
script_folder=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/ieQTL_nextflow_pipeline/

qtls_to_test=${script_folder}/data/sign_qtls_cistrans.txt.gz
chunk_file=${script_folder}/data/ChunkingFile.GRCh38.110.txt
exp_platform=RNAseq


# Command:
NXF_VER=21.10.6 nextflow run ${script_folder}/InteractionAnalysis.nf \
--vcf_dir $vcf_dir_path \
--raw_expfile ${raw_exp_path} \
--norm_expfile ${norm_exp_path} \
--gte ${gte_path} \
--covariates $covariate_path \
--exp_platform ${exp_platform} \
--cohort_name ${c} \
--covariate_to_test $covariate_to_test \
--qtls_to_test $qtls_to_test \
--genotype_pcs $genotype_pcs_path \
--chunk_file $chunk_file \
--outdir ${output_path}  \
--run_stratified false \
--preadjust true \
--cell_perc_interactions true \
-resume \
-profile singularity,slurm 

