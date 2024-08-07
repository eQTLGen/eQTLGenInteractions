#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

Channel.from(21..22)
  .map { chr -> tuple("$chr", file("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/test_nextflow/test_data/output/chr${chr}*bgen")) }
  .ifEmpty { exit 1, "Input .vcf.gz files not found!" } 
  .set { chr_vcf_pairs }
  

Channel
    .fromPath("$projectDir/data/ChunkingFile_test2.txt")
    .splitCsv( header: false )
    .map { row -> tuple(row[0].split(':')[0], row[0]) }
    .set { chunk_ch }
    

chunk_ch.join(chr_vcf_pairs).view()
    
    
