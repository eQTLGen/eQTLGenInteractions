library(optparse)
library(valr)
library(snpStats)
library(tidyverse, lib.loc = "/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/R/x86_64-pc-linux-gnu-library/4.0.3/")

#option_list <- list(
#  make_option(c("-e", "--expression"), type = "character",
#              help = "Normalized expression data without INT"),
#  make_option(c("-b", "--bfile"), type = "character",
#              help = "Plink genotype data prefix"),
#  make_option(c("-c", "--covariates"), type = "character",
#              help = "Covariate data. Should include RNA quality ,column name: AvgExprCorrelation",
#  make_option(c("-r", "--chunk"), type = "character",
#              help = "Genomic chunk chrom:start-end"),
#  make_option(c("-a", "--annot"), type = "character",
#              help = "Gene annotation file in Limix_QTL format"),
#  make_option(c("-l", "--eqtls"), type = "character",
#              help = "eQTLs to test. Columns should be called feature_id and snp_id"),
#  make_option(c("-g", "--egenes"), type = "character",
#              help = "eGenes to test. All SNPs within 1Mb will be tested."),
#  make_option(c("-n", "--cov_name"), type = "character",
#              help = "Covairate to test (E.g. sex)"),
#  make_option(c("-p", "--pcs"), type = "character", default = NULL,
#              help = "Expression PCs to regress out from expression data prior to analysis."),
#  make_option(c("-m", "--model"), type = "numeric",
#              help = "model=0: gene ~ ., model=1: gene ~ . + covariate * genotype; model=2: gene ~ . + covariate * genotype + cell_perc * genotype; model=3: gene * eQTL interaction analysis"),
#  make_option(c("-x", "--chr"), type = "character",
#              help = "Chromosome"),
#  make_option(c("-o", "--output"), type = "character",
#              help = "Output folder where to put the filtered data matrix. The interaction analysis results are printed to stdout.")
#)

option_list <- list(
  make_option(c("-e", "--expression"), type = "character",
              help = "Normalized expression data without INT"),
  make_option(c("-b", "--bfile"), type = "character",
              help = "Plink genotype data prefix"),
  make_option(c("-c", "--covariates"), type = "character",
              help = "Covariate data. Should include RNA quality in column AvgExprCorrelation"),
  make_option(c("-r", "--chunk"), type = "character",
              help = "Genomic chunk chrom:start-end"),
  make_option(c("-a", "--annot"), type = "character",
              help = "Gene annotation file in Limix_QTL format"),
  make_option(c("-l", "--eqtls"), type = "character",
              help = "eQTLs to test. Columns should be called feature_id and snp_id"),
  make_option(c("-g", "--egenes"), type = "character",
              help = "eGenes to test. All SNPs within 1Mb will be tested."),
  make_option(c("-n", "--cov_name"), type = "character",
              help = "Covairate to test (E.g. sex)"),
  make_option(c("-p", "--pcs"), type = "character", default = NULL,
              help = "Expression PCs (or other covariates) to regress out from expression data prior to analysis."),
  make_option(c("-m", "--model"), type = "numeric",
              help = "model=0: gene ~ ., model=1: gene ~ . + covariate * genotype; model=2: gene ~ . + covariate * genotype + cell_perc * genotype; model=3: gene * eQTL interaction analysis"),
  make_option(c("-x", "--chr"), type = "character",
              help = "Chromosome"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output folder where to put the filtered data matrix. The interaction analysis results are printed to stdout.")
)
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

write(paste(names(args), ": ", unlist(args)), stderr())


# Create a list of eQTLs to test in the specified chunk.
create_eqtl_list <- function() {
  egenes_path <- args$egenes
  eqtlts_path <- args$eqtls
  gene_path = args$annot
  plink_path = args$bfile
  chunk = args$chunk
  genome = as_tibble(matrix(c(1,248956422,2,242193529,3,198295559,4,190214555,5,181538259,6,170805979,7,159345973,8,145138636,9,138394717,10,133797422,11,135086622,12,133275309,13,114364328,14,107043718,15,101991189,16,90338345,17,83257441,18,80373285,19,58617616,20,64444167,21,46709983,22,50818468,23,156040895,24,57227415), byrow = T, ncol = 2))
  colnames(genome) <- c("chrom", "size")
  
  bim <- read.delim(paste0(plink_path, ".bim"), check.names = F, header = F)
  genotyped_snps <- tibble(chrom = as.numeric(bim[,1]), start = as.numeric(bim[,4]) - 1, end = as.numeric(bim[,4]), id = bim[,2])
  
  annot <- read.delim(gene_path, check.names = F, header = T, sep = "\t")
  annot <- annot[annot$chromosome == args$chr,]
  gene_coord <- tibble(chrom = as.numeric(annot$chromosome), start = as.numeric(annot$start) - 1, end = as.numeric(annot$end), id = annot$feature_id)
  
  chunk_split <- unlist(strsplit(chunk, '[\\:\\-]'))
  chunk_bed <- tibble(chrom = as.numeric(chunk_split[1]), start = as.numeric(chunk_split[2]) - 1, end = as.numeric(chunk_split[3]))
  
  gene_coord <- bed_intersect(gene_coord, chunk_bed)[,1:4]
  colnames(gene_coord) <- c("chrom", "start", "end", "id")
  
  if (is.null(args$eqtls)){ # get all SNPs within 1Mb from gene start 
    egenes <- read.delim(args$egenes, as.is = T, check.names = F)[,1]
    
    gene_coord <- gene_coord[gene_coord$id %in% egenes,]
    
    intersection <- bed_window(gene_coord, genotyped_snps, genome, both = 1e6)
    
    snp_genes <- intersection[,c("id.x", "id.y")]
    colnames(snp_genes) <- c("gene", "snp")
    
  } else { # select SNPs that are present in the genotype data and genes that are located in the given chunk
    eqtls <- read.delim(args$eqtls, as.is = T, check.names = F, sep = "\t")
    eqtls <- eqtls[eqtls$snp_id %in% bim[,2] & eqtls$feature_id %in% gene_coord$id,]
    
    snp_genes <- eqtls[,c("feature_id", "snp_id")]
    colnames(snp_genes) <- c("gene", "snp")
  }
  
  return(snp_genes)
}

# Extract numeric SNP genotypes for the sample ids
get_snp_genotypes <- function(geno, snp, ids) {
  snp_num <- which(geno$map$snp.name == snp)
  geno_num <- as(geno$genotypes[ids, snp], "numeric")
  
  snp_geno_map <- geno$map[snp_num,]
  ea <- snp_geno_map$allele.2
  oa <- snp_geno_map$allele.1
  chr <- snp_geno_map$chromosome
  pos <- snp_geno_map$position
  return(list(genotypes = geno_num, ea = ea, oa = oa, chr = chr, pos = pos))
}

# Read expression, covariate and genotype data
read_data_plink <- function(exp_path, cov_path, plink_path, snp_genes){
  # Read in the expression and covariate data
  expr <- as.data.frame(t(read.delim(exp_path, check.names = F, header = T, row.names = 1)))
  covar <- read.delim(cov_path, check.names = F, header = T, row.names = 1)
  covar <- covar %>% 
    drop_na() %>% 
    mutate_if(~n_distinct(.) < 3, as.factor)
  
  if (!is.null(args$pcs) ){
    write("\tAdjusting expression for 25 PCs!", stderr())
    expr <- regress_linear_covariates(expr, args$pcs)
  }
  
  # INT expression
  expr <- apply(expr, 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
  
  geno <- read.plink(plink_path, select.snps = unique(snp_genes$snp))
  
  # get ids for samples shared between expression,  genotypes and covariates
  ids <- intersect(intersect(row.names(expr), geno$fam$member), row.names(covar))
  expr <- expr[ids,]
  covar <- covar[ids,]
  
  return(list(expression = expr, genotypes = geno, covariates = covar))
}

# Regress out expression PCs (or other covariates) from expression data
regress_linear_covariates <- function(expr, covar_path, write_table = T){
  covars <- read.delim(covar_path, as.is = T, check.names = F, sep = "\t", row.names = 1)
  covars <- na.omit(covars)
  ids <- intersect(row.names(expr), row.names(covars))
  expr <- expr[ids,]
  covars <- covars[ids,]
  
  new_expr <- data.frame(matrix(ncol = ncol(expr), nrow = nrow(expr)))
  colnames(new_expr) <- colnames(expr)
  row.names(new_expr) <- row.names(expr)
  
  for (col in 1:ncol(expr)){
    tmp_data <- cbind(expr[,col], covars)
    colnames(tmp_data)[1] <- "gene"
    lm_fit <- lm(gene ~ . ,data = tmp_data)
    new_expr[,col] <- lm_fit$residuals
  }
  
  if (write_table & !is.null(args$output)) {
    write.table(new_expr, file = paste0(args$output, "/expression_corrected.noINT.txt"), sep = "\t", quote = F, col.names = NA)
  }
  return(new_expr)
}

################
## Analysis
################

# Create a list of SNP-gene pairs to test
snp_gene_pairs <- create_eqtl_list()

#read in data
data_list <- read_data_plink(args$expression, args$covariates, args$bfile, snp_gene_pairs)
sample_ids <- row.names(data_list$expression)
snp_map <- data_list$genotypes$map

# keep only SNPs and genes that are present in the data
snp_gene_pairs <- as.data.frame(snp_gene_pairs[snp_gene_pairs$gene %in% colnames(data_list$expression) & snp_gene_pairs$snp %in% data_list$genotypes$map$snp.name,])

write(paste0("Total number of SNP-gene pairs to test: ", nrow(snp_gene_pairs)), stderr())

# Determine the model formula
# No interactions:
if (args$model == 0) { 
  model_formula = "."
# Base model with interactions only for covariate of interest and RNA-quality 
} else if (args$model == 1){
  if ("AvgExprCorrelation" %in% colnames(data_list$covariates)) {
    model_formula = ". + covariate*genotype + AvgExprCorrelation*genotype"
  } else {
    model_formula = ". + covariate*genotype"
  }
# Model including genotype interactions with all covariates (except for genotype PCs)
} else if (args$model == 2){
  covars_for_interaction <- colnames(data_list$covariates)[!grepl("PC[1-9]", colnames(data_list$covariates))]
  covars_for_interaction <- gsub(paste0('\\b', args$cov_name, '\\b'), "covariate", covars_for_interaction)
  model_formula = paste(c(".", paste(covars_for_interaction, "*genotype", sep = "")), collapse = " + ")
# Gene x eQTL interaction analysis
} else if (args$model == 3){
  model_formula = "gene-gene"
} else {
  cat ("Wrong model formula!")
}

write(paste0("model formula:", model_formula), stderr())

if (model_formula != "gene-gene"){
  # Phenotype x eQTL interaction analysis
  for (e in 1:nrow(snp_gene_pairs)){
    gene <- snp_gene_pairs[e, "gene"]
    snp <- snp_gene_pairs[e, "snp"]
    
    snp_genotypes <- get_snp_genotypes(data_list$genotypes, snp, sample_ids)
    
    m <- as.data.frame(cbind(data_list$expression[, gene], snp_genotypes[[1]], data_list$covariates))
    colnames(m) <- c("gene", "genotype", colnames(data_list$covariates))
    colnames(m) <- gsub(paste0('\\b', args$cov_name, '\\b'), "covariate", colnames(m))
    
    if (min(table(m$genotype)) > 10){ # skip cases when # alt homo is less than 10
      lm_fit <- lm(as.formula(paste0("gene ~ ", model_formula)), data = m)
      coef <- summary(lm_fit)$coefficients
      row.names(coef) <- gsub("covariate2", "covariate", row.names(coef))
      cat(gene, snp, snp_genotypes$ea, snp_genotypes$oa, snp_genotypes$chr, snp_genotypes$pos, coef["genotype:covariate",c(1,2,4)], "\n", sep = "\t")
    }
    if (e %% 100 == 0){
      write(paste0("processed", e, "eQTLs"), stderr())
    }
  }
  
  
} else {
  # Gene x eQTL interaction analysis
  for (e in 1:nrow(snp_gene_pairs)){
    
    gene <- snp_gene_pairs[e, "gene"]
    snp <- snp_gene_pairs[e, "snp"]
    
    snp_genotypes <- get_snp_genotypes(data_list$genotypes, snp, sample_ids)
    
    for (covar_gene in colnames(data_list$expression)){

      # if covariate gene is not the same as eGene and if covariate is not in the list of eQTLs with the same SNP      
      if (covar_gene != gene && nrow(snp_gene_pairs[snp_gene_pairs$gene == covar_gene & snp_gene_pairs$snp == snp,]) == 0){
        m <- as.data.frame(cbind(data_list$expression[, gene], data_list$expression[, covar_gene], snp_genotypes[[1]], data_list$covariates))
        colnames(m) <- c("gene", "covariate","genotype", colnames(data_list$covariates))
        
        if (min(table(m$genotype)) > 10){ # skip cases when # alt homo is less than 10
          lm_fit <- lm(gene ~ genotype*covariate + genotype*AvgExprCorrelation, data = m)
          coef <- summary(lm_fit)$coefficients
          cat(gene, snp, covar_gene,snp_genotypes$ea, snp_genotypes$oa, snp_genotypes$chr, snp_genotypes$pos, coef["genotype:covariate",c(1,2,4)], "\n", sep = "\t")
        }     
      }
    }
  }
}

