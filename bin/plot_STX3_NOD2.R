library(dplyr)
library(tidyverse)
library(optparse)
library(patchwork)

option_list <- list(
  make_option(c("-e", "--expression"), type = "character",
              help = ""),
  make_option(c("-b", "--pedfile"), type = "character",
              help = ""),
  make_option(c("-c", "--covariates"), type = "character",
              help = ""),
  make_option(c("-p", "--pcs"), type = "character", default = NULL,
              help = "")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

write(paste(names(args), ": ", unlist(args)), stderr())


nod2 <- "ENSG00000167207"
stx3 <- "ENSG00000166900"
snp = "rs1981760"

#args <- list(expression = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/work/d6/7ed5547de410d23cab97894c3770c7/exp_data_preprocessed.txt",
#    covariates = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/work/d6/7ed5547de410d23cab97894c3770c7/covariates.combined.txt",
#    pedfile = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/work/d6/7ed5547de410d23cab97894c3770c7/snp_geno.ped",
#    pcs = "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/work/d6/7ed5547de410d23cab97894c3770c7/exp_PCs.txt")

read_data <- function(INT = F){ 
  exp_path <- args$expression
  cov_path <- args$covariates
  pcs_path <- args$pcs
  ped_path <- args$pedfile
  
  
  # Read in the expression and covariate data
  expr <- as.data.frame(t(data.table::fread(exp_path, data.table = F)))
  colnames(expr) <- expr[1,]
  expr <- mutate_all(expr[-1,c(nod2,stx3)], function(x) as.numeric(as.character(x)))
  expr2 <- expr[,nod2, drop = FALSE]
  stx3_expr <- expr[,stx3, drop = FALSE]
  
  # PC correction
  if (!is.null(args$pcs) ){
    cat("Running PC correction\n")
    pcs <- read.delim(pcs_path, as.is = T, check.names = F, sep = "\t", row.names = 1)
    
    ids <- intersect(row.names(expr2), row.names(pcs))
    expr2 <- expr2[ids,]
    pcs <- pcs[ids,1:25] 
    
    tmp_data <- cbind(expr2, pcs)
    
    colnames(tmp_data)[1] <- "gene"
    lm_fit <- lm(gene ~ . ,data = tmp_data)
    expr3 <- as.data.frame(lm_fit$residuals)
    colnames(expr3) <- "NOD2"
  } else {
    expr3 <- expr2
    colnames(expr3) <- "NOD2"
  }
  
  # INT expression
  if (INT) {
    cat("Running INT on expression\n")
    expr_INT <- apply(expr3, 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
  } else {
    expr_INT <- expr3
  }

  # Read covariate data
  covar <- read.delim(cov_path, check.names = F, header = T, row.names = 1)  
  stx3_expr  <- apply(stx3_expr, 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
  tmp_ids <- intersect(row.names(covar), row.names(stx3_expr))
  covar2 <- cbind(covar[tmp_ids, c("AvgExprCorrelation")], stx3_expr[tmp_ids,])
  colnames(covar2) <- c("AvgExprCorrelation", "STX3")
  
  
  # Read genotype data
  geno <- read.table(ped_path)
  geno$genotypes <- (geno$V7 + geno$V8) - 2
  geno[geno$genotypes < 0, "genotypes"] <- NA
  geno <- geno[,c("V2", "genotypes")]
  row.names(geno) <- geno$V2
  
  
  # get ids for samples shared between expr_noadjession,  genotypes and sex
  ids <- intersect(intersect(row.names(expr_INT), row.names(geno)), row.names(covar2))
  

  expr <- expr[ids,]
  expr_INT <- expr_INT[ids,, drop = F]
  covar2 <- covar2[ids,]
  geno <- geno[ids,]
  
  return(list(expression = expr_INT, expression_raw = expr, genotypes = geno, covariates = covar2))
}

#### no INT ####
data_list <- read_data(INT = F)
m <- as.data.frame(cbind(data_list$expression[, "NOD2"], data_list$genotypes[,"genotypes"], data_list$covariates))
colnames(m) <- c("gene", "genotype", colnames(data_list$covariates))

lm_fit0 <- lm(gene ~ AvgExprCorrelation , data = m)
m$gene_resid <- residuals(lm_fit0)

lm_fit <- lm(gene ~ genotype * STX3 + AvgExprCorrelation , data = m)
pval <- formatC(summary(lm_fit)$coefficients["genotype:STX3",4], digits = 3)

lm_fit_snp <- lm(gene ~ genotype , data = m)
pval_snp <- formatC(summary(lm_fit_snp)$coefficients["genotype",4], digits = 3)

p1 <- ggplot(data = m , aes(y = gene_resid, color = as.factor(genotype), x = STX3))  +
    geom_point() +
    geom_smooth(aes(x = STX3, y = gene_resid, color = as.factor(genotype)), method = 'lm', formula = y ~ x) +
    ylab("NOD2") + xlab("STX3") +
    theme_bw() +
    theme(plot.title = element_text(size=12)) +
    labs(color=snp) + ggtitle(paste0("\ninteraction P = ", pval, "\n SNP P = ", pval_snp))

#### with INT ####
data_list <- read_data(INT = T)
m <- as.data.frame(cbind(data_list$expression[, "NOD2"], data_list$genotypes[,"genotypes"], data_list$covariates))
colnames(m) <- c("gene", "genotype", colnames(data_list$covariates))

lm_fit0 <- lm(gene ~ AvgExprCorrelation , data = m)
m$gene_resid <- residuals(lm_fit0)

lm_fit <- lm(gene ~ genotype * STX3 + AvgExprCorrelation , data = m)
pval <- formatC(summary(lm_fit)$coefficients["genotype:STX3",4], digits = 3)

lm_fit_snp <- lm(gene ~ genotype , data = m)
pval_snp <- formatC(summary(lm_fit_snp)$coefficients["genotype",4], digits = 3)

p2 <- ggplot(data = m , aes(y = gene_resid, color = as.factor(genotype), x = STX3))  +
    geom_point() +
    geom_smooth(aes(x = STX3, y = gene_resid, color = as.factor(genotype)), method = 'lm', formula = y ~ x) +
    ylab("NOD2 (inverse-rank transformed)") + xlab("STX3") +
    theme_bw() +
    theme(plot.title = element_text(size=12)) +
    labs(color=snp) + ggtitle(paste0("\ninteraction P = ", pval, "\n SNP P = ", pval_snp))
  
  
if (!is.null(args$pcs) ) {
    plot_path = "NOD2_STX3.adj_PCs.pdf"
} else {
    plot_path = "NOD2_STX3.pdf"
}

pdf(plot_path, width = 12, height = 6)
p1 + p2
dev.off()