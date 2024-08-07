library(optparse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(caret)

option_list <- list(
  make_option(c("-s", "--general_covs"), type = "character",
              help = "Path to the covariate file with sex and age. Sample ids = genotype ids."),
  make_option(c("-c", "--cell_counts"), type = "character", default = NULL,
              help = "Path to the covariate file with (deconvoluted) cell counts. Sample ids - expression ids."),
  make_option(c("-g", "--genotype_pcs"), type = "character", default = NULL,
              help = "Path to the covariate file with genotype PCs. Sample ids = genotype ids."),
  make_option(c("-r", "--rna_qual"), type = "character", default = NULL,
              help = "Path to the file with RNA quality (for RNA-seq). Sample ids = expression ids."),
  make_option(c("-i", "--gte"), type = "character",
              help = "Path to the genotype to expression id conversion file."),
  make_option(c("-o", "--out"), type = "character",
              help = "Output file name."),
  make_option(c("-n", "--rename"), type = "boolean",
              help = "Are raw expression sample ids different from genotype sample ids?.")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)


rename_samples <- function(d, gte, key_col = 1, value_col = 2){
  if (length(intersect(row.names(d), gte[,key_col])) < length(intersect(row.names(d), gte[,value_col]))){
    cat("Sample ids are already converted to expression ids! Skipping id conversion.\n")
    return(d)
  }
  d_m <- d[match(gte[,key_col], row.names(d), nomatch = 0),, drop = F]
  gte_m <- gte[match(row.names(d_m), gte[,key_col], nomatch = 0),]
  row.names(d_m) <- gte_m[,value_col]
  cat("Coverted genotype sample ids to expression sample ids.\nN samples before: ", nrow(d), "\nN samples after: ", nrow(d_m), "\n")
  return(d_m)
}


run_INT <- function(d){
  contin_covars <- apply(d, 2, function(x) length(unique(x)) > 3)
  d_int <- cbind(d[,!contin_covars], apply(d[,contin_covars], 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) ))
  colnames(d_int) <- c(colnames(d)[!contin_covars], colnames(d)[contin_covars])
  return(d_int)
}

# Combine general covariates with cell counts if they are provided

covar_main <- as.data.frame(read.delim(args$general_covs, check.names = F, header = T, row.names = 1, as.is = T, sep ="\t"))
gte <- read.delim(args$gte, check.names = F, header = F, as.is = T, sep ="\t")
#covar_main <- rename_samples(covar_main, gte)

if (! is.null(args$cell_counts) && args$cell_counts != "NA"){
  covar_cell_counts <- read.delim(args$cell_counts, check.names = F, header = T, row.names = 1, as.is = T, sep ="\t")
  covar_cell_counts <- rename_samples(covar_cell_counts, gte, key_col = 2, value_col = 1)
  covar <- merge(covar_main, covar_cell_counts, by = 0)
  row.names(covar) <- covar$Row.names
  covar$Row.names = NULL
  cat("Combined major covariates with cell counts\n")

} else {
  covar <- covar_main
  cat("No cell counts provided!\n")

}


covar <- covar[,colSums(is.na(covar)) < nrow(covar), drop = F]
covar <- na.omit(covar)

# Remove covariates with near zero variance
near_zero_var <- nearZeroVar(covar,  uniqueCut = 5)
if (length(near_zero_var) > 0) covar <- covar[, -near_zero_var]


#contin_covars <-  apply(covar, 2, function(x) length(unique(x)) > 3)
covar_names_for_INT <- names(which(apply(covar, 2, function(x) length(unique(x)) > 3)))
#covar_int <- cbind(covar[,!contin_covars], apply(covar[,contin_covars], 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) ))
#colnames(covar_int) <- c(colnames(covar)[!contin_covars], colnames(covar)[contin_covars])

cat("Combined major covariates with cell counts\n")
cat ("Number of covariates in the combined file:", ncol(covar), "\n")
cat("Number of samples with covariate data available:", nrow(covar), "\n")

# add genotype PCs
geno_pcs <- read.delim(args$genotype_pcs, check.names = F, header = T, row.names = 1, as.is = T)
geno_pcs <- geno_pcs[,1:4]
#geno_pcs <- rename_samples(geno_pcs[,1:4], gte)

covar_merged <- merge(covar, geno_pcs, by = 0)
row.names(covar_merged) <- covar_merged$Row.names
covar_merged$Row.names = NULL

cat("Added genotype PCs\n")
cat ("Number of covariates in the combined file:", ncol(covar_merged), "\n")
cat("Number of samples with covariate data available:", nrow(covar_merged), "\n")


if (! is.null(args$rna_qual)){
  rna_qual <- read.delim(args$rna_qual, check.names = F, header = T, row.names = 1, as.is = T)
  #rna_qual <- rename_samples(rna_qual, gte)
  covar_names_for_INT <- c(covar_names_for_INT, colnames(rna_qual)) 
  covar_merged <- merge(covar_merged, rna_qual, by = 0)
  row.names(covar_merged) <- covar_merged$Row.names
  covar_merged$Row.names = NULL
  
  cat("Added RNA quality\n")
  cat ("Number of covariates in the combined file:", ncol(covar_merged), "\n")
  cat("Number of samples with covariate data available:", nrow(covar_merged), "\n")
  
}

# run INT on general covariates, cell counts and RNA quality
covar_merged_int <- cbind(covar_merged[, !colnames(covar_merged) %in% covar_names_for_INT], apply(covar_merged[,covar_names_for_INT, drop = F], 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) ))
colnames(covar_merged_int) <- c(colnames(covar_merged[, !colnames(covar_merged) %in% covar_names_for_INT]), covar_names_for_INT)

# Plot covariate distributions prior to INT
pdf(paste0(args$out, ".distributions.pdf"), height = 5, width = 10)
#covar_long <- pivot_longer(covar_merged, cols = everything())
for (covar_name in colnames(covar_merged)){
   p1 <- ggplot(covar_merged, aes(get(covar_name))) + geom_histogram(bins = 50, color = "black", alpha = 0.5, fill = "dodgerblue4") + theme_bw() + xlab(covar_name) + ggtitle(paste0("Raw ", covar_name))
   p2 <- ggplot(covar_merged_int, aes(get(covar_name))) + geom_histogram(bins = 50, color = "black", alpha = 0.5, fill = "dodgerblue4") + theme_bw() + xlab(covar_name) + ggtitle(paste0("INT ", covar_name))
   print(p1 + p2)
}
dev.off()

# Write INT covariates 
write.table(covar_merged_int, file = args$out, sep = "\t", quote = F, col.names = NA)

