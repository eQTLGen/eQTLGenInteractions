library(data.table)
#library(preprocessCore, lib.loc = '/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/R/x86_64-pc-linux-gnu-library/4.0.3/')
library(preprocessCore)
library(edgeR)
library(ggplot2)
library(optparse)
library(MASS)
library(dplyr)

setDTthreads(8)

# Argument parser
option_list <- list(
    make_option(c("-e", "--expression_matrix"), type = "character",
    help = "Unprocessed gene expression matrix from array or RNA-seq experiment. Samples in columns, genes/probes in the rows."),
    make_option(c("-n", "--norm_expression_matrix"), type = "character",
    help = "Normalized and filtered gene expression matrix, output from DataQC step. Samples in columns, genes/probes in the rows."),
    make_option(c("-l", "--genotype_to_expression_linking"), type = "character",
    help = "Genotype-to-expression linking file. No header. First column- sample IDs in the genotype file, second column- corresponding sample IDs in the gene expression matrix."),
    make_option(c("-p", "--platform"), type = "character",
    help = "Gene expression platform. This determines the normalization method and replaces probes with best-matching genes based on empirical probe mapping. One of: HT12v3, HT12v4, HuRef8, RNAseq, AffyU219, AffyHumanExon, RNAseq_HGNC."),
    make_option(c("-m", "--emp_probe_mapping"), type = "character",
    help = "Empirical probe matching file. Used to link the best array probe to each blood-expressed gene."),
    make_option(c("-o", "--output"), type = "character",
    help = "Output folder where to put the filtered data matrix.")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Debug
print(args$expression_matrix)
print(args$genotype_to_expression_linking)
print(args$norm_expression_matrix)
print(args$platform)
print(args$emp_probe_mapping)
print(args$output)

# Make output folder structure
dir.create(args$output)

#############
# functions #
#############

Z_transform <- function(x){
    z <- (x - mean(x))/sd(x)
    return(z)
    }

center_data <- function(x){
    z <- (x - mean(x))
    return(z)
    }

INT_transform <- function(x){
    int <- qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
    return(int)
    }

comp_cv <- function(x){sd(x) / mean(x)}

shap_test <- function(x){
    if (length(unique(x)) > 1) {
        return(shapiro.test(x)$p.value)
    } else {
        return(NA)
    }
}

illumina_array_preprocess <- function(exp, gte, normalize = TRUE){
    # Leave in only probes for which there is empirical probe mapping info and convert data into matrix
    emp <- fread(args$emp_probe_mapping, keepLeadingZeros = TRUE)
    emp <- emp[, c(1, 2), with = FALSE]
    colnames(exp)[1] <- "Probe"

    exp$Probe <- as.character(exp$Probe)
    emp$Probe <- as.character(emp$Probe)

    #fwrite(exp, "1_debug_raw_readin.txt", sep = "\t", quote = FALSE)

    exp <- merge(exp, emp, by = "Probe")
    exp <- as.data.frame(exp)
    rownames(exp) <- exp[, ncol(exp)]
    exp <- exp[, -ncol(exp)]
    exp <- exp[, -1]
    exp <- as.matrix(exp)

    #fwrite(exp, "2_debug_raw_ProbesReplaced.txt", sep = "\t", quote = FALSE)

    # Remove samples which are not in the gte or in genotype data
    gte <- fread(args$genotype_to_expression_linking, header = FALSE, keepLeadingZeros = TRUE, colClasses = "character")
    gte$V1 <- as.character(gte$V1)
    gte$V2 <- as.character(gte$V2)

    gte <- gte[gte$V2 %in% as.character(colnames(exp)),]

    message(paste(nrow(gte), "overlapping samples in the gte file AND genotype data."))

    exp <- exp[, colnames(exp) %in% gte$V2]

    #fwrite(exp, "3_debug_raw_ProbesReplaced_SamplesFilteredBasedOnGte.txt", sep = "\t", quote = FALSE)

    exp <- exp[, base::order(colnames(exp))]
    gte <- gte[base::order(gte$V2), ]

    if(!all(colnames(exp) == gte$V2)){stop("Something went wrong in matching genotype and expression IDs. Please debug!")}
    colnames(exp) <- gte$V1

    gene_variance <- data.frame(gene = rownames(exp), gene_variance = apply(exp, 1, var))
    exp <- exp[!rownames(exp) %in% gene_variance[gene_variance$gene_variance == 0, ]$gene, ]

    if (normalize == TRUE){
        # quantile normalization
        exp_n <- normalize.quantiles(exp, copy = FALSE)
        colnames(exp_n) <- colnames(exp)
        rownames(exp_n) <- rownames(exp)
    } else {exp_n <- exp}

    #fwrite(exp, "4_debug_raw_ProbesReplaced_SamplesFilteredBasedOnGteNormalised.txt", sep = "\t", quote = FALSE)

    message(paste(ncol(exp_n), "samples in normalised expression matrix."))

    return(exp_n)
}

RNAseq_preprocess <- function(exp, gte, normalize = TRUE, gene_inclusion=NULL){
    colnames(exp)[1] <- "Probe"

    exp$Probe <- as.character(exp$Probe)

    exp <- as.data.frame(exp)
    rownames(exp) <- exp[, 1]
    exp <- exp[, -1]
    exp <- abs(as.matrix(exp))

    # Remove samples which are not in the gte or in genotype data
    gte <- fread(args$genotype_to_expression_linking, header = FALSE, colClasses = "character")
    gte$V1 <- as.character(gte$V1)
    gte$V2 <- as.character(gte$V2)

    gte <- gte[gte$V2 %in% as.character(colnames(exp)),]

    message(paste(nrow(gte), "overlapping samples in the gte file AND genotype data."))

    exp <- exp[, colnames(exp) %in% gte$V2]
    exp <- exp[, base::order(colnames(exp))]
    gte <- gte[base::order(gte$V2), ]

    if(!all(colnames(exp) == gte$V2)){stop("Something went wrong in matching genotype and expression IDs. Please debug!")}

    colnames(exp) <- gte$V1

    # Remove genes with no variance
    gene_variance <- data.frame(gene = rownames(exp), gene_variance = apply(exp, 1, var))
    exp <- exp[!rownames(exp) %in% gene_variance[gene_variance$gene_variance == 0, ]$gene, ]

    # Remove genes with CPM>0.5 in less than 1% of samples
    exp_keep <- DGEList(counts = exp)
    n_samples_with_cpm_threshold_passed <- rowSums(cpm(exp_keep, log = FALSE) > 0.5)
    keep <- n_samples_with_cpm_threshold_passed >= round(ncol(exp) / 100, 0)

    message(paste(sum(keep), "genes has CPM>0.5 in at least 1% of samples."))

    if (all(!is.null(gene_inclusion)) && length(gene_inclusion) > 0) {
      missed_genes <- gene_inclusion[(gene_inclusion %in% names(keep[keep == FALSE]))]

      if (length(missed_genes) > 0) {
        warning(sprintf(
          "Including %d genes that do not pass cpm cutoff. printing number of samples passing cpm cutoff below.",
          length(missed_genes)))

        print(n_samples_with_cpm_threshold_passed[missed_genes])

      }
      keep[missed_genes] <- TRUE
    }

    exp <- exp[rownames(exp) %in% names(keep[keep == TRUE]), ]

    message(sprintf("Total number of genes that are included: %d", nrow(exp)))

    if (normalize == TRUE){
        # TMM-normalized counts
        exp_n <- DGEList(counts = exp)
        exp_n <- calcNormFactors(exp_n, method = "TMM")
        exp_n <- cpm(exp_n, log = FALSE)
    } else {exp_n <- exp}

    message(paste(ncol(exp_n), "samples in normalised expression matrix."))

    return(exp_n)
}

Affy_preprocess <- function(exp, gte){
    # Leave in only probes for which there is empirical probe mapping info and convert data into matrix
    message("For Affymetrix arrays we assume that input expression matrix is already appropriately normalised and transformed.")
    emp <- fread(args$emp_probe_mapping, keepLeadingZeros = TRUE)
    emp <- emp[, c(1, 2), with = FALSE]
    colnames(exp)[1] <- "Probe"

    exp$Probe <- as.character(exp$Probe)
    emp$Probe <- as.character(emp$Probe)

    exp <- merge(exp, emp, by = "Probe")
    exp <- as.data.frame(exp)
    rownames(exp) <- exp[, ncol(exp)]
    exp <- exp[, -ncol(exp)]
    exp <- exp[, -1]
    exp <- as.matrix(exp)

    # Remove samples which are not in the gte or in genotype data
    gte <- fread(args$genotype_to_expression_linking, header = FALSE, keepLeadingZeros = TRUE, colClasses = "character")
    gte$V1 <- as.character(gte$V1)
    gte$V2 <- as.character(gte$V2)

    gte <- gte[gte$V2 %in% as.character(colnames(exp)),]

    message(paste(nrow(gte), "overlapping samples in the gte file AND genotype data."))

    exp <- exp[, colnames(exp) %in% gte$V2]
    exp <- exp[, base::order(colnames(exp))]
    gte <- gte[base::order(gte$V2), ]

    if(!all(colnames(exp) == gte$V2)){stop("Something went wrong in matching genotype and expression IDs. Please debug!")}

    colnames(exp) <- gte$V1

    # Remove genes with no variance
    gene_variance <- data.frame(gene = rownames(exp), gene_variance = apply(exp, 1, var))
    exp <- exp[!rownames(exp) %in% gene_variance[gene_variance$gene_variance == 0, ]$gene, ]

    # No normalisation done, it assumes that you have already done this.
    message(paste(ncol(exp_n), "samples in normalised expression matrix."))

    return(exp_n)
}

RNAseq_HGNC_preprocess <- function(exp, gte,  normalize = TRUE, gene_inclusion = NULL){

    # Leave in only probes for which there is empirical probe mapping info and convert data into matrix
    emp <- fread(args$emp_probe_mapping, keepLeadingZeros = TRUE)
    emp <- emp[, c(1, 2), with = FALSE]
    colnames(exp)[1] <- "Probe"

    exp$Probe <- as.character(exp$Probe)
    emp$Probe <- as.character(emp$Probe)

    exp <- merge(exp, emp, by = "Probe")
    exp <- as.data.frame(exp)
    rownames(exp) <- exp[, ncol(exp)]
    exp <- exp[, -ncol(exp)]
    exp <- exp[, -1]
    exp <- abs(as.matrix(exp))

    # Remove samples which are not in the gte or in genotype data
    gte <- fread(args$genotype_to_expression_linking, header = FALSE, colClasses = "character")
    gte$V1 <- as.character(gte$V1)
    gte$V2 <- as.character(gte$V2)

    message(paste(nrow(gte), "overlapping samples in the gte file AND genotype data."))

    exp <- exp[, colnames(exp) %in% gte$V2]
    exp <- exp[, base::order(colnames(exp))]
    gte <- gte[base::order(gte$V2), ]

    if(!all(colnames(exp) == gte$V2)){stop("Something went wrong in matching genotype and expression IDs. Please debug!")}

    colnames(exp) <- gte$V1

    # Remove genes with no variance
    gene_variance <- data.frame(gene = rownames(exp), gene_variance = apply(exp, 1, var))
    exp <- exp[!rownames(exp) %in% gene_variance[gene_variance$gene_variance == 0, ]$gene, ]

    # Remove genes with CPM>0.5 in less than 1% of samples
    exp_keep <- DGEList(counts = exp)
    n_samples_with_cpm_threshold_passed <- rowSums(cpm(exp_keep, log = FALSE) > 0.5)
    keep <- n_samples_with_cpm_threshold_passed >= round(ncol(exp) / 100, 0)

    message(paste(sum(keep), "genes has CPM>0.5 in at least 1% of samples."))

    if (all(!is.null(gene_inclusion)) && length(gene_inclusion) > 0) {
      missed_genes <- gene_inclusion[(gene_inclusion %in% names(keep[keep == FALSE]))]

      if (length(missed_genes) > 0) {
        warning(sprintf(
          "Including %d genes that do not pass cpm cutoff. printing number of samples passing cpm cutoff below.",
          length(missed_genes)))

        print(n_samples_with_cpm_threshold_passed[missed_genes])

      }
      keep[missed_genes] <- TRUE
    }

    exp <- exp[rownames(exp) %in% names(keep[keep == TRUE]), ]

    message(sprintf("Total number of genes that are included: %d", nrow(exp)))

    if (normalize == TRUE){
        # TMM-normalized counts
        exp_n <- DGEList(counts = exp)
        exp_n <- calcNormFactors(exp_n, method = "TMM")
        exp_n <- cpm(exp_n, log = FALSE)
    } else {exp_n <- exp}

    message(paste(ncol(exp_n), "samples in normalised expression matrix."))

    return(exp_n)
}

exp_summary <- function(x){

    per_gene_mean <- apply(x, 1, mean)
    per_gene_median <- apply(x, 1, median)
    per_gene_min <- apply(x, 1, min)
    per_gene_max <- apply(x, 1, max)
    per_gene_sd <- apply(x, 1, sd)
    nr_values <- apply(x, 1, function(x) length(x))
    unique_values <- apply(x, 1, function(x) length(unique(x)))
    zero_values <- apply(x, 1, function(x) length(x[x == 0]))
    per_gene_shapiro <- apply(x, 1, shap_test)

    gene_summary <- data.table(gene = rownames(x), 
    mean = per_gene_mean,
    median = per_gene_median,
    min = per_gene_min,
    max = per_gene_max,
    sd = per_gene_sd,
    nr_values = nr_values,
    nr_unique_values = unique_values,
    nr_of_zero_values = zero_values,
    shapiro_P = per_gene_shapiro)

    return(gene_summary)
}


############
# Analysis #
############
# Read in raw expression matrix
and <- fread(args$expression_matrix, header = TRUE,
             keepLeadingZeros = TRUE)
colnames(and) <- as.character(colnames(and))
colnames(and)[1] <- "Feature"
and$Feature <- as.character(and$Feature)
message(paste("Initially:", nrow(and), "genes/probes and ", ncol(and), "samples"))

summary_table <- data.table(Stage = "Unprocessed matrix", Nr_of_features = nrow(and), Nr_of_samples = ncol(and))

# Remove samples which are not in the gte or in genotype data
gte <- fread(args$genotype_to_expression_linking, header = FALSE,
             keepLeadingZeros = TRUE, colClasses = "character")
gte$V1 <- as.character(gte$V1)
gte$V2 <- as.character(gte$V2)

filtered_expr <- read.delim(args$norm_expression_matrix, header = TRUE, row.names = 1, sep = "\t", as.is = T, check.names = F)

samples_to_include <- row.names(filtered_expr)
if (nrow(filtered_expr) > ncol(filtered_expr)) { # if samples in columns
   samples_to_include <- colnames(filtered_expr)
}

gte <- gte[gte$V1 %in% samples_to_include, ]


and <- and[, colnames(and) %in% c("Feature", gte$V2), with = FALSE]
message(paste("After keeping only the same samples as in the eQTLGen normalized matrix:", nrow(and), "genes/probes and ", ncol(and), "samples"))

if (nrow(and) < 100){stop("Less than 100 samples overlap with QCd genotype data!")}
if (!args$platform %in% c("HT12v3", "HT12v4", "HuRef8", "RNAseq", "AffyU219", "AffyHumanExon", "RNAseq_HGNC")){stop("Platform has to be one of HT12v3, HT12v4, HuRef8, RNAseq, AffyU291, AffyHuEx, RNAseq_HGNC")}

message("Running normalization")
if (args$platform %in% c("HT12v3", "HT12v4", "HuRef8")){
  and_p <- illumina_array_preprocess(and, args$genotype_to_expression_linking)
}
if (args$platform %in% c("RNAseq")){
  and_p <- RNAseq_preprocess(and, args$genotype_to_expression_linking)
}
if (args$platform %in% c("RNAseq_HGNC")){
  and_p <- RNAseq_HGNC_preprocess(and, args$genotype_to_expression_linking)
}
if (args$platform %in% c("AffyU219", "AffyHumanExon")){
  and_p <- Affy_preprocess(and, args$genotype_to_expression_linking)
}

summary_table_temp <- data.table(Stage = "Samples which overlap filtered expression table from eQTLGen", Nr_of_features = nrow(and), Nr_of_samples = ncol(and) - 1)
summary_table <- rbind(summary_table, summary_table_temp)

fwrite(as.data.frame(and_p), paste0(args$output, "/exp_data_preprocessed.txt"), sep = "\t", quote = FALSE, row.names = T)


# Summary statistics ----
message("Calculating descriptive summary statistics for every gene...")
# Per gene, calculate mean, median, min, max, sd, and shapiro test P-value
# before preprocessing
message("Before normalisation.")
# Remove sample outliers

message(paste("tmp ", dim(and)))
if (args$platform %in% c("HT12v3", "HT12v4", "HuRef8")){
and_pp <- illumina_array_preprocess(and, args$genotype_to_expression_linking,  normalize = FALSE)
} else if (args$platform %in% c("RNAseq")){
and_pp <- RNAseq_preprocess(and, args$genotype_to_expression_linking, normalize = FALSE)
} else if (args$platform %in% c("RNAseq_HGNC")){
and_pp <- RNAseq_HGNC_preprocess(and, args$genotype_to_expression_linking, normalize = FALSE)
} else if (args$platform %in% c("AffyU219", "AffyHumanExon")){
and_pp <- Affy_preprocess(and, args$genotype_to_expression_linking)
}

message(paste("tmp ", dim(and_pp)))

gene_summary <- exp_summary(and_pp)

fwrite(gene_summary, paste0(args$output, "/raw_gene_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# on fully processed expression matrix
message("After normalization.")

gene_summary <- exp_summary(and_p)
fwrite(gene_summary, paste0(args$output, "/processed_gene_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Write out summary table about features and samples
colnames(summary_table) <- c("Stage", "Nr. of features", "Nr. of samples")
fwrite(summary_table, paste0(args$output, "/summary_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
