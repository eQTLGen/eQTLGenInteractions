library(data.table)
library(preprocessCore)
library(edgeR)
library(ggplot2)
library(optparse)
library(patchwork)
library(MASS)
library(dplyr)

setDTthreads(8)

# Argument parser
option_list <- list(
    make_option(c("-e", "--expression_matrix"), type = "character",
    help = "Unprocessed gene expression matrix from array or RNA-seq experiment. Samples in columns, genes/probes in the rows."),
    make_option(c("-l", "--genotype_to_expression_linking"), type = "character",
    help = "Genotype-to-expression linking file. No header. First column- sample IDs in the genotype file, second column- corresponding sample IDs in the gene expression matrix."),
    make_option(c("-p", "--platform"), type = "character",
    help = "Gene expression platform. This determines the normalization method and replaces probes with best-matching genes based on empirical probe mapping. One of: HT12v3, HT12v4, HuRef8, RNAseq, AffyU219, AffyHumanExon, RNAseq_HGNC."),
    make_option(c("-m", "--emp_probe_mapping"), type = "character",
    help = "Empirical probe matching file. Used to link the best array probe to each blood-expressed gene."),
    make_option(c("-s", "--sd"), type = "double", default = 4,
    help = "Standard deviation threshold for removing expression samples. By default, samples away 4 SDs from the median of PC1 are removed."),
    make_option(c("-c", "--contamination_area"), type = "double", default = 30,
    help = "Area that marks likely contaminated samples based on sex-chromosome gene expression. Must be an angle between 0 and 90. The angle represents the total area around the y = x function."),
    make_option(c("-a", "--contamination_slope"), type = "double", default = 45,
    help = "Angle that is used to discriminate based on expression inferred sex. Increase to make it less steep, decrease to make it steeper"),
    make_option(c("-i", "--sex_info"), type = "character",
    help = "File with sex information. Plink2 --check-sex filtered output."),
    make_option(c("-f", "--geno_filter"), type = "character",
    help = "File with filtered genotype samples. Plink fam file."),
    make_option(c("-o", "--output"), type = "character",
    help = "Output folder where to put preprocessed data matrix, expression PCs, etc.")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Debug
print(args$expression_matrix)
print(args$genotype_to_expression_linking)
print(args$geno_filter)
print(args$sex_info)
print(args$platform)
print(args$emp_probe_mapping)
print(args$sd)
print(args$output)

# Make output folder structure
dir.create(args$output)
dir.create(paste0(args$output, "/exp_plots"))
dir.create(paste0(args$output, "/exp_data_QCd"))
dir.create(paste0(args$output, "/exp_PCs"))
dir.create(paste0(args$output, "/exp_data_summary"))

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

illumina_array_preprocess <- function(exp, gte, gen, normalize = TRUE){
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

    geno_fam <- fread(args$sex_info, header = TRUE, keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))
    gte <- gte[gte$V1 %in% geno_fam$IID,]
    gte <- gte[gte$V2 %in% as.character(colnames(exp)),]

    message(paste(nrow(gte), "overlapping samples in the gte file AND genotype data."))

    exp <- exp[, colnames(exp) %in% gte$V2]

    #fwrite(exp, "3_debug_raw_ProbesReplaced_SamplesFilteredBasedOnGte.txt", sep = "\t", quote = FALSE)

    exp <- exp[, base::order(colnames(exp))]
    gte <- gte[base::order(gte$V2), ]

    if(!all(colnames(exp) == gte$V2)){stop("Something went wrong in matching genotype and expression IDs. Please debug!")}
    colnames(exp) <- gte$V1

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

RNAseq_preprocess <- function(exp, gte, gen, normalize = TRUE, gene_inclusion=NULL){
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

    geno_fam <- fread(args$sex_info, header = TRUE, keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))
    gte <- gte[gte$V1 %in% geno_fam$IID,]
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

Affy_preprocess <- function(exp, gte, gen){
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

    geno_fam <- fread(args$sex_info, header = TRUE, keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))
    gte <- gte[gte$V1 %in% geno_fam$IID,]
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

RNAseq_HGNC_preprocess <- function(exp, gte, gen, normalize = TRUE, gene_inclusion = NULL){

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

    geno_fam <- fread(args$sex_info, header = TRUE, keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))
    gte <- gte[gte$V1 %in% geno_fam$IID,]
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


rnaseq_cpm_summary <- function(exp, gte){
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

  geno_fam <- fread(args$sex_info, header = TRUE, keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))
  gte <- gte[gte$V1 %in% geno_fam$IID,]
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

  # Remove genes with CPM<0.5 in less than 1% of samples
  exp_keep <- DGEList(counts = exp)
  n_samples_with_cpm_threshold_passed <- rowSums(cpm(exp_keep, log = FALSE) > 0.5)
  exp_keep <- calcNormFactors(exp_keep)
  n_samples_with_cpm_threshold_passed_normalized <- rowSums(cpm(exp_keep, log = FALSE) > 0.5)
  threshold <- round(ncol(exp) / 100, 0)

  genes_cpm_summary <-
    data.frame(n_samples_with_cpm_threshold_passed=n_samples_with_cpm_threshold_passed,
               n_samples_with_cpm_threshold_passed_normalized=n_samples_with_cpm_threshold_passed_normalized)

  genes_cpm_summary$n_samples_threshold <- threshold

  return(genes_cpm_summary)
}


IterativeOutlierDetection <- function(
  input_exp, sd_threshold = 1,
  platform = c("HT12v3", "HT12v4", "HuRef8", "RNAseq", "AffyU291", "AffyHuEx", "RNAseq_HGNC"),
  gene_inclusion = NULL) {

  and <- input_exp
  list_ggplots <- list()
  summary_pcs <- list()
  outliers <- c()
  platform <- match.arg(platform)
  message(paste0("Expression platform is ", platform, ", preprocessing the data..."))

  it_round <- 0
  plot_it_round <- 0
  nr_outliers <- 1

  message(paste("Removing samples which deviate more than", sd_threshold, "SDs from the mean values of PC1 or PC2."))
  message("Starting iterative outlier detection...")

  while (nr_outliers > 0) {

    if (platform %in% c("HT12v3", "HT12v4", "HuRef8")){
      and_p <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
      #fwrite(and_p, "5_debug_normalised_exp.txt", sep = "\t", quote = FALSE)
      and_p <- log2(and_p)
      #fwrite(and_p, "6_debug_normalised_log2_exp.txt", sep = "\t", quote = FALSE)
      #and_p <- apply(and_p, 1, INT_transform)
      #and_p <- t(and_p)
      #and_p <- apply(and_p, 1, center_data)
      #and_p <- t(and_p)
    } else if(platform %in% c("RNAseq")){
      and_p <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples,
                                 gene_inclusion = gene_inclusion)
      and_p <- log2(and_p + 0.25)
      #and_p <- apply(and_p, 1, INT_transform)
      #and_p <- t(and_p)
      #and_p <- apply(and_p, 1, center_data)
      #and_p <- t(and_p)
    } else if(platform %in% c("RNAseq_HGNC")){
      and_p <- RNAseq_HGNC_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples,
                                 gene_inclusion = gene_inclusion)
      and_p <- log2(and_p + 0.25)
      #and_p <- apply(and_p, 1, INT_transform)
      #and_p <- t(and_p)
      #and_p <- apply(and_p, 1, center_data)
      #and_p <- t(and_p)
    } else if(platform %in% c("AffyU219", "AffyHumanExon")){
      and_p <- Affy_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
      #and_p <- apply(and_p, 1, INT_transform)
      #and_p <- t(and_p)
      #and_p <- apply(and_p, 1, center_data)
      #and_p <- t(and_p)
    }
    message("Data preprocessed!")

    it_round <- it_round + 1
    pcs <- prcomp(t(and_p), center = FALSE, scale. = FALSE)
    message("PCA calculation finished!")

    PCs <- as.data.table(pcs$x)
    PCs$sample <- rownames(pcs$x)

    importance <- pcs$sdev^2 / sum(pcs$sdev^2)
    summary_pcs[[it_round]] <- data.table(PC = paste0("PC", 1:50), explained_variance = importance[1:50])

    PCs$outlier <- "no"
    PCs1_mean <- mean(PCs$PC1)
    PCs1_sd <- sd(PCs$PC1)
    PCs[(PCs$PC1 > PCs1_mean + sd_threshold * PCs1_sd) | (PCs$PC1 < PCs1_mean - sd_threshold * PCs1_sd)]$outlier <- "yes"

    PCs2_mean <- mean(PCs$PC2)
    PCs2_sd <- sd(PCs$PC2)
    PCs[(PCs$PC2 > PCs2_mean + sd_threshold * PCs2_sd) | (PCs$PC2 < PCs2_mean - sd_threshold * PCs2_sd)]$outlier <- "yes"

    plot_it_round <- plot_it_round + 1
    list_ggplots[[plot_it_round]] <- ggplot(PCs, aes(x = PC1, y = PC2, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))
    
    plot_it_round <- plot_it_round + 1
    list_ggplots[[plot_it_round]] <- ggplot(PCs, aes(x = PC3, y = PC4, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))
    
    plot_it_round <- plot_it_round + 1
    list_ggplots[[plot_it_round]] <- ggplot(PCs, aes(x = PC5, y = PC6, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))
    
    plot_it_round <- plot_it_round + 1
    list_ggplots[[plot_it_round]] <- ggplot(PCs, aes(x = PC7, y = PC8, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))

    non_outliers <- PCs[!PCs$outlier %in% c("yes"), ]$sample

    non_outliers <- gte[as.character(gte$V1) %in% non_outliers, ]$V2
  
    # Remove outlier samples from the unprocessed data
    and <- and[, colnames(and) %in% c(non_outliers, "Feature"), with = FALSE]

    nr_outliers <- length(PCs[PCs$outlier %in% c("yes"), ]$sample)

    if (nr_outliers > 0) {
      message(paste0("Iteration round ", it_round, ". Removed ", nr_outliers, " outlier(s). Re-processing the data and running another round of PCA."))
      summary_table_temp <- data.table(Stage = paste("After removal of expression outliers in ", it_round, " iteration."), Nr_of_features = nrow(and), Nr_of_samples = ncol(and))
      summary_table <- rbind(summary_table, summary_table_temp)

      outliers <- c(outliers, PCs[PCs$outlier %in% c("yes"), ]$sample)
    } else if (nr_outliers == 0) {
      message(paste0("Iteration round ", it_round, ". No outliers detected. Finalizing interactive outlier detection."))
    }
    if (it_round == 1){
      PCs_it1 <- PCs
    } else {PCs_it1[PCs_it1$sample %in% outliers, ]$outlier <- "yes"}
  }

  return(list(exp_mat = and_p, plots = list_ggplots, summary_pcs = summary_pcs, PCs_first = PCs_it1))
}

ExpressionBasedSampleSwapIdentification <- function(and, summary_table) {
  # Visualise the expression of X-specific and Y-specific genes
  # Read in emp probe mapping file with chr information
  emp_probe_mapping <- fread(args$emp_probe_mapping, keepLeadingZeros = TRUE)

  sex_specific_genes <- c(
    "ENSG00000229807",
    emp_probe_mapping[emp_probe_mapping$chromosome_name == "Y", ]$Ensembl)

  # Final re-process, re-calculate PCs, re-visualise and write out
  if (args$platform %in% c("HT12v3", "HT12v4", "HuRef8")){
    and_pp <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
  } else if (args$platform %in% c("RNAseq")){
    and_pp <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, gene_inclusion = sex_specific_genes)
  } else if (args$platform %in% c("RNAseq_HGNC")){
    and_pp <- RNAseq_HGNC_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, gene_inclusion = sex_specific_genes)
  } else if (args$platform %in% c("AffyU219", "AffyHumanExon")){
    and_pp <- Affy_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
  }

  xist <- and_pp[rownames(and_pp) == "ENSG00000229807", ]
  y_genes <- emp_probe_mapping[emp_probe_mapping$chromosome_name == "Y", ]$Ensembl
  y_genes <- and_pp[rownames(and_pp) %in% y_genes, ]

  #fwrite(xist, "7_debug_xist_genes.txt")
  #fwrite(y_genes, "7_debug_y_chr_genes.txt")


  nr_of_y_genes <- nrow(y_genes)
  message(paste("Nr. of Y genes:", nr_of_y_genes))

  if (nrow(y_genes) > 0){

    y_mean <- apply(y_genes, 2, mean)
    min_y_mean <- min(y_mean)

    message(print(min_y_mean))

    if (length(xist) == 0){
      message("XIST not detected, plotting the mean of X chr genes instead!")
      xist <- and_pp[rownames(and_pp) %in% emp_probe_mapping[emp_probe_mapping$chromosome_name == "X", ]$Ensembl, ]
      xist_mean <- apply(xist, 2, mean)
      min_xist <- min(xist_mean)
      nr_of_x_genes <- nrow(xist)
      y_genes <- data.frame(sample = colnames(y_genes), xist = xist_mean - min_xist, y_genes = y_mean - min_y_mean)
      xist_missing <- TRUE
    } else {
      min_xist <- min(xist)
      y_genes <- data.frame(sample = colnames(y_genes), xist = xist - min_xist, y_genes = y_mean - min_y_mean)
      xist_missing <- FALSE
    }

    geno_fam_f <- geno_fam[, c(2, 4), with = FALSE]
    colnames(geno_fam_f) <- c("sample", "Sex")
    geno_fam_f$Sex <- as.character(geno_fam_f$Sex)

    y_genes <- merge(y_genes, geno_fam_f, by = "sample")
    max_exp <- max(y_genes$y_genes, y_genes$xist)

    lower_slope <- tan((args$contamination_slope - args$contamination_area / 2) / 180*pi)
    upper_slope <- tan((args$contamination_slope + args$contamination_area / 2) / 180*pi)
    middle_slope <- tan(args$contamination_slope / 180*pi)

    y_genes$expressionSexNaive <- case_when(
      y_genes$y_genes > y_genes$xist * middle_slope ~ 1,
      y_genes$y_genes < y_genes$xist * middle_slope ~ 2
    )

    x_expression_median <- median(y_genes[y_genes$Sex == 1 & y_genes$expressionSexNaive == 1, "xist"])
    y_expression_median <- median(y_genes[y_genes$Sex == 2 & y_genes$expressionSexNaive == 2, "y_genes"])

    x_expression_min <- min(y_genes[y_genes$Sex == 2 & y_genes$expressionSexNaive == 2, "xist"])
    y_expression_min <- min(y_genes[y_genes$Sex == 1 & y_genes$expressionSexNaive == 1, "y_genes"])

    y_max <- max(y_genes$y_genes)
    x_max <- max(y_genes$xist)

    x_anchor <- min(x_expression_median, x_expression_min - 0.01 * x_max)
    y_anchor <- min(y_expression_median, y_expression_min - 0.01 * y_max)

    y_genes$xist_corrected <- y_genes$xist - x_anchor

    y_genes$contaminated <- case_when(
      (y_genes$y_genes > ((y_genes$xist_corrected) * lower_slope + y_anchor)
        & y_genes$y_genes < ((y_genes$xist_corrected) * upper_slope + y_anchor)) ~ "yes",
      TRUE ~ "no"
    )

    y_genes$expressionSex <- case_when(
      (y_genes$y_genes < ((y_genes$xist_corrected) * middle_slope + y_anchor)) ~ 2,
      (y_genes$y_genes > ((y_genes$xist_corrected) * middle_slope + y_anchor)) ~ 1
    )

    y_genes$mismatch <- case_when(
      y_genes$Sex == 0 ~ "unknown",
      y_genes$expressionSex == y_genes$Sex ~ "no",
      y_genes$expressionSex != y_genes$Sex ~ "yes"
    )

    y_genes$status <- case_when(
      y_genes$contaminated == "yes" & y_genes$mismatch == "yes" ~ "Contaminated and\nsex mismatch",
      y_genes$contaminated == "yes" ~ "Likely contaminated",
      y_genes$mismatch == "yes" ~ "Sex mismatch",
      TRUE ~ "Passed"
    )

    y_genes_zoom <- max(y_genes %>% ungroup() %>%
      filter(contaminated == "yes" | mismatch == "yes") %>%
      summarise(max_y = max(y_genes), max_x = max(xist)))

    exclusion_zone <- tibble(x = c(x_expression_median, max_exp)) %>%
      mutate(lower_bound = (x - x_expression_median) * lower_slope + y_expression_median,
             upper_bound = (x - x_expression_median) * upper_slope + y_expression_median,
             middle_line = (x - x_expression_median) * middle_slope + y_expression_median)

    if (xist_missing == FALSE){
      base_plot <- ggplot(data = exclusion_zone, aes(x = x, ymin = lower_bound, ymax = upper_bound)) +
        geom_ribbon(alpha = 0.2) +
        geom_line(aes(x = x, y = middle_line), linetype = 2, colour = "blue") +
        geom_point(data = y_genes, inherit.aes = F, aes(col = status, shape = Sex, x = xist, y = y_genes)) +
        scale_colour_manual(
          values = alpha(c("Passed" = "black",
                           "Likely contaminated" = "orange",
                           "Sex mismatch" = "red",
                           "Contaminated and\nsex mismatch" = "darkmagenta"),
                         0.5),
          name = "Passed checks") +
        coord_cartesian(ylim = c(0, max_exp), xlim = c(0, max_exp)) +
        theme_bw() + ylab(paste0("mean of Y genes - min(mean of Y genes)\n(n=", nr_of_y_genes, ")")) + xlab("XIST - min(XIST)")

      ggsave(paste0(args$output, "/exp_plots/SexSpecificGenesXIST.png"), height = 5, width = 7, units = "in", dpi = 300, type = "cairo")
      ggsave(paste0(args$output, "/exp_plots/SexSpecificGenesXIST.pdf"), height = 5, width = 7, units = "in", dpi = 300)

      naive_plot <- ggplot(data = exclusion_zone, aes(x = x, ymin = lower_bound, ymax = upper_bound)) +
        geom_ribbon(alpha = 0.2) +
        geom_line(aes(x = x, y = middle_line), linetype = 2, colour = "blue") +
        geom_point(data = y_genes, inherit.aes = F, aes(col = expressionSexNaive, shape = Sex, x = xist, y = y_genes)) +
        coord_cartesian(ylim = c(0, y_genes_zoom), xlim = c(0, y_genes_zoom)) +
        theme_bw() + ylab(paste0("mean of Y genes - min(mean of Y genes)\n(n=", nr_of_y_genes, ")")) + xlab("XIST - min(XIST)")

      ggsave(paste0(args$output, "/exp_plots/SexSpecificGenesXIST_naive.png"), height = 5, width = 7, units = "in", dpi = 300, type = "cairo")
      ggsave(paste0(args$output, "/exp_plots/SexSpecificGenesXIST_naive.pdf"), height = 5, width = 7, units = "in", dpi = 300)

      zoomed_plot <- ggplot(data = exclusion_zone, aes(x = x, ymin = lower_bound, ymax = upper_bound)) +
        geom_ribbon(alpha = 0.2) +
        geom_line(aes(x = x, y = middle_line), linetype = 2, colour = "blue") +
        geom_point(data = y_genes, inherit.aes = F, aes(col = status, shape = Sex, x = xist, y = y_genes)) +
        scale_colour_manual(
          values = alpha(c("Passed" = "black",
                           "Likely contaminated" = "orange",
                           "Sex mismatch" = "red",
                           "Contaminated and\nsex mismatch" = "darkmagenta"),
                         0.5),
          name = "Passed checks") +
        coord_cartesian(ylim = c(0, y_genes_zoom), xlim = c(0, y_genes_zoom)) +
        theme_bw() + ylab(paste0("mean of Y genes - min(mean of Y genes)\n(n=", nr_of_y_genes, ")")) + xlab("XIST - min(XIST)")

      ggsave(paste0(args$output, "/exp_plots/SexSpecificGenesXIST_zoomed.png"), height = 5, width = 7, units = "in", dpi = 300, type = "cairo")
      ggsave(paste0(args$output, "/exp_plots/SexSpecificGenesXIST_zoomed.pdf"), height = 5, width = 7, units = "in", dpi = 300)

    } else {
      base_plot <- ggplot(data = exclusion_zone, aes(x = x, ymin = lower_bound, ymax = upper_bound)) +
        geom_ribbon(alpha = 0.2) +
        geom_line(aes(x = x, y = middle_line), linetype = 2, colour = "blue") +
        geom_point(data = y_genes, inherit.aes = F, aes(col = status, shape = Sex, x = xist, y = y_genes)) +
        scale_colour_manual(
          values = alpha(c("Passed" = "black",
                           "Likely contaminated" = "orange",
                           "Sex mismatch" = "red",
                           "Contaminated and\nsex mismatch" = "darkmagenta"),
                         0.5),
          name = "Passed checks") +
        coord_cartesian(ylim = c(0, max_exp), xlim = c(0, max_exp)) +
        theme_bw() + ylab(paste0("mean of Y genes - min(mean of Y genes)\n(n=", nr_of_y_genes, ")")) + xlab(paste0("mean of X genes - min(mean of X genes)\n(n=", nr_of_x_genes, ")"))

      ggsave(paste0(args$output, "/exp_plots/SexSpecificGenes.png"), height = 5, width = 7, units = "in", dpi = 300, type = "cairo")
      ggsave(paste0(args$output, "/exp_plots/SexSpecificGenes.pdf"), height = 5, width = 7, units = "in", dpi = 300)

      zoomed_plot <- ggplot(data = exclusion_zone, aes(x = x, ymin = lower_bound, ymax = upper_bound)) +
        geom_ribbon(alpha = 0.2) +
        geom_line(aes(x = x, y = middle_line), linetype = 2, colour = "blue") +
        geom_point(data = y_genes, inherit.aes = F, aes(col = status, shape = Sex, x = xist, y = y_genes)) +
        scale_colour_manual(
          values = alpha(c("Passed" = "black",
                           "Likely contaminated" = "orange",
                           "Sex mismatch" = "red",
                           "Contaminated and\nsex mismatch" = "darkmagenta"),
                         0.5),
          name = "Passed checks") +
        coord_cartesian(ylim = c(0, y_genes_zoom), xlim = c(0, y_genes_zoom)) +
        theme_bw() + ylab(paste0("mean of Y genes - min(mean of Y genes)\n(n=", nr_of_y_genes, ")")) + xlab("XIST - min(XIST)")

      ggsave(paste0(args$output, "/exp_plots/SexSpecificGenes_zoomed.png"), height = 5, width = 7, units = "in", dpi = 300, type = "cairo")
      ggsave(paste0(args$output, "/exp_plots/SexSpecificGenes_zoomed.pdf"), height = 5, width = 7, units = "in", dpi = 300)
    }

    # Filter out potential sex mismatches
    and_pp <- and_pp[, colnames(and_pp) %in% y_genes[y_genes$mismatch != "yes", ]$sample]

    # debug mismatched samples
    comb <- merge(gte, y_genes, by.x = "V1", by.y = "sample")
    colnames(comb)[c(1, 2)] <- c("Genotype_SampleID", "Expression_SampleID")
    fwrite(comb, "SexMismatchStatus.txt", sep = "\t", quote = FALSE)

    summary_table_temp <- data.table(Stage = "Samples after removal of sex errors", Nr_of_features = nrow(and_pp), Nr_of_samples = ncol(and_pp))
    summary_table <- rbind(summary_table, summary_table_temp)

    and_pp <- and_pp[, colnames(and_pp) %in% y_genes[y_genes$contaminated != "yes", ]$sample]

    summary_table_temp <- data.table(Stage = "Samples after removal of likely contaminated samples", Nr_of_features = nrow(and_pp), Nr_of_samples = ncol(and_pp))
    summary_table <- rbind(summary_table, summary_table_temp)

  } else {
    message("There is no Y chromosome genes in the expression data, omitting sex check.")
  }
  return(list(and_pp=and_pp, summary_table=summary_table))
}

ExpressionBasedSampleSwapIdentificationNb <- function(sex_specific_genes, xist_col="xist", y_col = "y_genes") {
  # fit Gaussian mixture model
  fit <- Mclust(sex_specific_genes[,c(y_col, xist_col)], G=2, modelNames="VVV")

  # Get means from parameters
  means <- fit$parameters$mean
  sigmas <- fit$parameters$sigma

  # Mean ratios
  mean_ratios <- means['y_genes',] / means['xist',]

  # Create factor of expression derived sex
  expression_derived_sex <- factor(fit$classification, levels = c(1, 2), labels = c("F", "M")[order(mean_ratios)])

  out <- tibble(sex_specific_genes)
  out$uncertainty <- fit$uncertainty
  out$expression_derived_sex <- expression_derived_sex
  out$unconfident <- out$uncertainty > 0.01
  out$Sex <- factor(out$Sex, levels = c(2, 1), labels = c("F", "M"))
  out$mismatch <- out$Sex != out$expression_derived_sex
  out$status <- case_when(out$mismatch & out$unconfident ~ "Contaminated and\nsex mismatch",
                          out$mismatch ~ "Sex mismatch",
                          out$unconfident ~ "Likely contaminated",
                          TRUE ~ "Passed"
  )

  max_exp <- max(sex_specific_genes[,c(y_col, xist_col)])

  base_plot <- ggplot(data = out, aes(x = xist, y = y_genes)) +
    geom_point(aes(col = status, shape = Sex)) +
    scale_colour_manual(
      values = alpha(c("Passed" = "black",
                       "Likely contaminated" = "yellow",
                       "Sex mismatch" = "orange",
                       "Contaminated and\nsex mismatch" = "red"),
                     0.5),
      name = "Passed checks") +
    coord_cartesian(ylim = c(0, max_exp), xlim = c(0, max_exp)) +
    theme_bw() + ylab(paste0("mean of Y genes - min(mean of Y genes)\n(n=", nr_of_y_genes, ")")) + xlab("XIST - min(XIST)")

  ggsave(paste0(args$output, "/exp_plots/SexSpecificGenesXIST2.png"), height = 5, width = 7, units = "in", dpi = 300, type = "cairo")
  ggsave(paste0(args$output, "/exp_plots/SexSpecificGenesXIST2.pdf"), height = 5, width = 7, units = "in", dpi = 300)


  # data.grid <- expand.grid(s.1 = seq(0, max_expression, length.out=100),
  #                          s.2 = seq(0, max_expression, length.out=100))
  #
  # q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = means, sigma = sigmas))
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

geno_fam <- fread(args$sex_info, header = TRUE,
                  keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))
geno_fam$IID <- as.character(geno_fam$IID)

gen_filter <- fread(args$geno_filter, header = FALSE,
                    keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))
gen_filter$V2 <- as.character(gen_filter$V2)
gte <- gte[gte$V1 %in% geno_fam$IID, ]
gte <- gte[gte$V1 %in% gen_filter$V2, ]

and <- and[, colnames(and) %in% c("Feature", gte$V2), with = FALSE]

if (nrow(and) < 100){stop("Less than 100 samples overlap with QCd genotype data!")}

summary_table_temp <- data.table(Stage = "Samples which overlap with QCd genotype info", Nr_of_features = nrow(and), Nr_of_samples = ncol(and) - 1)
summary_table <- rbind(summary_table, summary_table_temp)

if (!args$platform %in% c("HT12v3", "HT12v4", "HuRef8", "RNAseq", "AffyU219", "AffyHumanExon", "RNAseq_HGNC")){stop("Platform has to be one of HT12v3, HT12v4, HuRef8, RNAseq, AffyU291, AffyHuEx, RNAseq_HGNC")}

iterative_outliers <- IterativeOutlierDetection(
  and, sd_threshold = args$sd, platform = args$platform)

# Keep in the original data only non-outlier samples
exp_non_outliers <- colnames(iterative_outliers$exp_mat)
exp_non_outliers <- gte[gte$V1 %in% exp_non_outliers, ]$V2
# Remove outlier samples from MDS
and <- and[, colnames(and) %in% c("Feature", exp_non_outliers), with = FALSE]

summary_table_temp <- data.table(Stage = "After removal of all expression outliers (PCA)", Nr_of_features = nrow(and), Nr_of_samples = ncol(and) - 1)
summary_table <- rbind(summary_table, summary_table_temp)

message("Running MDS")
if (args$platform %in% c("HT12v3", "HT12v4", "HuRef8")){
  and_p <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
  and_p <- log2(and_p)
}
if (args$platform %in% c("RNAseq")){
  and_p <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
  and_p <- log2(and_p + 0.25)
}
if (args$platform %in% c("RNAseq_HGNC")){
  and_p <- RNAseq_HGNC_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
  and_p <- log2(and_p + 0.25)
}
if (args$platform %in% c("AffyU219", "AffyHumanExon")){
  and_p <- Affy_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
}

dist <- cor(and_p, method = "pearson")
mds <- isoMDS(1 - dist, k = 2)
mds <- as.data.frame(mds$points)
colnames(mds) <- paste("MDS coordinate", 1:2)
mds$Sample <- colnames(and_p)

# Find samples deviating from the mean MDS1 and MDS2
mds$outlier <- "no"
mean_mds1 <- mean(mds$`MDS coordinate 1`)
sd_mds1 <- sd(mds$`MDS coordinate 1`)
mean_mds2 <- mean(mds$`MDS coordinate 2`)
sd_mds2 <- sd(mds$`MDS coordinate 2`)
if (nrow(mds[mds$`MDS coordinate 1` > mean_mds1 + args$sd * sd_mds1 |
    mds$`MDS coordinate 1` < mean_mds1 - args$sd * sd_mds1 |
    mds$`MDS coordinate 2` > mean_mds2 + args$sd * sd_mds2 |
    mds$`MDS coordinate 2` < mean_mds2 - args$sd * sd_mds2, ]) > 0){

    mds[mds$`MDS coordinate 1` > mean_mds1 + args$sd * sd_mds1 |
    mds$`MDS coordinate 1` < mean_mds1 - args$sd * sd_mds1 |
    mds$`MDS coordinate 2` > mean_mds2 + args$sd * sd_mds2 |
    mds$`MDS coordinate 2` < mean_mds2 - args$sd * sd_mds2, ]$outlier <- "yes"

}

# Add sex info
sex <- fread(args$sex_info, header = FALSE, keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))
sex <- sex[, c(2, 4), with = FALSE]
colnames(sex) <- c("Sample", "Sex")
sex$Sample <- as.character(sex$Sample)
mds <- merge(mds, sex, by = "Sample")

p <- ggplot(mds, aes(x = `MDS coordinate 1`, `MDS coordinate 2`, colour = outlier, shape = Sex)) +
geom_point(alpha = 0.3) +
theme_bw() +
scale_colour_manual(values = c("no" = "black", "yes" = "red"))

ggsave(paste0(args$output, "/exp_plots/MDS_before.png"), height = 5, width = 6, units = "in", dpi = 300, type = "cairo")
ggsave(paste0(args$output, "/exp_plots/MDS_before.pdf"), height = 5, width = 6, units = "in", dpi = 300)

exp_non_outliers <- gte[gte$V1 %in% mds[mds$outlier == "no", ]$Sample, ]$V2
and <- and[, colnames(and) %in% c("Feature", exp_non_outliers), with = FALSE]

summary_table_temp <- data.table(Stage = "After removal of all expression outliers (MDS)", Nr_of_features = nrow(and), Nr_of_samples = ncol(and) - 1)
summary_table <- rbind(summary_table, summary_table_temp)

#fwrite(and, "7_debug_before_sampleswap.txt", sep = "\t", quote = FALSE)

expression_based_sample_swap_out <- ExpressionBasedSampleSwapIdentification(and, summary_table)

exp_non_outliers <- colnames(expression_based_sample_swap_out$and_pp)
exp_non_outliers <- gte[gte$V1 %in% exp_non_outliers, ]$V2
and <- and[, colnames(and) %in% c("Feature", exp_non_outliers), with = FALSE]

summary_table <- expression_based_sample_swap_out$summary_table

# Final re-process, re-calculate PCs, re-visualise and write out
if (args$platform %in% c("HT12v3", "HT12v4", "HuRef8")){
  and_pp <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
} else if (args$platform %in% c("RNAseq")){
  gene_cpm_summary <- rnaseq_cpm_summary(and, args$genotype_to_expression_linking)
  fwrite(gene_cpm_summary, paste0(args$output, "/exp_data_summary/", "cpm_summary.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
  and_pp <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
} else if (args$platform %in% c("RNAseq_HGNC")){
  gene_cpm_summary <- rnaseq_cpm_summary(and, args$genotype_to_expression_linking)
  fwrite(gene_cpm_summary, paste0(args$output, "/exp_data_summary/", "cpm_summary.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
  and_pp <- RNAseq_HGNC_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
} else if (args$platform %in% c("AffyU219", "AffyHumanExon")){
  and_pp <- Affy_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
}

summary_table_temp <- data.table(Stage = "After filtering genes (CPM>0.5 in >= 1% of samples)", Nr_of_features = nrow(and_pp), Nr_of_samples = ncol(and_pp))
summary_table <- rbind(summary_table, summary_table_temp)

and_p <- and_pp
fwrite(as.data.frame(and_p), paste0(args$output, "/exp_data_QCd/exp_data_preprocessed.txt"), sep = "\t", quote = FALSE, row.names = T)


# Summary statistics ----
message("Calculating descriptive summary statistics for every gene...")
# Per gene, calculate mean, median, min, max, sd, and shapiro test P-value
# before preprocessing
message("Before normalisation.")
# Remove sample outliers
sample_non_outliers <- colnames(and_p)
sample_non_outliers <- gte[gte$V1 %in% sample_non_outliers, ]$V2

and <- and[, colnames(and) %in% c("Feature", sample_non_outliers), with = FALSE]
message(paste("tmp ", dim(and)))
if (args$platform %in% c("HT12v3", "HT12v4", "HuRef8")){
and_pp <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, normalize = FALSE)
} else if (args$platform %in% c("RNAseq")){
and_pp <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, normalize = FALSE)
} else if (args$platform %in% c("RNAseq_HGNC")){
and_pp <- RNAseq_HGNC_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, normalize = FALSE)
} else if (args$platform %in% c("AffyU219", "AffyHumanExon")){
and_pp <- Affy_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
}

message(paste("tmp ", dim(and_pp)))

gene_summary <- exp_summary(and_pp)

fwrite(gene_summary, paste0(args$output, "/exp_data_summary/", "raw_gene_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# on fully processed expression matrix
message("After normalization.")
if (args$platform %in% c("HT12v3", "HT12v4", "HuRef8")){
and_pp <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, normalize = TRUE)
and_pp <- log2(and_pp)
} else if (args$platform %in% c("RNAseq")){
and_pp <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, normalize = TRUE)
and_pp <- log2(and_pp + 0.25)
} else if (args$platform %in% c("RNAseq_HGNC")){
and_pp <- RNAseq_HGNC_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, normalize = TRUE)
and_pp <- log2(and_pp + 0.25)
} else if (args$platform %in% c("AffyU291", "AffyHuEx")){
and_pp <- Affy_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
}

gene_summary <- exp_summary(and_pp)
fwrite(gene_summary, paste0(args$output, "/exp_data_summary/", "processed_gene_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Write out summary table about features and samples
colnames(summary_table) <- c("Stage", "Nr. of features", "Nr. of samples")
fwrite(summary_table, paste0(args$output, "/exp_data_summary/", "summary_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
