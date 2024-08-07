library(data.table)
library(caret)

args <- commandArgs(trailingOnly = TRUE)

covars <- read.delim(args[1], sep = "\t", row.names = 1, as.is = T, check.names = F)
covar_name <- args[2]
norm_expression <- read.delim(args[3], sep = "\t", row.names = 1, as.is = T, check.names = F)
output_prefix <- args[4]

pheno_is_factor <- F
if (length(unique(covars[,covar_name])) < 3) pheno_is_factor <- T

cat ("covariate ", covar_name, " is a factor: ", pheno_is_factor, "\n")

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

shap_test <- function(x){
  if (length(unique(x)) > 1) {
    return(shapiro.test(x)$p.value)
  } else {
    return(NA)
  }
}

if (pheno_is_factor){
  level_num = 1
  for (level in unique(covars[,covar_name])){
    covar_for_bin <- covars[covars[,covar_name] == level,]
    expr_summary_for_bin <- exp_summary(norm_expression[,row.names(covar_for_bin) ])
    cat ("level_num ", level_num, "corresponds to ", level, "\n")

    # Remove covariates with near zero variance
    near_zero_var <- nearZeroVar(covar_for_bin,  uniqueCut = 5)
    covar_for_bin <- covar_for_bin[, -near_zero_var]

    write.table(covar_for_bin, file = paste0(output_prefix, "covariates.", covar_name, "_", level_num, ".txt"), sep = "\t", quote = F, col.names = NA)
    write.table(expr_summary_for_bin, file = paste0(output_prefix, "expression_summary.", covar_name, "_", level_num, ".txt"), sep = "\t", quote = F, col.names = NA)
    
    level_num <- level_num + 1
  }
} else {
  q1 <- quantile(covars[,covar_name], probs = 0.25)
  q3 <- quantile(covars[,covar_name], probs = 0.75)
  covars_q1 <- covars[covars[,covar_name] < q1,]
  covars_q3 <- covars[covars[,covar_name] > q3,]
  
  #Remove covariates with near zero variance
  near_zero_var_q1 <- nearZeroVar(covars_q1,  uniqueCut = 5)
  covars_q1 <- covars_q1[, -near_zero_var_q1]
  near_zero_var_q3 <- nearZeroVar(covars_q3,  uniqueCut = 5)
  covars_q3 <- covars_q3[, -near_zero_var_q3]

  expr_summary_q1 <- exp_summary(norm_expression[,row.names(covars_q1) ])
  expr_summary_q3 <- exp_summary(norm_expression[,row.names(covars_q3) ])
  
  cat ("level_num 1 corresponds to the first quartile, level_num 2 - to the third quartile \n")

  write.table(covars_q1, file = paste0(output_prefix, "covariates.", covar_name, "_1.txt"), sep = "\t", quote = F, col.names = NA)
  write.table(covars_q3, file = paste0(output_prefix, "covariates.", covar_name, "_2.txt"), sep = "\t", quote = F, col.names = NA)

  write.table(expr_summary_q1, file = paste0(output_prefix, "expression_summary.", covar_name, "_1.txt"), sep = "\t", quote = F, col.names = NA)
  write.table(expr_summary_q3, file = paste0(output_prefix, "expression_summary.", covar_name, "_2.txt"), sep = "\t", quote = F, col.names = NA)
}
                                                                                                                                                                                                                   
