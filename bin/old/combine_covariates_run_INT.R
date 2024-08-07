args <- commandArgs(trailingOnly = TRUE)

library(dplyr)

in_file <- args[1]
in_file2 <- args[2]
out_file <- args[3]

covar1 <- read.delim(in_file, check.names = F, header = T, row.names = 1, as.is = T, sep ="\t")
covar2 <- read.delim(in_file2, check.names = F, header = T, row.names = 1, as.is = T, sep ="\t")

covar <- merge(covar1, covar2, by = 0)
row.names(covar) <- covar$Row.names
covar$Row.names = NULL

covar <- covar[,colSums(is.na(covar)) < nrow(covar)]
covar <- na.omit(covar)

contin_covars <- apply(covar, 2, function(x) length(unique(x)) > 3)
covar_int <- cbind(covar[,!contin_covars], apply(covar[,contin_covars], 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) ))
colnames(covar_int) <- c(colnames(covar)[!contin_covars], colnames(covar)[contin_covars])

cat ("Number of covariates in the combined file:", ncol(covar), "\n")
cat("Number of samples with covariate data available:", nrow(covar), "\n")

write.table(covar_int, file = out_file, sep = "\t", quote = F, col.names = NA)

