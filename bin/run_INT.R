args <- commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]

covar <- read.delim(in_file, check.names = F, header = T, row.names = 1, as.is = T, sep ="\t")

contin_covars <- apply(covar, 2, function(x) length(unique(x)) > 3)
covar_int <- cbind(covar[,!contin_covars], apply(covar[,contin_covars], 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) ))
colnames(covar_int) <- c(colnames(covar)[!contin_covars], colnames(covar)[contin_covars])

write.table(covar_int, file = out_file, sep = "\t", quote = F, col.names = NA)
