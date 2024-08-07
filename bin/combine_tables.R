args <- commandArgs(trailingOnly = TRUE)
library(data.table)

table_list <- list()
cnt <- 1
for (path in args[2:length(args)]){
    cat("Reading", path, "\n")
    table_list[[cnt]] <- fread(path, data.table = F)
    colnames(table_list[[cnt]])[1] <- "sample"
    cnt <- cnt + 1
}

merged <- na.omit(Reduce(function(...) merge(..., all = TRUE, by = "sample"), table_list))

cat ("Number of covariates in resulting file:", ncol(merged), "\n")
cat("\tCovariate names:", paste(colnames(merged), collapse = ","), "\n")
cat("Number of samples with data available:", nrow(merged), "\n")

write.table(merged, file = args[1], sep = "\t", quote = F, row.names = F)
