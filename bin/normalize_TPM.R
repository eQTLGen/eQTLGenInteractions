args <- commandArgs(trailingOnly = TRUE)

bulk_expr <- as.matrix(read.delim(args[1], row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
gene_lengths <- read.delim(args[2], row.names = 1, header = F, as.is = T, check.names = F, sep = "\t")
outfile=args[3]

# TPM normalize
bulk_expr <- bulk_expr[which(rowSums(bulk_expr) != 0),]
shared_genes <- intersect(row.names(bulk_expr), row.names(gene_lengths))
cat("Number of genes in expr matrix: ", nrow(bulk_expr), "\nNumber of genes shared with gene length table: ", length(shared_genes), "\n")
gene_lengths <- gene_lengths[shared_genes,]
bulk_expr <- bulk_expr[shared_genes,]

bulk_expr_tpm <- do.call(cbind, lapply(1:ncol(bulk_expr), function(i) {
  rate = log(bulk_expr[,i]) - log(gene_lengths)
  denom = log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}))
colnames(bulk_expr_tpm) <- colnames(bulk_expr)


write.table(na.omit(bulk_expr_tpm), file = outfile, sep = "\t", quote = F, col.names = NA)



