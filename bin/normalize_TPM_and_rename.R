args <- commandArgs(trailingOnly = TRUE)

bulk_expr <- as.matrix(read.delim(args[1], row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
gene_lengths <- read.delim(args[2], row.names = 1, header = F, as.is = T, check.names = F, sep = "\t")
probe_annotation <- read.delim(args[3], row.names = 1, header = T, as.is = T, check.names = F, sep = "\t")
outfile=args[4]

rename = T
if (length(args) > 4){
  if (args[5] == "RNAseq_HGNC") rename = F
}

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

if (rename){
  # Convert Ensembl ids to gene names
  new_names <-  probe_annotation[match(row.names(bulk_expr_tpm), row.names(probe_annotation), nomatch = "NA"),"gene_name"]
  new_names <- gsub(" ", "", new_names, fixed = TRUE)
  row.names(bulk_expr_tpm) <- new_names

  n_occur <- data.frame(table(new_names))
  dups <- n_occur[n_occur$Freq > 1,"new_names"]
  bulk_expr_tpm <- bulk_expr_tpm[!row.names(bulk_expr_tpm) %in% dups,]
}

write.table(na.omit(bulk_expr_tpm), file = outfile, sep = "\t", quote = F, col.names = NA)



