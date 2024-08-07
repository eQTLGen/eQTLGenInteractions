args <- commandArgs(trailingOnly = TRUE)

deconvolve_nnls <- function(y, x){
  library(nnls)
  fit <- apply(y, 2, function(z) {
    nnls(A = x, b = as.vector(z))})
  
  coeff <- data.frame(lapply(fit, function(z) z$x))
  colnames(coeff) <- colnames(y)
  rownames(coeff) <- colnames(x)
  coeff <- t(coeff)
  
  return(coeff)
}

deconvolve_dtangle <- function(y, x, dtype = "rna-seq"){
  library(dtangle)
  res <- dtangle(
    Y = log2(t(y)+1), 
    references = log2(t(x)+1),
    marker_method = "ratio",
    data_type = dtype)
  
  return(res$estimates)
}

expr_fname = args[1]
sign_fname = args[2]
method= args[3]
out_fname = args[4]

datatype="rna-seq"
if (length(args) > 4) datatype = args[5]

if (datatype %in% c("rna-seq", "rnaseq", "RNAseq", "RNAseq_HGNC")) datatype="rna-seq"
if (datatype %in% c("microarray", "array", "HT12v3", "HT12v4", "HuRef8", "AffyU219", "AffyHumanExon")) datatype = "microarray-gene"


bulk_expr_tpm <- as.matrix(read.delim(expr_fname, row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
signature_matrix <- as.matrix(read.delim(sign_fname, row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))

shared_genes <- intersect(row.names(bulk_expr_tpm), row.names(signature_matrix))
cat("Running deconvolution using method:", method, "; signature matrix:", sign_fname, "\n")
cat(length(shared_genes), "genes shared between expression table and signature matrix\n")

if (method == "nnls"){
  coef <- deconvolve_nnls(x = signature_matrix[shared_genes,], y = bulk_expr_tpm[shared_genes,])
} else if (method == "dtangle") {
  coef <- deconvolve_dtangle(x = signature_matrix[shared_genes,], y = bulk_expr_tpm[shared_genes,], dtype = datatype)
} else {
  cat ("Wrong deconvolution method!\n")
}
proportions <- round((100/max(rowSums(coef), na.rm=TRUE))*coef, 4)

# remove zero variance cell types
mode(proportions) <- "numeric"
proportions <- data.frame(proportions)
proportions <- proportions[,which(!unlist(lapply(proportions, function(x) 0 == var(x)))) ]

write.table(proportions, file = out_fname, sep = "\t", quote = F, col.names = NA)

