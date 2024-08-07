args <- commandArgs(trailingOnly = TRUE)
library(dplyr)

exp_path <- args[1]
covar_path <- args[2]
out_dir <- args[3]

pc_path <- NULL
if (length(args) > 3){
  pc_path <- args[4]
}

expr <- as.data.frame(t(data.table::fread(exp_path, data.table = F)))
colnames(expr) <- expr[1,]
expr <- mutate_all(expr[-1,], function(x) as.numeric(as.character(x)))

covars <- data.table::fread(covar_path, data.table = F)
colnames(covars)[1] <- "id"

if (! is.null(pc_path) ){
  pcs <- data.table::fread(pc_path, data.table = F)
  colnames(pcs)[1] <- "id"
  #pcs <- pcs[,1:25]
  colnames(pcs) <- gsub("PC", "expPC", colnames(pcs))
  covars <- inner_join(covars, pcs, by = c("id"))
  row.names(covars) <- covars[,1]
  covars[,1] <- NULL
} else {
  row.names(covars) <- covars$id
  covars[,"id"] <- NULL
}

covars <- na.omit(covars)

ids <- intersect(row.names(expr), row.names(covars))
expr <- expr[ids,]
covars <- covars[ids,]

new_expr <- data.frame(matrix(ncol = ncol(expr), nrow = nrow(expr)))
colnames(new_expr) <- colnames(expr)
row.names(new_expr) <- row.names(expr)

for (col in 1:ncol(expr)){
  tmp_data <- cbind(expr[,col], covars)
  colnames(tmp_data)[1] <- "gene"
  print(str(tmp_data))
  lm_fit <- lm(gene ~ . ,data = tmp_data)
  new_expr[,col] <- lm_fit$residuals
}

write.table(t(new_expr), file = paste0(out_dir, "/expression_corrected.noINT.txt"), sep = "\t", quote = F, col.names = NA)
 

