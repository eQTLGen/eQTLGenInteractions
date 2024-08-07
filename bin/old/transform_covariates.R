library(optparse)
library(ggplot2)
library(caret)

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Path to the covariate file with sex and age. Sample ids = genotype ids."),
  make_option(c("-g", "--gte"), type = "character",
              help = "Path to the genotype to expression id conversion file."),
  make_option(c("-t", "--int"),  action="store_true", default=FALSE, type = "logical",
              help = "INT transform the table?"),
  make_option(c("-o", "--out"), type = "character",
              help = "Output file name."),
  make_option(c("-r", "--rename"), action="store_true", default=FALSE, type = "logical",
              help = "Are raw expression sample ids different from genotype sample ids?."),
  make_option(c("-c", "--cut"), type = "integer", default = NA,
              help = "Number of columns to keep")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)


rename_samples <- function(d, gte, key_col = 2, value_col = 1){
  if (length(intersect(row.names(d), gte[,key_col])) < length(intersect(row.names(d), gte[,value_col]))){
    cat("Sample ids are already converted to expression ids! Skipping id conversion.\n")
    return(d)
  }
  d_m <- d[match(gte[,key_col], row.names(d), nomatch = 0),, drop = F]
  gte_m <- gte[match(row.names(d_m), gte[,key_col], nomatch = 0),]
  row.names(d_m) <- gte_m[,value_col]
  cat("Coverted genotype sample ids to expression sample ids.\nN samples before: ", nrow(d), "\nN samples after: ", nrow(d_m), "\n")
  return(d_m)
}

run_INT <- function(d){
  contin_covars <- apply(d, 2, function(x) length(unique(x)) > 3)
  d_int <- cbind(d[,!contin_covars, drop = F], apply(d[,contin_covars, drop = F], 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) ))
  colnames(d_int) <- c(colnames(d)[!contin_covars], colnames(d)[contin_covars])
  return(d_int)
}

in_data <- as.data.frame(read.delim(args$input, check.names = F, header = T, row.names = 1, as.is = T, sep ="\t"))
gte <- read.delim(args$gte, check.names = F, header = F, as.is = T, sep ="\t")

cat("Number of covariates in the input file:", ncol(in_data), "\n")
cat("Number of samples in the input file:", nrow(in_data), "\n") 

if (args$rename){
  in_data <- rename_samples(in_data, gte)
}

# keep only samples present in gte
in_data <- in_data[match(gte[,1], row.names(in_data), nomatch = 0),, drop = F]

# keep the required number of columns
if (!(is.na(args$cut) | args$cut == "NA" | args$cut == 0)) in_data <- in_data[,1:args$cut]

# remove colums with all zeros and rows with NAs
if (ncol(in_data) > 1) in_data <- in_data[,colSums(is.na(in_data)) < nrow(in_data), drop = F]
in_data <- na.omit(in_data)

# remove columns with near zero variance
near_zero_var <- nearZeroVar(in_data,  uniqueCut = 5)
if (length(near_zero_var) > 0) in_data <- in_data[, -near_zero_var, drop = F]

# Plot covariate distributions prior to INT
pdf(paste0(args$out, ".distributions.pdf"), height = 5, width = 5)
for (covar_name in colnames(in_data)){
   print(ggplot(in_data, aes(get(covar_name))) + geom_histogram(bins = 50, color = "black", alpha = 0.5, fill = "dodgerblue4") + theme_bw() + xlab(covar_name) + ggtitle(paste0("Raw ", covar_name)))
}
dev.off()

if (args$int){
  in_data <- run_INT(in_data)
  # Plot covariate distributions after INT
    pdf(paste0(args$out, ".distributions.INT.pdf"), height = 5, width = 5)
    for (covar_name in colnames(in_data)){
        print(ggplot(in_data, aes(get(covar_name))) + geom_histogram(bins = 50, color = "black", alpha = 0.5, fill = "dodgerblue4") + theme_bw() + xlab(covar_name) + ggtitle(paste0("Raw ", covar_name)))
    }
    dev.off()
}

cat ("Number of covariates in resulting file:", ncol(in_data), "\n")
cat("Number of samples with data available:", nrow(in_data), "\n")

write.table(in_data, file = args$out, sep = "\t", quote = F, col.names = NA)




