#!/usr/bin/Rscript
args <- commandArgs(TRUE)

require(edgeR)

countsFile = args[1]
normOut = args[2]
samplePerc = args[3] # (0.01)


table<-read.delim(gzfile(countsFile),row.names=1, check.names=FALSE)

# Remove genes without any counts
table <- table[rowMeans(table) > 0, ]

D<-DGEList(counts=table)
cpm.D <-cpm(D)

samplePerc <- as.numeric(samplePerc)
keep <- rowSums(cpm.D > 0.5) >= ncol(table)*samplePerc
write.table(keep, file="keep_genes.txt", sep=",", row.names=TRUE)

D <- D[keep, ,keep.lib.sizes=FALSE]
d <- calcNormFactors(D)
scalar <- d$samples$lib.size*d$samples$norm.factors/exp(mean(log(d$samples$lib.size*d$samples$norm.factors)))
scal.mat <- outer(rep(1,nrow(d$counts)), scalar)
scaled.counts <- d$counts/scal.mat

write.table(scaled.counts, file = normOut, sep = "\t", row.names=TRUE, quote=FALSE, col.names = NA)