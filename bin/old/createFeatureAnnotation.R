#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--in_gtf"), action="store", default=NA, type='character',
              help="Input GTF file. (required)"),
  make_option(c("--n_genes"), action="store", default=100, type='integer',
              help="Number of genes to test in one QTL map job. "),
  make_option(c("--feature_name"), action="store", default="gene_id",
              help="The feature_id used in the expression matrix."),
  make_option(c("--biotype_flag"), action="store", default="gene_biotype",
              help="The feature_id used in the expression matrix."),
  make_option(c("--autosomes_only"), action="store_true", default=FALSE, type='character',
              help="Flag to only keep autosomal genes, default false. "),
  make_option(c("--out_dir"), action="store", default=NA, type='character',
              help="Output main directory. (required)"))
opt = parse_args(OptionParser(option_list=option_list))

# check date is provided
if ((is.na(opt$in_gtf) || is.na(opt$out_dir))) {
  stop("required parameters (in_gtf & out_dir) must be provided. ")
}

#### Read gene annotation
gtfInfo = read.delim(opt$in_gtf,as.is=T,comment="#",header=F)

##S1. make annotation file.

##Select only genes.
gtfInfo = gtfInfo[which(gtfInfo$V3=="gene"),]
##Make sure genes are mapped to numeric
gtfInfo$V1 = gsub("chr","",gtfInfo$V1)
##Select only main chromosome mappings.
gtfInfo = gtfInfo[which(gtfInfo$V1 %in% c(1:22,"X","MT", "M","Y")),]

##Extend matrix with relevant entries.
partialMatrix = strsplit(gtfInfo$V9,split = "; ")
ensgIdLoc = which(startsWith(partialMatrix[[1]],"gene_id"))
gtfInfo["ENSG"] = gsub( "gene_id ","",unlist(lapply(partialMatrix,"[[",ensgIdLoc)))
nameLoc = which(startsWith(partialMatrix[[1]],"gene_name"))
gtfInfo["Gene_Name"] = gsub( "gene_name ","",unlist(lapply(partialMatrix,"[[",nameLoc)))

# If there is no gene_name field, put gene id instead
gtfInfo[startsWith(gtfInfo$Gene_Name, "gene_source"), "Gene_Name"] <- gtfInfo[startsWith(gtfInfo$Gene_Name, "gene_source"), "ENSG"]


# If not working use:
# library(rtracklayer)
# gtfInfo <- readGFF(opt$in_gtf)
#

if(opt$biotype_flag!="NA"){
	btLoc = which(startsWith(partialMatrix[[1]],opt$biotype_flag))
	gtfInfo["biotype"] = gsub(paste(opt$biotype_flag," ",sep=""),"",unlist(lapply(partialMatrix,"[[",btLoc)))
} else {
	gtfInfo["biotype"] = NA
}
print(head(gtfInfo))

if(opt$feature_name=="ENSG"){
  geneInfo = gtfInfo[,c(10,1,4,5,11,12)]
  colnames(geneInfo)=c("feature_id","chromosome","start","end","GeneName","biotype")
} else {
  geneInfo = gtfInfo[,c(11,1,4,5,10,12)]
  colnames(geneInfo)=c("feature_id","chromosome","start","end","ENSG","biotype")
}

if(opt$autosomes_only){
  geneInfo = geneInfo[which(geneInfo$chromosome %in% c(1:22)),]
}

##drop all genes that have the same names, to avoid having issues with matching and mapping after.
toDrop = geneInfo$feature_id[which(duplicated(geneInfo$feature_id))]
geneInfo = geneInfo[which(!(geneInfo$feature_id %in% toDrop)),]

write.table(geneInfo,paste0(opt$out_dir,"/LimixAnnotationFile.txt"),quote = F,sep="\t",row.names = F)

##S2. make chunk file.
testCombinations = NULL
#
nGenes = opt$n_genes
startPos = 0
endOffset = 1000000000

lines = NULL
for(chr in unique(geneInfo$chromosome)){
  #print(chr)
  annotationRel = geneInfo[which(geneInfo$chromosome==chr),]
  annotationRel = annotationRel[order(annotationRel$start,annotationRel$end),]
  ##First go through the list to fix genes to they all touch.
  annotationRel$start[1] = startPos
  for(i in 2:nrow(annotationRel)){
    if(i == nrow(annotationRel)){
      annotationRel$end[i] = annotationRel$end[i]+endOffset
    }
    #If "overlapping" than we don't need to do anything.
    if((annotationRel$start[i]>annotationRel$end[i-1])){
      #print(i)
      distance = (annotationRel$start[i]-annotationRel$end[i-1])/2
      annotationRel$start[i] = ceiling(annotationRel$start[i]-distance)
      annotationRel$end[i-1] = floor(annotationRel$end[i-1]+distance)
    }
  }
  chunks = seq(1, nrow(annotationRel),nGenes)
  #Need to add the last one as a separate entry.
  if(chunks[length(chunks)] < nrow(annotationRel)){
    chunks = c(chunks,(nrow(annotationRel)+1))
  }
  for(i in 1:(length(chunks)-1)){
    lines = c(lines,(paste(chr,":",annotationRel$start[chunks[i]],"-",annotationRel$end[(chunks[i+1]-1)],sep="")))
  }
}
write.table(lines,paste0(opt$out_dir,"/ChunkingFile.txt"),quote = F,sep="\t",row.names = F,col.names=F)


##S3. make gene lengths file
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(opt$in_gtf,format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))

write.table(exonic.gene.sizes,paste0(opt$out_dir,"/GeneLengths_ensg.txt"),quote = F,sep="\t",col.names=F)
