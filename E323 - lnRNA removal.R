### Long Non-coding RNA removal - Andrew R Gross - 2020-09-08
### A standard exectution of the DESeq2 pipeline.  
### INPUT: A table of transcript counts
### OUTPUT: A new table of transcript counts without lnRNA transcripts.

####################################################################################################################################################
### Header
library("pasilla")
library("DESeq2")
library("biomaRt")
#library("heatmap2")
#library("VennDiagram")
library (ggplot2)


ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
listMarts(host="www.ensembl.org")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

round.DESeq.results <- function(dataframe) {
  dataframe$baseMean <- round(dataframe$baseMean, 2)
  dataframe$log2FoldChange <- round(dataframe$log2FoldChange, 2)
  dataframe$lfcSE <- round(dataframe$lfcSE, 3)
  dataframe$stat <- round(dataframe$stat, 2)
  #dataframe$pvalue <- formatC(dataframe$pvalue, format = "e", digits = 2)
  #dataframe$padj <- formatC(dataframe$padj, format = "e", digits = 2)
  return(dataframe)
}
convert.ids <- function(dataframe, add.gene.name.column = TRUE) {
  ### This function will convert a row name consisting of a contactenated ensembl ID and gene to one or the other,
  ### based on the users instruction (2018-10-04)
  ensemblIDs <- c()                                           # Empty lists are initialized to receive IDs as they're created
  gene.names <- c()
  for (rowName in row.names(dataframe)) {                     # Loops through all rows in the data frame
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]                 # Splits the row name and declares the ensembl ID
    gene.name <- strsplit(rowName,"\\_")[[1]][2]                 # Splits the row name, declares the gene name
    ensemblIDs <- c(ensemblIDs, ensemblID)                       # Adds ensembl ID and gene name to appropriate lists
    gene.names <- c(gene.names, gene.name)
  }
  row.names(dataframe) <- make.unique(ensemblIDs)                          # assigns the new row names
  if(add.gene.name.column == TRUE) {
    dataframe$Gene <- gene.names
  }
  return(dataframe)                                           # Returns the data frame with new rows
}

addGene <- function(dataframe) {
  genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  genes <- genes[match(row.names(dataframe),genes[,1]),]
  Gene <- c()
  for (rowNumber in 1:length(genes[,1])) {
    newGene <- genes[rowNumber,][,2]
    Gene <- c(Gene, newGene)
  }
  dataframe[length(dataframe)+1] <- Gene
  names(dataframe)[ncol(dataframe)] <- "Gene"
  return(dataframe)
}
add.description <- function(dataframe) {
  descr <- getBM(attributes=c('ensembl_gene_id','description'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  descr <- descr[match(row.names(dataframe),descr[,1]),]
  descriptions <- c()
  for (rowNumber in 1:length(descr[,1])) {
    newDescr <- descr[rowNumber,][,2]
    newDescr <- strsplit(newDescr, " \\[")[[1]][1]
    descriptions <- c(descriptions, newDescr)
  }
  dataframe[length(dataframe)+1] <- descriptions
  names(dataframe)[ncol(dataframe)] <- "Description"
  return(dataframe)
}

add.type <- function(dataframe) {
  descr <- getBM(attributes=c('ensembl_gene_id','transcript_biotype'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  descr <- descr[match(row.names(dataframe),descr[,1]),]
  dataframe <- cbind(dataframe, descr[c(2,3,4,5)])
  return(dataframe)
}

getBM(attributes=c('ensembl_gene_id','transcript_biotype'), filters='ensembl_gene_id', values= row.names(dataframe), mart=ensembl)
dataframe <- head(results.cov)
#dataframe <- data.frame(dataframe, test = '', test2 = '')
####################################################################################################################################################
### Input
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/')
counts.cov <- read.csv('VW-10228--08--04--2020_COUNTS.csv', row.names = 1)

sample.names <- c('mock-d1-r1', 'mock-d1-r2', 'mock-d1-r3', 
                  'cov-d1-r1', 'cov-d1-r2', 'cov-d1-r3', 
                  'mock-d3-r1', 'mock-d3-r2', 'mock-d3-r3', 
                  'cov-d3-r1', 'cov-d3-r2', 'cov-d3-r3')

####################################################################################################################################################
### Format
##########################################################################
### Reassign names
names(counts.cov) <- sample.names

### Filter low expression genes
summary(counts.cov)
counts.max <- apply(counts.cov, 1, max)
quantile(counts.max, seq(0,1,0.1))
rows.to.keep <- which(counts.max >5)
length(rows.to.keep)

############################################################################################
### Filter data

### Filter rows
results.cov <- counts.cov[rows.to.keep,]

############################################################################################
### Annotate rows
############################################################################################

### Annotate genes
results.cov <- convert.ids(results.cov)

##########################################################################
### Filter out Long non-coding transcripts
transcript.type.table <- getBM(attributes=c('ensembl_gene_id','transcript_biotype'), filters='ensembl_gene_id', values= row.names(results.cov), mart=ensembl)
transcript.type.table <- transcript.type.table
transcript.type.table$transcript_biotype <- as.factor(transcript.type.table$transcript_biotype)

summary(transcript.type.table)

lnrna.transcripts <- grep("lncRNA", transcript.type.table$transcript_biotype)
lnrna.table <- transcript.type.table[lnrna.transcripts,]

length(row.names(results.cov))
length(lnrna.table$ensembl_gene_id)

#test3 <- results.cov[lnrna.table$ensembl_gene_id,]
lnrna.rows <- match(lnrna.table$ensembl_gene_id, row.names(results.cov))
results.wo.lnrna <- results.cov[-lnrna.rows,]

##########################################################################
### Rename rows with genes
### Find the redundant genes
unique.genes <- unique(results.wo.lnrna$Gene)
unique.genes.pos <- match(unique.genes, results.wo.lnrna$Gene)
redundant.gene.list <- unique(results.wo.lnrna$Gene[-unique.genes.pos])
test <- redundant.gene.list %in% results.wo.lnrna$Gene
test2 <- results.wo.lnrna$Gene %in% redundant.gene.list 
redundant.gene.table <- results.wo.lnrna[test2,]

test3 <- match(redundant.gene.list, results.wo.lnrna$Gene)
test3 <- which(results.wo.lnrna$Gene)
test3 <- intersect(redundant.gene.list, results.wo.lnrna$Gene)
#redundant.gene.table[order(redundant.gene.table$`mock-d1-r1`, decreasing = TRUE),]

new.summed.rows <- redundant.gene.table[0,]


for(counter in 1:length(redundant.gene.list)){
  current.gene <- redundant.gene.list[counter]
#grep(redundant.gene.list[1], redundant.gene.table$Gene)
  rows.to.sum <- redundant.gene.table[grep(current.gene, redundant.gene.table$Gene),]
  rows.to.sum$Gene=1
  new.row <- apply(rows.to.sum,2,sum)
  new.summed.rows[counter,] <- new.row
  new.summed.rows$Gene[counter] = current.gene

}
new.summed.rows

### Remove the redundant genes, then add back in their sums

unique.genes.minus.red <- setdiff(unique.genes,redundant.gene.list)
results.wo.lnrna.wo.red <- results.wo.lnrna[match(unique.genes.minus.red, results.wo.lnrna$Gene),]
results.wo.lnrna.wo.red <- rbind(results.wo.lnrna.wo.red, new.summed.rows)
row.names(results.wo.lnrna.wo.red) <- results.wo.lnrna.wo.red$Gene

results.wo.lnrna <- results.wo.lnrna.wo.red[1:12]

#results.wo.lnrna <- rbind(results.wo.lnrna, new.summed.rows)
#nrow(results.wo.lnrna)
#length(results.wo.lnrna.wo.red$Gene)
#length(unique(results.wo.lnrna.wo.red$Gene))

#results.wo.lnrna <- results.wo.lnrna[1:12]
##########################################################################
### Output the data

getwd()

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/Analysis/Input files formatted for analysis/')

write.csv(results.wo.lnrna, 'Counts wo lncRNA transcripts - all.csv')
write.csv(results.wo.lnrna[c(1,2,3,4,5,6)], 'Counts wo lncRNA transcripts - day 1.csv')
write.csv(results.wo.lnrna[c(7,8,9,10,11,12)], 'Counts wo lncRNA transcripts - day 3.csv')

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/Analysis/ALTERNATIVE ANALYSIS - R3 sample reversed/Input Files formatted for analysis - R3 switched/')
results.wo.lnrna.B <- results.wo.lnrna[c(1,2,9,4,5,6,7,8,4,10,11,12)]
names(results.wo.lnrna.B) <- names(results.wo.lnrna)
write.csv(results.wo.lnrna.B, 'Counts wo lncRNA transcripts - all - ALT.csv')
write.csv(results.wo.lnrna.B[c(1,2,3,4,5,6)], 'Counts wo lncRNA transcripts - day 1 - ALT.csv')
write.csv(results.wo.lnrna.B[c(7,8,9,10,11,12)], 'Counts wo lncRNA transcripts - day 3- ALT.csv')
