### DESeq2 Pipeline - Andrew R Gross - 2018-08-06
### A standard exectution of the DESeq2 pipeline.  
### INPUT: This script requires a counts table.  Normalized TPM data and sample data is recommended.
### OUTPUT: This script generates a table normalized expression with fold change and p-values between sample groups.

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
  descr <- getBM(attributes=c('ensembl_gene_id','transcript_biotype','name_1006','namespace_1003','phenotype_description'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  descr <- descr[match(row.names(dataframe),descr[,1]),]
  #for (rowNumber in 1:length(descr[,1])) {
  #  newDescr <- descr[rowNumber,][,c(2,3,4,5)]
  #  descriptions <- c(descriptions, newDescr)
  #}
  dataframe <- cbind(dataframe, descr[c(2,3,4,5)])
  return(dataframe)
}
#getBM(attributes=c('transcript_biotype','name_1006','namespace_1003','phenotype_description'), filters='ensembl_gene_id', values= id, mart=ensembl)
#dataframe <- head(results.cov)
#dataframe <- data.frame(dataframe, test = '', test2 = '')
####################################################################################################################################################
### Input
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/')
#counts.cov <- read.csv('VW-10228--08--04--2020_COUNTS.csv', row.names = 1)
counts.cov <- read.csv('Counts wo lncRNA transcripts - all.csv', row.names = 1)


#sample.names <- c('mock-d1-r1', 'mock-d1-r2', 'mock-d1-r3', 'cov-d1-r1', 'cov-d1-r2', 'cov-d1-r3', 'mock-d3-r1', 'mock-d3-r2', 'mock-d3-r3', 'cov-d3-r1', 'cov-d3-r2', 'cov-d3-r3')

genes.of.interest <- read.table('Genes of interest Mock vs Infected.txt', sep = '\t', header = TRUE)[,1]

####################################################################################################################################################
### Format
##########################################################################
### Reassign names
#names(counts.cov) <- sample.names

### Filter low expression genes
summary(counts.cov)
counts.max <- apply(counts.cov, 1, median)
quantile(counts.max, seq(0,1,0.1))
rows.to.keep <- which(counts.max >5)
length(rows.to.keep)

### Output reformatted data for other programs
#counts.cov.out <- convert.ids(counts.cov)
#setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/Analysis/Input files/')
#write.csv(counts.cov.out, 'counts - reformatted.csv')

############################################################################################
### Filter data

### Filter rows
results <- counts.cov[rows.to.keep,]

### Select Columns
results <- results ; subtitle = 'All samples'
results.d1 <- results[c(1,2,3,4,5,6)] ;subtitle = 'Day 1 only'
results.d3 <- results[c(7,8,9,10,11,12)] ;subtitle = 'Day 3 only'

##########################################################################
### Reorder columns

### Make a column data data frame
columndata <- data.frame(row.names = names(results), disease = factor(c(rep('mock',3),rep('cov',3),rep('mock',3),rep('cov',3)), levels = c('mock','cov')))
columndata.d1 <- data.frame(row.names = names(results.d1), disease = factor(c(rep('mock',3),rep('cov',3)), levels =  c('mock','cov')))
columndata.d3 <- data.frame(row.names = names(results.d3), disease = factor(c(rep('mock',3),rep('cov',3)), levels =  c('mock','cov')))

# Convert counts into a matrix
results.m <- as.matrix(results)
results.d1.m <- as.matrix(results.d1)
results.d3.m <- as.matrix(results.d3)

####################################################################################################################################################
### Differential Expression
### Make our DESeq data sets
dds.cov <- DESeqDataSetFromMatrix(countData = results, colData = columndata, design = ~ disease)
dds.cov1 <- DESeqDataSetFromMatrix(countData = results.d1, colData = columndata.d1, design = ~ disease)
dds.cov3 <- DESeqDataSetFromMatrix(countData = results.d3, colData = columndata.d3, design = ~ disease)

### Run DESeq
dds.cov <- DESeq(dds.cov)
dds.cov1 <- DESeq(dds.cov1)
dds.cov3 <- DESeq(dds.cov3)

####################################################################################################################################################
### Format Results
results.cov <- as.data.frame(results(dds.cov))
results.cov1 <- as.data.frame(results(dds.cov1))
results.cov3 <- as.data.frame(results(dds.cov3))

### Drop NAs
results.cov <- results.cov[!is.na(results.cov$padj),]
results.cov1 <- results.cov1[!is.na(results.cov1$padj),]
results.cov3 <- results.cov3[!is.na(results.cov3$padj),]

##########################################################################
### Filter results
### Round figures

### Keep p-adjusted
results.cov <- round.DESeq.results(results.cov)[-c(3,4,5)]
results.cov1 <- round.DESeq.results(results.cov1)[-c(3,4,5)]
results.cov3 <- round.DESeq.results(results.cov3)[-c(3,4,5)]

##########################################################################
### Annotate genes
#results.cov <- convert.ids(results.cov) ; nrow(results.cov)
#results.cov1 <- convert.ids(results.cov1) ; nrow(results.cov1)
#results.cov3 <- convert.ids(results.cov3) ; nrow(results.cov3)

### Filter out LINCs
#results.cov <- results.cov[-grep('LINC', results.cov$Gene),] ; nrow(results.cov)
#results.cov1 <- results.cov1[-grep('LINC', results.cov1$Gene),] ; nrow(results.cov1)
#results.cov3 <- results.cov3[-grep('LINC', results.cov3$Gene),] ; nrow(results.cov3)

##########################################################################
### Volcano Plots

### Format #############################################################################################
volcano.data.cov <- results.cov3
volcano.data.cov$Gene = rownames(volcano.data.cov)
subtitle = 'Day 3 Infected vs. Mock'

### Converting FC from log2 to log10
volcano.data.cov$log10FC <- log10(2^volcano.data.cov$log2FoldChange)

volcano.data.cov.sig <- volcano.data.cov[volcano.data.cov$padj<0.05,]
volcano.data.cov.sig <- volcano.data.cov.sig[abs(volcano.data.cov.sig$log2FoldChange)>0.5,]
volcano.data.cov.sig <- volcano.data.cov.sig[order(volcano.data.cov.sig$baseMean),]


### Selecting genes to list
volcano.data.cov.sig$Edge.up <- -log(volcano.data.cov.sig$padj) * volcano.data.cov.sig$log2FoldChange
volcano.data.cov.sig$Edge.down <- -log(volcano.data.cov.sig$padj) * -volcano.data.cov.sig$log2FoldChange

volcano.text.up <- volcano.data.cov.sig[volcano.data.cov.sig$Edge.up > 12,]; nrow(volcano.text.up)
volcano.text.down <- volcano.data.cov.sig[volcano.data.cov.sig$Edge.down > 20,]; nrow(volcano.text.down)

### Manually filtering
#genes.to.omit <- c('ANKRD36', 'AC022400.7', 'LINC00342' )
#volcano.text.up <- volcano.text.up[!volcano.text.up$Gene %in% genes.to.omit,]
#volcano.text.down <- volcano.text.down[!volcano.text.down$Gene %in% genes.to.omit,]



ggplot( ) +
  geom_point(data = volcano.data.cov, aes(x=log2FoldChange, y = log10(padj) ),color = 'grey', size = 0.9) +
  geom_point(data = volcano.data.cov.sig, aes(x=log2FoldChange, y = log10(padj), color = log10(baseMean), size = log10(baseMean))) +
  geom_text(data = volcano.text.up, aes(x=log2FoldChange + 0.2, y = log10(padj) - 0, label = Gene), 
            hjust = 0, vjust = 0.1, size = 3, check_overlap = TRUE) +
  geom_text(data = volcano.text.down, aes(x=log2FoldChange - 0.2, y = log10(padj) - 0, label = Gene), 
            hjust = 1, vjust = 0.1, size = 3, check_overlap = TRUE) +
  scale_color_gradient(low="pink", high="red") +
  scale_size('Log10 Expression', range = c(0.5,4)) +
  ylim(c(0, min(log10(volcano.data.cov$padj)))) +
  xlim(c(min(volcano.data.cov$log2FoldChange-0.5), max(volcano.data.cov$log2FoldChange)+1.5)) +
  labs(title="P-Value vs. Log Fold Change Volcano Plot",
       subtitle = subtitle,
       x = 'Log2 Fold Change', 
       y = 'Log10 P-value (Adjusted)') +
  theme(plot.title = element_text(hjust = 0.5, color="black", face="bold", size=18, margin=margin(0,0,5,0)),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_text(face="bold", size=12,margin =margin(5,0,5,0)),
        axis.title.y = element_text(face="bold", size=12,margin =margin(0,5,0,5)),
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12),
        legend.position = 'none') 

      

volcano.text.down 
volcano.text.up

