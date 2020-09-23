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

####################################################################################################################################################
### Input
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/')
counts.cov <- read.csv('VW-10228--08--04--2020_COUNTS.csv', row.names = 1)

sample.names <- c('mock-d1-r1', 'mock-d1-r2', 'mock-d1-r3', 
                  'cov-d1-r1', 'cov-d1-r2', 'cov-d1-r3', 
                  'mock-d3-r1', 'mock-d3-r2', 'mock-d3-r3', 
                  'cov-d3-r1', 'cov-d3-r2', 'cov-d3-r3')

#counts.cov <- counts.cov[c(1,2,9,4,5,6,7,8,3,10,11,12)] ; sample.names <- c('mock-d1-r1 ALT', 'mock-d1-r2 ALT', 'mock-d1-r3 ALT', 
#                  'cov-d1-r1 ALT', 'cov-d1-r2 ALT', 'cov-d1-r3 ALT', 
#                  'mock-d3-r1 ALT', 'mock-d3-r2 ALT', 'mock-d3-r3 ALT', 
#                  'cov-d3-r1 ALT', 'cov-d3-r2 ALT', 'cov-d3-r3 ALT')

genes.of.interest <- read.table('Genes of interest Mock vs Infected.txt', sep = '\t', header = TRUE)[,1]

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

head(results)

### Convert to matrix
#results <- as.matrix(results)

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
pvalue.cutoff <- 0.05

### Keep p-adjusted
results.cov <- round.DESeq.results(results.cov)[-c(3,4,5)]
results.cov1 <- round.DESeq.results(results.cov1)[-c(3,4,5)]
results.cov3 <- round.DESeq.results(results.cov3)[-c(3,4,5)]

### Keep raw p-value
#results.cov <- round.DESeq.results(results.cov)[-c(3,4,6)]

### Apply cutoff
results.cov <- results.cov[results.cov$padj <= pvalue.cutoff,]
results.cov1 <- results.cov1[results.cov1$padj <= pvalue.cutoff,]
results.cov3 <- results.cov3[results.cov3$padj <= pvalue.cutoff,]

#results.cov <- results.cov[results.cov$pvalue <= pvalue.cutoff,]

### Reorder by p-value
results.cov <- results.cov[order(results.cov$padj),]
results.cov1 <- results.cov1[order(results.cov1$padj),]
results.cov3 <- results.cov3[order(results.cov3$padj),]

#results.cov <- results.cov[order(results.cov$pvalue),]

head(results.cov3, 20)

### Reorder by Fold-change
#results.cov <- results.cov[order(results.cov$log2FoldChange),]
#results.ko <- results.ko[order(results.ko$log2FoldChange),]

##########################################################################
### Apply cutoffs
#summary(counts.cov)
#counts.cutoff = 1
#results.cov <- results.cov[results.cov$basemean >= counts.cutoff,]

##########################################################################
### Add additional columns for context
### Add in normalized expression values
counts.cov2 <- counts.cov[row.names(results.cov),][c(1,2,3,7,8,9,4,5,6,10,11,12)]
results.cov <- cbind(results.cov, counts.cov2)
counts.cov2 <- counts.cov[row.names(results.cov1),][c(1,2,3,7,8,9,4,5,6,10,11,12)]
results.cov1 <- cbind(results.cov1, counts.cov2)
counts.cov2 <- counts.cov[row.names(results.cov3),][c(1,2,3,7,8,9,4,5,6,10,11,12)]
results.cov3 <- cbind(results.cov3, counts.cov2)

### Annotate genes
results.cov <- convert.ids(results.cov)
#results.cov <- addGene(results.cov)
results.cov <- add.description(results.cov)
results.cov$GOI.check <- match(results.cov$Gene, genes.of.interest, nomatch = ' ')

results.cov1 <- convert.ids(results.cov1)
#results.cov1 <- addGene(results.cov1)
results.cov1 <- add.description(results.cov1)
results.cov1$GOI.check <- match(results.cov1$Gene, genes.of.interest, nomatch = '')


results.cov3 <- convert.ids(results.cov3)
#results.cov3 <- addGene(results.cov3)
results.cov3 <- add.description(results.cov3)
results.cov3$GOI.check <- match(results.cov3$Gene, genes.of.interest, nomatch = '')

head(results.cov3)

##########################################################################
### Generate Log-Fold Change tables in each direction

gen.lfc.df <- function(dataframe) {
  sorted.df <- dataframe[order(dataframe$log2FoldChange, decreasing = TRUE),]
  up <- sorted.df[sorted.df$log2FoldChange >0,]
  down <- sorted.df[sorted.df$log2FoldChange <0,]
  down <- down[seq(nrow(down),1),]
  return(list(up,down))
}

lfc.sort <- gen.lfc.df(results.cov)
results.cov.up <- lfc.sort[[1]]
results.cov.down <- lfc.sort[[2]]

lfc.sort <- gen.lfc.df(results.cov1)
results.cov1.up <- lfc.sort[[1]]
results.cov1.down <- lfc.sort[[2]]

lfc.sort <- gen.lfc.df(results.cov3)
results.cov3.up <- lfc.sort[[1]]
results.cov3.down <- lfc.sort[[2]]

### Add type
results.cov.up <- add.type(results.cov.up)
results.cov.down <- add.type(results.cov.down)

results.cov1.up <- add.type(results.cov1.up)
results.cov1.down <- add.type(results.cov1.down)

results.cov3.up <- add.type(results.cov3.up)
results.cov3.down <- add.type(results.cov3.down)

##########################################################################
### Filter out Long non-coding transcripts

rm.lnc <- function(dataframe){
  rows.to.remove <- which(dataframe$transcript_biotype == 'lncRNA')
  print(paste(length(rows.to.remove),'Rows removed'))
  return(dataframe[-rows.to.remove,])
}

results.cov.up <- rm.lnc(results.cov.up)
results.cov.down <- rm.lnc(results.cov.down)

results.cov1.up <- rm.lnc(results.cov1.up)
results.cov1.down <- rm.lnc(results.cov1.down)

results.cov3.up <- rm.lnc(results.cov3.up)
results.cov3.down <- rm.lnc(results.cov3.down)

##########################################################################
### Output the data

getwd()

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/Analysis/Differential Expression/')
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/Analysis/ALTERNATIVE ANALYSIS - R3 sample reversed/Differential Expression - R3 switched/')

write.csv(results.cov, 'DEGs - COV all 0.05 padj - ALT.csv')
write.csv(results.cov.up, 'DEGs - COV all 0.05 padj - UP - ALT.csv')
write.csv(results.cov.down, 'DEGs - COV all 0.05 padj - DOWN - ALT.csv')

write.csv(results.cov1, 'DEGs - COV Day 1 0.05 padj - ALT.csv')
write.csv(results.cov1.up, 'DEGs - COV Day 1 0.05 padj - UP - ALT.csv')
write.csv(results.cov1.down, 'DEGs - COV Day 1 0.05 padj - DOWN - ALT.csv')

write.csv(results.cov3, 'DEGs - COV Day 3 0.05 padj - ALT.csv')
write.csv(results.cov3.up, 'DEGs - COV Day 3 0.05 padj - UP - ALT.csv')
write.csv(results.cov3.down, 'DEGs - COV Day 3 0.05 padj - DOWN - ALT.csv')

results.cov3.down <- results.cov3.down[order(results.cov3.down$padj),]
results.cov3.up <- results.cov3.up[order(results.cov3.up$padj),]

results.cov3.500 <- head(results.cov3.up,500)

results.cov3.500 <- head(results.cov3.down,500)

write.csv(results.cov3.500, 'top 500 day 3 upreg.csv')
write.csv(results.cov3.500, 'top 500 day 3 downreg.csv')

##########################################################################
### Heatmaps

### Format #############################################################################################
### SWAPPED COLUMNS
### Create color palette
pal <- colorRampPalette(c('yellow','black','blue'))(100)

heatmap.data.all <- results.cov[order(results.cov$log2FoldChange),] ; nrow(heatmap.data.all)
heatmap.data.all <- heatmap.data.all[heatmap.data.all$baseMean > 100,] ; nrow(heatmap.data.all)
heatmap.data.all <- heatmap.data.all[heatmap.data.all$padj < 0.01,][c(4,5,6,7,8,9,10,11,12,13,14,15)]  ; nrow(expression.df)
names(heatmap.data.all) <- c('Mock D1 #1','Mock D1 #2','Mock D1 #3','Mock D3 #1','Mock D3 #2','Mock D3 #3','COV D1 #1','COV D1 #2','COV D1 #3','COV D3 #1','COV D3 #2','COV D3 #3')

heatmap.data.d1 <- results.cov1[order(results.cov1$log2FoldChange),] ; nrow(heatmap.data.d1)
heatmap.data.d1 <- heatmap.data.d1[heatmap.data.d1$baseMean > 100,] ; nrow(heatmap.data.d1)
heatmap.data.d1 <- heatmap.data.d1[heatmap.data.d1$padj < 0.01,][c(4,5,6,10,11,12)]  ; nrow(expression.df)
names(heatmap.data.d1) <- c('Mock D1 #1','Mock D1 #2','Mock D1 #3','COV D1 #1','COV D1 #2','COV D1 #3')

heatmap.data.d3 <- results.cov3[order(results.cov3$log2FoldChange),] ; nrow(heatmap.data.d3)
heatmap.data.d3 <- heatmap.data.d3[heatmap.data.d3$baseMean > 100,] ; nrow(heatmap.data.d3)
heatmap.data.d3 <- heatmap.data.d3[heatmap.data.d3$padj < 0.01,][c(7,8,9,13,14,15)]  ; nrow(expression.df)
names(heatmap.data.d3) <- c('Mock D3 #1','Mock D3 #2','Mock D3 #3','COV D3 #1','COV D3 #2','COV D3 #3')


### Plot Heatmap

heatmap(as.matrix(log(heatmap.data.all)), Rowv = NA, Colv = NA, col = pal, scale = 'row', labRow = NA, margins = c(8,1))

heatmap(as.matrix(log(heatmap.data.d1)), Rowv = NA, Colv = NA, col = pal, scale = 'row', labRow = NA, margins = c(8,1))

heatmap(as.matrix(log(heatmap.data.d3)), Rowv = NA, Colv = NA, col = pal, scale = 'row', labRow = NA, margins = c(8,1))




### Formate data for heatmap
dataframe <- results.cov

sorted.df <- dataframe[order(dataframe$log2FoldChange, decreasing = TRUE),]

dataframe <- sorted.df[c(4,5,6,7,8,9,10,11,12,13,14,15)]

heatmap(as.matrix(dataframe))


ggplot(test, aes(x = variable, y = gene)) +
  geom_tile(aes(fill = value)) +
  #scale_fill_gradientn(colors = c('blue', 'white', 'red'), limits = c(-3.3,3.3)) +
  #scale_fill_gradientn(colors = c('green', 'black', 'red'), limits = c(-3.1,3.1)) +
  scale_fill_gradientn(colors = c('green', 'black', 'red')) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(position = "top") +
  labs(title = 'Differential Expression') +
  theme(#axis.text.x = element_text(size = 12),
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 11),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    panel.border = element_rect(color = "black", fill = NA)) 

##########################################################################
### Compare KO data to ALS line data

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Input files/")

results.cov <- read.csv('DEGs - ALS 0.05 pv.csv', row.names= 1)
results.ko <- read.csv('DEGs - KO 0.05 pv.csv', row.names = 1)

### Find the overlap

shared.genes <- intersect(results.cov$Gene, results.ko$Gene)

length(shared.genes)

### Draw Venn Diagram
draw.pairwise.venn(area1 = nrow(results.cov), area2 = nrow(results.ko), cross.area = length(shared.genes), 
                   category = c("Mutant vs. Wildtype", "Knockout vs. Wildtype"),
                   lty = rep("blank", 2), fill = c("light blue", "pink"), 
                   alpha = rep(0.5, 2), cat.pos = c(0, 0), 
                   cat.dist = rep(0.025, 2))
                                                                                                                                                        0), cat.dist = rep(0.025, 2))

### Find expression values
shared.pos <- match(shared.genes, results.cov$Gene)
results.cov.shrd <- results.cov[shared.pos,]

shared.pos <- match(shared.genes, results.ko$Gene)
results.ko.shrd <- results.ko[shared.pos,]

### Create data frame of shared genes
names(results.cov.shrd)[1:3] <- c('Mut: Mean', 'Mut: lFC', 'Mut: pv')
names(results.ko.shrd)[1:3] <- c('KO: Mean', 'KO: lFC', 'KO: pv')

results.cov.shrd <- results.cov.shrd[order(row.names(results.cov.shrd)),]
results.ko.shrd <- results.ko.shrd[order(row.names(results.ko.shrd)),]

results.shared.genes <- cbind(results.cov.shrd, results.ko.shrd)
results.shared.genes <- results.shared.genes[c(12,13,1,2,3,14,15,16)]
shared.genes.order <- order(results.shared.genes$`Mut: pv`,decreasing = TRUE)
results.shared.genes <- results.shared.genes[shared.genes.order,]
dir.check <- sign(results.shared.genes$`Mut: lFC`*results.shared.genes$`KO: lFC`)==1
results.shared.genes$same.dir <- dir.check

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/DEG analyses/")

write.csv(results.shared.genes, 'Shared genes - 0.05 padj.csv')

### Sort by direction

up.reg <- results.cov.shared$log2FoldChange > 0
als.up <- results.cov.shared[up.reg,]
als.down <- results.cov.shared[!up.reg,]

up.reg <- results.ko.shared$log2FoldChange > 0
ko.up <- results.ko.shared[up.reg,]
ko.down <- results.ko.shared[!up.reg,]

### Find Shared by direction

shared.up <- intersect(als.up$Gene, ko.up$Gene)
shared.down <- intersect(als.down$Gene, ko.down$Gene)

shared.up.pos <- match(shared.up, als.up$Gene)
shared.up.df <- als.up[shared.up.pos,]

shared.down.pos <- match(shared.down, als.down$Gene)
shared.down.df <- als.down[shared.down.pos,]

### Save results

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Prasanthi/")

write.csv(shared.up.df, "Shared UP reg - ALS-KO.csv")
write.csv(shared.down.df, "Shared DOWN reg - ALS-KO.csv")
write.csv(results.cov.shared, 'Shared DEGs - ALS-KO.csv')

