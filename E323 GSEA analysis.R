### GAGE Pathway Analysis -- Andrew R Gross -- 2018-08-17
### Application of the GAGE pathway analysis pipeline
### This analysis is downstream of the DESeq2 Pipeline
### INPUT: Differential Expression table with ENTREZ IDs
### OUTPUT: Predicted pathways tables, predicted pathways bar graphs

####################################################################################################################################################
### Header
library(GSEABase)
library(hgu95av2.db)
library(GO.db)
library(ggplot2)
library(gage)

####################################################################################################################################################
### Functions
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

####################################################################################################################################################
### Input
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/Analysis/Gene Enrichment Analysis/')
go.d3.up <- read.csv('GO - Day 3 - up.csv')
go.cov.d3.up <- read.csv('COVID-19_Related_Gene_Sets_table - d3 upalt.txt', sep = '\t', header = TRUE)
go.d3.down <- read.table('GO - Day 3 - down.txt', sep = '\t', header = TRUE)
go.cov.d3.down <- read.csv('COVID-19_Related_Gene_Sets_table - d3 down alt2.txt', sep = '\t', header = TRUE)


(go.d3.up <- go.d3.up[order(go.d3.up$Adjusted.P.value),][1:10,])
(go.d3.down<-go.d3.down[order(go.d3.down$Adjusted.P.value),][1:10,])
(go.cov.d3.up <- go.cov.d3.up[order(go.cov.d3.up$Adjusted.P.value),][1:7,])
(go.cov.d3.down <- go.cov.d3.down[order(go.cov.d3.down$Adjusted.P.value),][1:7,])


### Format
names(go.d3.up)[1] <- 'Term' 
go.d3.up <- go.d3.up[1:10,]

### Split rows
split.positions = c(100,100,100,100,134,129,130)
dataframe <- go.cov.d3.down
for(row.num in 1:nrow(dataframe)){
  dataframe$Row1[row.num] = substr(dataframe$Term[row.num], 1, split.positions[row.num])
  dataframe$Row2[row.num] = substr(dataframe$Term[row.num], split.positions[row.num]+1, 200)
}

break.across.rows <- function(dataframe,split.positions){
  for(row.num in 1:nrow(dataframe)){
    dataframe$Row1[row.num] = substr(dataframe$Term[row.num], 1, split.positions[row.num])
    dataframe$Row2[row.num] = substr(dataframe$Term[row.num], split.positions[row.num]+1, 200)
  }
  return(dataframe)
}
  
go.cov.d3.down <- break.across.rows(go.cov.d3.down,split.positions = c(54,100,100,100,134,129,130))  

break.across.rows(go.d3.up, c(100,100,100,100,44,129,130,129,130,129))$Row1
go.d3.up <- break.across.rows(go.d3.up, c(100,100,100,100,44,129,130,129,130,129))

go.d3.down <- break.across.rows(go.d3.down,  c(100,100,100,100,44,129,130,129,130,129))


### Plot

ggplot(data = go.d3.up, aes(x = Term, y = -log10(Adjusted.P.value))) +
  geom_bar(stat = 'identity', fill = 'red4') +
  geom_text(aes(label=Overlap), hjust=-0.3, vjust = 1.5, color="Blue", size=4) +
  geom_text(aes(label=paste("E",round(log10(Adjusted.P.value)))), hjust=-0.5, vjust = -1.4, color="black", size=4) +
  geom_text(aes(label=Row1, y = 0.5, hjust=0), vjust = -0.5, color = "White", fontface = "bold") +
  geom_text(aes(label=Row2, y = 0.5, hjust=0), vjust = 1.5, color = "white", fontface = "bold") +
  scale_x_discrete(limits=rev(go.d3.up$Term)) +
  ylim(0,24) +
  coord_flip() + 
  labs(title = 'Enriched Gene Ontology Biological Processes',
       subtitle = 'Day 3, Up-regulated in Infected', 
       y = 'P-value, Adjusted') +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(color = "black", fill = NA)) 


ggplot(data = go.d3.down, aes(x = Term, y = -log10(Adjusted.P.value))) +
  geom_bar(stat = 'identity', fill = 'blue4') +
  geom_text(aes(label=Overlap), hjust=1.3, vjust = 1.5, color="Blue", size=4) +
  geom_text(aes(label=paste("E",round(log10(Adjusted.P.value)))), hjust=1.5, vjust = -1.2, color="black", size=4) +
  geom_text(aes(label=Row1, y = 24, hjust=0), vjust = -0.5, color = "black", fontface = "bold") +
  geom_text(aes(label=Row2, y = 24, hjust=0), vjust = 1.5, color = "black", fontface = "bold") +
  scale_x_discrete(limits=rev(go.d3.down$Term)) +
  ylim(24,0) +
  coord_flip() + 
  labs(title = 'Enriched Gene Ontology Biological Processes',
       subtitle = ' Day 3, Down-regulated in Infected',
       y = 'P-value, Adjusted') +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(color = "black", fill = NA)) 

  
## COV, Up

ggplot(data = go.cov.d3.up, aes(x = Term, y = -log10(Adjusted.P.value))) +
  geom_bar(stat = 'identity', fill = 'red4') +
  geom_text(aes(label=Overlap), hjust=-0.3, vjust = 1.5, color="Blue", size=4) +
  geom_text(aes(label=paste("E",round(log10(Adjusted.P.value)))), hjust=-0.5, vjust = -1.4, color="black", size=4) +
  #geom_label(aes(label=Term, y = 0.5, hjust=0), vjust = -0.1, color = "black", fill = 'white', alpha = 0.8, label.size = 0, fontface = "bold") +
  geom_text(aes(label=Row1, y = 0.5, hjust=0), vjust = -0.5, color = "White", fontface = "bold") +
  geom_text(aes(label=Row2, y = 0.5, hjust=0), vjust = 1.5, color = "white", fontface = "bold") +
  scale_x_discrete(limits=rev(go.cov.d3.up$Term)) +
  ylim(0,29) +
  coord_flip() + 
  labs(title = 'Enriched COVID19-associated transcripts',
       subtitle = 'Day 3, Upregulated in Infected',
       y = 'P-value, Adjusted') +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(color = "black", fill = NA)) 

## COV, Down

ggplot(data = go.cov.d3.down, aes(x = Term, y = -log10(Adjusted.P.value))) +
  geom_bar(stat = 'identity', fill = 'blue4') +
  geom_text(aes(label=Overlap), hjust=1.3, vjust = 1.5, color="Blue", size=4) +
  geom_text(aes(label=paste("E",round(log10(Adjusted.P.value)))), hjust=1.5, vjust = -1.2, color="black", size=4) +
  geom_text(aes(label=Row1, y = 29, hjust=0), vjust = -0.5, color = "black", fontface = "bold") +
  geom_text(aes(label=Row2, y = 29, hjust=0), vjust = 1.5, color = "black", fontface = "bold") +
  scale_x_discrete(limits=rev(go.cov.d3.down$Term)) +
  ylim(29,0) +
  coord_flip() + 
  labs(title = 'Enriched COVID19-associated transcripts',
       subtitle = 'Day 3, Down-regulated in Infected',
       y = 'P-value, Adjusted') +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(color = "black", fill = NA)) 



#########################################################
### Scratchwork

### Text above bars
ggplot(data = go.cov.d3.down, aes(x = Term, y = -log10(Adjusted.P.value))) +
  geom_bar(stat = 'identity', fill = 'blue4', width = 0.5) +
  geom_text(aes(label=Overlap), hjust=-0.3, color="Blue", size=4.5) +
  geom_text(aes(label=paste("E",round(log10(Adjusted.P.value)))), hjust=1.3, color="white", size=4.5) +
  geom_text(aes(label=Term, y = 0, hjust=0, vjust = -1.9), color = "black") +
  scale_x_discrete(limits=rev(go.cov.d3.down$Term)) +
  ylim(0,40) +
  coord_flip() + 
  labs(title = 'Enriched COVID19-associated transcripts: Day 3, Downregulated in Infected') +
  theme(#axis.text.x = element_text(size = 12),
    axis.text.y = element_blank(),
    #axis.text.y = element_text(size = 11),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    panel.border = element_rect(color = "black", fill = NA)) 



#  scale_fill_gradientn(colors = c('green', 'black', 'red')) +
#  scale_y_discrete(position = "right") +
#  scale_x_discrete(position = "top") +
#  labs(title = 'Differential Expression') 

### Composing list of relevant transcripts
go.d3.up$Genes[1]
go.d3.down$Genes

strsplit(go.d3.up$Genes[1],';')
transcript.list <- c()

for (transcripts in go.d3.up$Genes) {
  transcript.list [length(transcript.list)+1] <- strsplit(transcripts, ';')
}



summary(transcript.list)
intersect(transcript.list[[1]],transcript.list[[7]])
setdiff(transcript.list[[10]],transcript.list[[9]])


transcript.list.concat <- c()

for (element in transcript.list) {
  transcript.list.concat <- union(transcript.list.concat, element)
}

#union(transcript.list[[1]],transcript.list[[2]],transcript.list[[7]],transcript.list[[8]],transcript.list[[9]])
transcript.desc <- results.cov3.500[match(transcript.list.concat,results.cov3.500$Gene),]


setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/Analysis/Gene Enrichment Analysis/')
write.csv(transcript.desc,'GO - D3 - Up -- Transcripts.csv')
write.csv(transcript.desc,'GO - D3 - Down -- Transcripts.csv')


#### Scratchwork

df <- data.frame(dose=c("D0.5", "D1", "D2"),
                 len=c(4.2, 10, 29.5))
head(df)

p<-ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity")
p

# Horizontal bar plot
p + coord_flip()

ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=len), vjust=-0.3, size=3.5)+
  theme_minimal()
# Inside bars
ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=len), vjust=1.6, color="white", size=3.5)+
  theme_minimal()

####################################################################################################################################################
### Input
data(gse16873)
data(kegg.gs)
data(go.gs)

### Differential Expression Results
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/Analysis/Differential Expression/')
results.cov3.up <- read.csv('DEGs - COV Day 3 0.05 padj - UP.csv')

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Input files/")

results.als.ez <- read.csv('DEGs - ALS 0.05 padj ENTREZ.csv', row.names= 1)
results.ko.ez <- read.csv('DEGs - KO 0.05 padj ENTREZ.csv', row.names = 1)

results.als.ez <- read.csv('DEGs - ALS 0.05 pv ENTREZ.csv', row.names= 1)
results.ko.ez <- read.csv('DEGs - KO 0.05 pv ENTREZ.csv', row.names = 1)

### Full expression tables
setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/Motor Neurons/")
results.als <- read.csv('E099 - PM - ALS v CTRL/PM-5119--07--02--2018_TPM.csv', row.names = 1)


setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/Motor Neurons/')
tpm.als <- read.csv('E099 - PM - ALS v CTRL/PM-5119--07--02--2018_TPM.csv', row.names = 1)
names(tpm.als) <- c('02iCTR','03iCTR','159iALS','172iCTR','372iALS_n1','372iALS_n2','372iALS_n3','395iCTR')

tpm.ko <- read.csv('E0x - PP - D10 KO v CTRL/TPM.csv', row.names = 1)
names(tpm.ko) <- read.csv('E0x - PP - D10 KO v CTRL/Sample names.txt', header = FALSE)[,1]

####################################################################################################################################################
### Formatting
### Reorder columns
tpm.als <- tpm.als[c(1,2,4,8,5,6,7,3)]
#head(results.als.ez)

####################################################################################################################################################
### Subsetting
### Remove universally low scoring genes

row.max <- apply(tpm.als,1,max)
#row.is.zero <- row.max==0
#quantile(row.max,c(0.1,0.2,0.3))
#tpm.als <- tpm.als[row.max>0,]
(cutoff = quantile(row.max,0.5))
rows.to.keep <- row.max > cutoff
tpm.als <- tpm.als[rows.to.keep,]
nrow(tpm.als)
summary(tpm.als)

### Renaming
test <- convertIDs(tpm.als)

test$DESCRIPTION <- NA
output <- test[c(9,10,1:8)]

### Generate Phenotype file
phenotype <- c('CTR','CTR','CTR','CTR','ALS','ALS','ALS','ALS')
phenotype.text <- '8\t2\t1\n#CTR\tALS\n0 0 0 0 1 1 1 1'

### KNOCKOUT

row.max <- apply(tpm.ko,1,max)
#row.is.zero <- row.max==0
#quantile(row.max,c(0.1,0.2,0.3))
#tpm.als <- tpm.als[row.max>0,]
(cutoff = quantile(row.max,0.5))
rows.to.keep <- row.max > cutoff
tpm.ko <- tpm.ko[rows.to.keep,]
nrow(tpm.ko)
summary(tpm.ko)

### Renaming
test <- convertIDs(tpm.ko)

test$DESCRIPTION <- NA
output.ko <- test[c(7,8,1:6)]

### Generate Phenotype file

phenotype.text.ko <- '6\t2\t1\n#KO\tCTR\n0 0 0 0 1 1 1 1'

####################################################################################################################################################
### Rank lists

output.deg.ko <- results.ko.ez[c(10,3)]
output.deg.ko$logpvalue <- -log(output.deg.ko$pvalue)
output.deg.ko <- output.deg.ko[c(1,3)]

####################################################################################################################################################
### Output

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/GSEA/Input to GSEA/")

write.table(output, 'ALS expression data.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(phenotype.text, 'ALS phenotype.cls', col.names = FALSE, quote = FALSE, row.names = FALSE)


write.table(output.ko, 'KO expression data.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(phenotype.text.ko, 'KO phenotype.cls', col.names = FALSE, quote = FALSE, row.names = FALSE)


write.table(output.deg.ko, 'KO deg ranked list.rnk', row.names = FALSE, sep = '\t', quote = FALSE)







#results.als.ez <- results.als.ez[results.als.ez$padj<0.05,]
for.gage.als <- as.matrix(results.als.ez[c(4,5,6,7,8,9,10,11)])

#results.ko.ez <- results.ko.ez[results.ko.ez$padj<0.05,]
for.gage.ko <- as.matrix(results.ko.ez[c(7,8,9,4,5,6)])

for.gage.loadings <- as.matrix(tpm.o.i.ez2[c(1,2,4,8,3,5,6,7)])
for.gage.loadings <- as.matrix(tpm.null[c(1,2,4,8,3,5,6,7)])

##########################################################################
### Confirm ID type for gene set
lapply(kegg.gs[1:3],head)
head(kegg.gs[1])

head(go.gs[1])    # They all appear to be Entrez from googling

##########################################################################
### Run analysis

### Example analysis
#cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = ref.idx, samp = samp.idx, compare ="unpaired")
#gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, ref = hn, samp = dcis)
##go.gs here only the first 1000 entries as a fast example.
#gse16873.go.p <- gage(gse16873, gsets = go.gs, ref = hn, samp = dcis)
##or we can do analysis on BP, MF, or CC subcategories if we've
##generated the data as above.
##gse16873.bp.p <- gage(gse16873, gsets = go.bp,
## ref = hn, samp = dcis)
#str(gse16873.kegg.p, strict.width='wrap')
#gse16873.kegg.2d.p <- gage(gse16873, gsets = kegg.gs, + ref = hn, samp = dcis, same.dir = F)

#results.kegg.2d.p <- gage(results.als.ez, gsets = kegg.gs, + ref = ctrl, samp = als, same.dir = F)

### Run real data
ctrl.index = 1:4
als.index = 5:8
pway.als.bd <- gage(for.gage.als, gsets = kegg.gs, ref = ctrl.index, samp = als.index, compare = "unpaired", same.dir=F)
pway.als.bd <- pway.als.bd$greater
pway.als <- gage(for.gage.als, gsets = kegg.gs, ref = ctrl.index, samp = als.index, compare = "unpaired")
pway.als.up <- pway.als$greater
pway.als.down <- pway.als$less

#head(pway.ALS$greater,20)
#head(pway.ALS$less)

ctrl.index = 1:3
ko.index = 4:6
pway.ko.bd <- gage(for.gage.ko, gsets = kegg.gs, ref = ctrl.index, samp = ko.index, compare = "unpaired", same.dir=F)
pway.ko.bd <- pway.ko.bd$greater
pway.ko <- gage(for.gage.ko, gsets = kegg.gs, ref = ctrl.index, samp = ko.index, compare = "unpaired")
pway.ko.up <- pway.ko$greater
pway.ko.down <- pway.ko$less

ctrl.index = 1:4
als.index = 5:8
pway.loadings.bd <- gage(for.gage.loadings, gsets = kegg.gs, ref = ctrl.index, samp = als.index, compare = "unpaired", same.dir=F)
pway.loadings.bd <- gage(for.gage.loadings, gsets = kegg.gs)

pway.loadings.bd <- pway.loadings.bd$greater
head(pway.loadings.bd)

head(pway.ko.bd)
head(pway.ko.up)
head(pway.ko.down)
######################################################################################################################################
### Output data

############################################################
### Plot function
gen.pway.bar.chart <- function(data.frame, experiment, direction) {
  if(direction == 'Up') {color = 'Red'}
  if(direction == 'Down') {color = 'Blue'}
  if(direction == 'Bidir') {color = 'Purple'; direction = 'Bidirectional'}
  plot.df <- as.data.frame(data.frame)[1:20,]
  title = paste('Disreg. Pathways in', experiment, ' - ', direction)
  
  plot.df <- data.frame(Pathways = row.names(plot.df), log.p.val = -log(plot.df$p.val))
  plot.df$Pathways <- factor(plot.df$Pathways, levels = rev(plot.df$Pathways))
  
  plot <- ggplot(plot.df, aes(Pathways,log.p.val)) +
    geom_col(fill = color) + coord_flip() + theme_bw() + labs(title = title)
  return(plot)
}

experiment = 'Mutant'
direction = 'Up'
plot <- gen.pway.bar.chart(pway.als.up, paste(experiment, 'vs. WT'), direction)
plot

direction = 'Down'
plot <- gen.pway.bar.chart(pway.als.down, paste(experiment, 'vs. WT'), direction)
plot

direction = 'Bidir'
plot <- gen.pway.bar.chart(pway.als.bd, paste(experiment, 'vs. WT'), direction)
plot

experiment = 'KO'
direction = 'Up'
plot <- gen.pway.bar.chart(pway.ko.up, paste(experiment, 'vs. WT'), direction)
plot

direction = 'Down'
plot <- gen.pway.bar.chart(pway.ko.down, paste(experiment, 'vs. WT'), direction)
plot

direction = 'Bidir'
plot <- gen.pway.bar.chart(pway.ko.bd, paste(experiment, 'vs. WT'), direction)
plot

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Pathway analyses/Pathway bar charts - 0.05 pv/')

output.name = paste0('Pway bar chrt - ', experiment, ' - ', direction, '.png')
png(filename = output.name,
    width = 1000, height = 800, units = "px", pointsize = 20,
    bg = "white",  res = 120)
plot
dev.off()

############################################################
### Examine significant gene sets
test <- sigGeneSet(pway.loadings.bd)

############################################################
### Manual plotting
plot.df <- as.data.frame(pway.als.up)[1:10,]
title = 'Disreg. Pathways in ALS lines - Up'; color = 'Red'
plot.df <- as.data.frame(pway.als.down)[1:10,]
title = 'Disreg. Pathways in ALS lines - Down'; color = 'Blue'

plot.df <- as.data.frame(pway.ko.up)[1:10,]
title = 'Disreg. Pathways in KO lines - Up'; color = 'Red'
plot.df <- as.data.frame(pway.ko.down)[1:10,]
title = 'Disreg. Pathways in KO lines - Down'; color = 'Blue'


plot.df <- data.frame(Pathways = row.names(plot.df), log.p.val = -log(plot.df$p.val))
### Set order
plot.df$Pathways <- factor(plot.df$Pathways, levels = rev(plot.df$Pathways))

ggplot(plot.df, aes(Pathways,log.p.val)) +
  geom_col(fill = color) + coord_flip() + theme_bw() + labs(title = title)

ggplot(plot.df, aes(Pathways,log.p.val)) +
  geom_col(fill = 'Purple') + coord_flip() + theme_bw() + labs(title = 'Disreg. Pathways in ALS lines - Bidirectional')


####################################
### Output
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Pathway analyses/')

write.csv(pway.als.bd, 'pways - ALS 0.05pv - Bidir.csv')
write.csv(pway.als.up, 'pways - ALS 0.05pv - Up.csv')
write.csv(pway.als.down, 'pways - ALS 0.05pv - Down.csv')
write.csv(pway.ko.bd, 'pways - ko 0.05pv - Bidir.csv')
write.csv(pway.ko.up, 'pways - ko 0.05pv - Up.csv')
write.csv(pway.ko.down, 'pways - ko 0.05pv - Down.csv')





##########################################################################
### Pathview Analysis
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Pathway analyses/Mapped pathways/')

fc.data <- as.matrix(results.als.ez)[, 2] ; suffix = 'als - 0.05pv' 
fc.data <- as.matrix(results.ko.ez)[, 2] ; suffix = 'ko - 0.05pv' 

pway.id = '04145'
pway.id = '00190'
pway.id = '04512'
pway.id = '04722'
pway.id = '04540'

pv.out <- pathview(gene.data = fc.data, pathway.id = pway.id, species = "hsa", out.suffix = suffix)

pv.out <- pathview(gene.data = fc.data, pathway.id = pway.id, species = "hsa", out.suffix = suffix)
