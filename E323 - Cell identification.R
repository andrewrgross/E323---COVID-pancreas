### E323 Cell Identification

### Header
library(QICFIT)

### Functions


### Import Data
query.df.full <- read.csv("C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/VW-10228--08--04--2020_FPKM.csv", row.names = 1)

ref.df <- read.csv("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/reference_data/gtex/GTEx_median_expression_scores_by_tissue.csv", row.names = 1)


### Format Data

query.df.full <- convert.ids(query.df.full)
query.df <- query.df.full[1]

### Edit Ensembl IDs and reorder
test <- head(ref.df[1:10])
test[names]

input.df <- ref.df
for(row.num in 1:nrow(input.df)){
  id = ref.df$Name[row.num]
  new.id = strsplit(id,"\\.")[[1]][1]
  input.df$Name[row.num] <- new.id
}

row.names(ref.df) <- input.df$Name
ref.df <- ref.df[3:55]
ref.df <- ref.df[order(row.names(ref.df)),]

### Run Spearman Correlation

spearman.results <- spearman.corr.for.single.sample(query.df = query.df[1], ref.df = ref.df)
spearman.results.ordered <- spearman.results[order(spearman.results[1], decreasing = TRUE),, drop=FALSE]

### Comparing to Pituitary

which(names(ref.df) == "Pancreas")
ref.df.panc <- ref.df[41]#[1:20,, drop = FALSE]

### Generate list of top 100 Genes

ref.panc.1000 <- ref.df.panc[order(ref.df.panc, decreasing = TRUE),,drop=FALSE][1:1000,, drop = FALSE]
query.1000 <- query.df2[order(query.df2, decreasing = TRUE),, drop=FALSE][1:1000,, drop=FALSE]

matches <- match(row.names(ref.panc.100), row.names(query.1000))
head(ref.panc.1000)
head(query.1000)
head(matches)
ref.panc.1000[15,,drop=FALSE]
query.1000[17,,drop = FALSE]

### Find the spearman score of the top 1000?
spearman.results.1000 <- spearman.corr.for.single.sample(query.df = query.1000, ref.df = ref.df)
spearman.results.1000.ordered <- spearman.results.1000[order(spearman.results.1000[1], decreasing = TRUE),, drop=FALSE]

### Find the top 1000 from each sample and compare to the top 1000 of the query
### For each sample, find the top 1000

spearman.corr.for.single.sample.10000 <- function(query.df, ref.df) {   # Calculate the spearman correlation between the query.df and the references
  spearman.results <- data.frame(rep(0,ncol(ref.df)))                             # Generate empty results table
  row.names(spearman.results) <- names(ref.df)                                # Name empty results table
  names(spearman.results) <- names(query.df)
  
  for(ref.tissue.num in 1:ncol(ref.df)) {                                     # Loop through each reference
    ref.tissue.data <- ref.df[ref.tissue.num]                                 # Call the current tissue from the references
    tissue <- names(ref.tissue.data)
    ### Subset to the top 10000 genes
    ref.tissue.data <- ref.tissue.data[order(ref.tissue.data, decreasing = TRUE),,drop=FALSE][1:10000,, drop=FALSE]
    ### Generate a data frame containing the two transcriptomes being compared
    ref.tissue.data <- ref.tissue.data[which(ref.tissue.data[,1]>0),,drop=FALSE]         # Filter out missing values from tissue
    genes.present.in.ref <- row.names(ref.tissue.data)                                   # Declare the genes present in the reference
    genes.missing.in.query <- setdiff(genes.present.in.ref, row.names(query.df))    # Declare genes in reference missing from query.df
    rows.to.add <- data.frame(rep(0,length(genes.missing.in.query)))                       # Generate a zero data frame the size of the missing rows
    row.names(rows.to.add) <- genes.missing.in.query                                       # Name rows after missing rows
    names(rows.to.add) <- names(query.df)                                           # Name column the same as the query
    query.df2 <- rbind(query.df,rows.to.add)                                # Use rbind to make a full data frame containing exactly the genes in the reference
    query.df2 <- query.df2[genes.present.in.ref, , drop = FALSE]           # Reorder query.df to match reference
    spearman.input <- cbind(query.df2, ref.tissue.data)                          # Bind query.df and reference
    result <- rcorr(as.matrix(spearman.input), type = 'spearman')[[1]]                     # Perform spearman calculation using rcorr
    result <- round(result[2], 5)                        # Round result
    spearman.results[tissue,] <- result                        # Add to results table
  }
  return(spearman.results)
}

query.10000 <- query.df2[order(query.df2, decreasing = TRUE),, drop=FALSE][1:10000,, drop=FALSE]

spearman.results.10000 <- spearman.corr.for.single.sample(query.df = query.10000, ref.df = ref.df)
spearman.results.10000.ordered <- spearman.results.10000[order(spearman.results.10000[1], decreasing = TRUE),, drop=FALSE]

spearman.results.1000 <- spearman.corr.for.single.sample(query.df = query.1000, ref.df = ref.df)
(spearman.results.1000.ordered <- spearman.results.1000[order(spearman.results.1000[1], decreasing = TRUE),, drop=FALSE])

spearman.results.1000 <- spearman.corr.for.single.sample.1000(query.df = query.1000, ref.df = ref.df)
(spearman.results.1000.ordered <- spearman.results.1000[order(spearman.results.1000[1], decreasing = TRUE),, drop=FALSE])

spearman.manual.results <- spearman.manual.for.single.sample(query.df,ref.df,weighted = TRUE)
(spearman.manual.results.o <- spearman.manual.results[order(spearman.manual.results[1], decreasing = TRUE),, drop=FALSE])


spearman.genes <- spearman.gene.assessor.single.sample(query.df, ref.df, weighted = FALSE)


spearman.gene.assessor.single.sample <- function(query.df, ref.df, weighted = FALSE) {   # Calculate the spearman correlation between the query.df and the references
  spearman.results <- data.frame(rep(0,ncol(ref.df)))                             # Generate empty results table
  row.names(spearman.results) <- names(ref.df)                                # Name empty results table
  names(spearman.results) <- names(query.df)
  
  for(ref.tissue.num in 1:ncol(ref.df)) {                                     # Loop through each reference
    ref.tissue.data <- ref.df[ref.tissue.num]                                 # Call the current tissue from the references
    tissue <- names(ref.tissue.data)
    ### Generate a data frame containing the two transcriptomes being compared
    ref.tissue.data <- ref.tissue.data[which(ref.tissue.data[,1]>0),,drop=FALSE]         # Filter out missing values from tissue
    genes.present.in.ref <- row.names(ref.tissue.data)                                   # Declare the genes present in the reference
    genes.missing.in.query <- setdiff(genes.present.in.ref, row.names(query.df))    # Declare genes in reference missing from query.df
    rows.to.add <- data.frame(rep(0,length(genes.missing.in.query)))                       # Generate a zero data frame the size of the missing rows
    row.names(rows.to.add) <- genes.missing.in.query                                       # Name rows after missing rows
    names(rows.to.add) <- names(query.df)                                           # Name column the same as the query
    query.df2 <- rbind(query.df,rows.to.add)                                # Use rbind to make a full data frame containing exactly the genes in the reference
    query.df2 <- query.df2[genes.present.in.ref, , drop = FALSE]           # Reorder query.df to match reference
    spearman.input <- cbind(query.df2, ref.tissue.data)                          # Bind query.df and reference
    
    #result <- rcorr(as.matrix(spearman.input), type = 'spearman')[[1]]                     # Perform spearman calculation using rcorr
    #result <- round(result[2], 5)                        # Round result
    #spearman.results[tissue,] <- result                        # Add to results table
    
    #spearman.input.o <- spearman.input[order(spearman.input[1], decreasing = TRUE),, drop=FALSE]
    
    ### Assign rank and delta rank
    spearman.input$sample_rank <- rank(-spearman.input[1])
    spearman.input$ref_rank <- rank(-spearman.input[2])
    spearman.input$d <- spearman.input$sample_rank - spearman.input$ref_rank
    spearman.input$d2 <- (spearman.input$d)^2
    
    if(weighted == TRUE) {
      ### Calculate weighted delta rank
      spearman.input$d2 <- spearman.input$d2 * ((-2*spearman.input$ref_rank/max(spearman.input$ref_rank))+2)
    }
    
    gene.report <- spearman.input[order(spearman.input$d2, decreasing = TRUE),]
    # Add to results table
  }
  return(gene.report)
}

spearman.genes <- spearman.gene.assessor.single.sample(query.1000, ref.df[41], weighted = FALSE)

spearman.genes.g <- addGene(spearman.genes[1:100,])
spearman.genes.g2 <- addGene(tail(spearman.genes,100))
spearman.genes.g2 <- spearman.genes.g2[order(spearman.genes.g2$sample_rank),]
spearman.genes.g2 <- add.description(spearman.genes.g2)
