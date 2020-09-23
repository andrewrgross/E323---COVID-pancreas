# E323 Boxplot Plotting

### Header
library (ggplot2)

### Import Data

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Vicky/E283-COVID Infection/Round 2 RNAseq/Analysis/")
viral.reads <- read.csv

mock.d1 = c(0,1,0)
infected.d1 = c(3530, 3751, 4114)
mock.d3 = c(0,0,0)
infected.d3 = c(20335,24292,25956)

reads = c(mock.d1, infected.d1, mock.d3, infected.d3)

condition = c("Day 1, Mock", "Day 1, Mock", "Day 1, Mock", 
                 "Day 1, Infected", "Day 1, Infected", "Day 1, Infected",
                 "Day 3, Mock", "Day 3, Mock", "Day 3, Mock", 
                 "Day 3, Infected", "Day 3, Infected", "Day 3, Infected")

data <- data.frame(mock.d1, infected.d1, mock.d3, infected.d3)
data <- data.frame(reads, condition)


ggplot(data = data, aes(x = condition, y = reads, color = condition)) +
  geom_boxplot(size = 1.5, fill = "Grey80") +
  scale_x_discrete(limits = c("Day 1, Mock", "Day 1, Infected", "Day 3, Mock", "Day 3, Infected")) +
  scale_y_continuous(breaks = c(0,5000,10000,15000,20000, 25000)) +
  scale_color_manual(values = c("#fab9b6", "#a6acf7", "#fe1c1c", "#000dc4")) +
  labs(title = "Sequenced Reads by Condition", 
       x = "Condition", 
       y = "Number of Reads") + 
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12))

  



(g <- ggplot(count.pca.df, aes(PC1, PC2, color = time, fill = dose, shape = shape)) + geom_point(size = 4, stroke = 3) +
    geom_text(aes(label=row.names(count.pca.df)),hjust=0,vjust=0) +
    xlab(paste0("PC1: ",round(eig.val$variance.percent[1]),"% variance")) +
    ylab(paste0("PC2: ",round(eig.val$variance.percent[2]),"% variance")) +
    xlim(min(count.pca.df$PC1)-5,max(count.pca.df$PC1)+20) +
    theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
          axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
          axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
          panel.background = element_rect(fill = 'white', color = 'black'),
          plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12)) +
    ggtitle(paste(paste(perterbagen.to.include, collapse = ', '), '-', paste(cell.type.to.include))) + 
    scale_color_gradient(low="blue", high="red") +
    scale_shape_manual(values = c(21, 25)) +
    scale_fill_gradient(low="white", high="orange"))


# Boxplot of MPG by Car Cylinders
boxplot(mpg~cyl,data=mtcars, main="Car Milage Data",
        xlab="Number of Cylinders", ylab="Miles Per Gallon") 

library(ggplot2)

ToothGrowth$dose <- as.factor(ToothGrowth$dose)
head(ToothGrowth)

# Basic box plot
p <- ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_boxplot()
p
# Rotate the box plot
p + coord_flip()
# Notched box plot
ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_boxplot(notch=TRUE)
# Change outlier, color, shape and size
ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)
