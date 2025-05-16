"
@Author: AmirMohammad Khosravi
@Title: Microarray Analysis of GSE115800 for just two T24 samples(GSM3190198,GSM3190199) 
"

#packages:
install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("limma")



#########################################
#load the data
#########################################
# Set your working directory where raw files are saved
setwd("C:/Users/khosravi/OneDrive/Desktop/MetaAnalysis_khosravi/Microarray_dual-channel/microarray-dual-channel-analysis")
library(limma)
setwd("./GSE115800_T24samples_RAW")
#prepare targets file manually
targets <- read.delim("targets.txt", sep = "\t")
head(targets)
#read the data:
RG <- read.maimages(files = targets$FileName,
                    source = "agilent",
                    green.only = FALSE)
#########################################
#QC & Normalization
#########################################
#QC:
boxplot(RG$R, main = "Raw Red Channel", las = 2)
boxplot(RG$G, main = "Raw Green Channel", las = 2)
#background correction:
RGb <- backgroundCorrect(RG, method = "normexp")
MA <- normalizeWithinArrays(RGb, method = "loess")

#normalization:
#MA <- normalizeWithinArrays(RG.bg, method = "loess")
#MA <- normalizeBetweenArrays(MA, method = "Aquantile")


#########################################
#Analysis
#########################################
fit <- lmFit(MA)
fit <- eBayes(fit)
top_genes <- topTable(fit, number = 50)
write.csv(top_genes, file = "top_genes.csv", row.names = TRUE)

