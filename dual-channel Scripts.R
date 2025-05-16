"
@Author: AmirMohammad Khosravi
@Title: Microarray Analysis GSE119639 
"

#packages:
install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("limma")
#BiocManager::install(version = "3.18")#for using getGEO we need 3.18 versioninstall.packages("BiocManager")


setwd("C:/Users/khosravi/OneDrive/Desktop/MetaAnalysis_khosravi/Microarray_dual-channel/microarray-dual-channel-analysis")

library(limma)
#########################################
#load the data
#########################################
# Set your working directory where raw files are saved

setwd("./GSE19712_RAW")
#prepare targets file manually
targets <- read.delim("targets.txt", sep = "\t")
head(targets)
#read the data:
RG <- read.maimages(files = targets$FileName,
                    source = "agilent",
                    green.only = FALSE)
#########################################
#Normalization
#########################################