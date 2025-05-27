# @Author: AmirMohammad Khosravi
# @Title: Microarray Analysis of GSE115800 for just two T24 samples (GSM3190198, GSM3190199)

# Install packages (only run once if not installed)
install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("BiocGenerics")

#########################################
# Load the data
#########################################

# Set working directory
setwd("C:/Users/Khosravi/Desktop/micro-dual/microarray-dual-channel-analysis-main")

# Load required libraries
library(limma)
library(BiocGenerics)

# Load the targets file (targets.txt)
targets <- read.delim("targets.txt", sep = "\t")
head(targets)

# Read the image files
RG <- read.maimages(files = targets$FileName,
                    source = "agilent",
                    green.only = FALSE)

#########################################
# Raw Data Quality Control (QC)
#########################################

# 1. Boxplot
boxplot(RG$R, main = "Raw Red Channel", las = 2)
boxplot(RG$G, main = "Raw Green Channel", las = 2)

# 2. Density plots
par(mfrow = c(1, 2))
plotDensities(RG$R, main = "Density of Raw Red Channel")
plotDensities(RG$G, main = "Density of Raw Green Channel")

# 3. MA plots before normalization
par(mfrow = c(1, 1))
for (i in 1:ncol(RG$R)) {
  ma <- limma::MA.RG(RG[, i])
  limma::plotMA(ma, main = paste("MA-Plot Array", i), ylim = c(-5, 5))
}

#########################################
# Normalization and Background Correction
#########################################

# Background correction
RGb <- backgroundCorrect(RG, method = "normexp")

# Within-array normalization
MA <- normalizeWithinArrays(RGb, method = "loess")

#########################################
# Display and Save Normalized Expression Values (M and A)
#########################################

# Extract M-values
m_df <- as.data.frame(MA$M)

# Add gene names if available
if (!is.null(MA$genes$Name)) {
  rownames(m_df) <- MA$genes$Name
}

# Add gene names from ProbeName column if available
if (!is.null(MA$genes$ProbeName)) {
  rownames(m_df) <- MA$genes$ProbeName
  rownames(a_df) <- MA$genes$ProbeName
}

View(m_df)

# Save M-values to CSV
write.csv(m_df, "normalized_M_values.csv")

# Extract A-values
a_df <- as.data.frame(MA$A)

# Add gene names if available
if (!is.null(MA$genes$Name)) {
  rownames(a_df) <- MA$genes$Name
}

View(a_df)

# Save A-values to CSV
write.csv(a_df, "normalized_A_values.csv")

#########################################
# Statistical Analysis
#########################################

# Fit linear model
fit <- lmFit(MA)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Get top 50 differentially expressed genes
top_genes <- topTable(fit, number = 50)

# Save the results
write.csv(top_genes, file = "top_genes.csv", row.names = TRUE)
