setwd('C:\\Users/riocx/Documents/Masters new/Spring 2023/Applied computational genomics/Applied comp project')
input <- "GSE112978_gene_read_counts.txt"
dfgene <- read.table(input,
                     header = T,
                     sep = '\t',
                     row.names = 1)

dfgene_exp1 <- dfgene[,c(6,7,8,12,13,14)]
dfgene_wildType <- dfgene[,c(6,7,8)]

# install.packages("DESeq2")
library(DESeq2)

table2 <- data.frame(name=colnames(dfgene_exp1), condition = c('control','control', 'control', 'treatment','treatment', 'treatment'))
dds_1 <- DESeqDataSetFromMatrix(dfgene_exp1, 
                                colData=table2,design= ~ condition)
# row.names(dds_1) <- dfExp
dds_1 <- DESeq(dds_1)

#normalized reads count
norCounts_1 <- counts(dds_1, normalized=TRUE)

#DEseq results
res_1 <- results(dds_1, alpha=0.01, contrast = c("condition", "treatment", "control"))
#split 2
input <- "GSE112978_gene_read_counts.txt"
dfgene <- read.table(input,
                     header = T,
                     sep = '\t',
                     row.names = 1)

dfgene_ok1485 <- dfgene[,c(12,13,14)]
#Combine wildtype with 0k1485
dfgene_exp2 <- cbind.data.frame(dfgene_wildType, dfgene_ok1485)
# install.packages("DESeq2")
library(DESeq2)

table3 <- data.frame(name=colnames(dfgene_exp2), condition = c('control','control', 'control', 'treatment','treatment', 'treatment'))
dds_2 <- DESeqDataSetFromMatrix(dfgene_exp2, 
                                colData=table3,design= ~ condition)
# row.names(dds_1) <- dfExp
dds_2 <- DESeq(dds_2)

#normalized reads count
norCounts_2 <- counts(dds_2, normalized=TRUE)

#DEseq results
res_2 <- results(dds_2, alpha=0.01, contrast = c("condition", "treatment", "control"))

#PQM.txt file
input_pqm <- "GSE112978_gene_read_counts_pqm1.txt"
dfgene2 <- read.table(input_pqm,
                     header = T,
                     sep = '\t',
                     row.names = 1)

dfgene_ok485 <- dfgene[,c(6,7,8)]
#combine wildtype with ok_485
dfgene_exp3 <- cbind.data.frame(dfgene_wildType, dfgene_ok485)

dfgene_rhd90 <- dfgene[,c(9,10,11)]
#combine wildtype with rhd90
dfgene_exp4 <- cbind.data.frame(dfgene_wildType, dfgene_rhd90)

dfgene_ok485_ok1485 <- dfgene[,c(12,13,14)]
#combine wildtype with ok485_ok1485
dfgene_exp5 <- cbind.data.frame(dfgene_wildType, dfgene_ok485_ok1485)

dfgene_rhd90_ok1485 <- dfgene[,c(15,16,17)]
#combine wildtype with rhd90_ok1485
dfgene_exp6 <- cbind.data.frame(dfgene_wildType, dfgene_rhd90_ok1485)

# install.packages("DESeq2")
library(DESeq2)

table4 <- data.frame(name=colnames(dfgene_exp3), condition = c('control','control', 'control', 'treatment','treatment', 'treatment'))
dds_3 <- DESeqDataSetFromMatrix(dfgene_exp3, 
                                colData=table4,design= ~ condition)
# row.names(dds_1) <- dfExp
dds_3 <- DESeq(dds_3)

#normalized reads count
norCounts_3 <- counts(dds_3, normalized=TRUE)

#DEseq results
res_3 <- results(dds_3, alpha=0.01, contrast = c("condition", "treatment", "control"))
head(res_3)

#for 4th
library(DESeq2)

table5 <- data.frame(name=colnames(dfgene_exp4), condition = c('control','control', 'control', 'treatment','treatment', 'treatment'))
dds_4 <- DESeqDataSetFromMatrix(dfgene_exp4, 
                                colData=table5,design= ~ condition)
# row.names(dds_1) <- dfExp
dds_4 <- DESeq(dds_4)

#normalized reads count
norCounts_4 <- counts(dds_4, normalized=TRUE)

#DEseq results
res_4 <- results(dds_4, alpha=0.01, contrast = c("condition", "treatment", "control"))
head(res_4)
#for 5th
library(DESeq2)

table6 <- data.frame(name=colnames(dfgene_exp5), condition = c('control','control', 'control', 'treatment','treatment', 'treatment'))
dds_5 <- DESeqDataSetFromMatrix(dfgene_exp5, 
                                colData=table6,design= ~ condition)
# row.names(dds_1) <- dfExp
dds_5 <- DESeq(dds_5)

#normalized reads count
norCounts_5 <- counts(dds_5, normalized=TRUE)

#DEseq results
res_5 <- results(dds_5, alpha=0.01, contrast = c("condition", "treatment", "control"))
head(res_5)

#for 6th
library(DESeq2)

table7 <- data.frame(name=colnames(dfgene_exp6), condition = c('control','control', 'control', 'treatment','treatment', 'treatment'))
dds_6 <- DESeqDataSetFromMatrix(dfgene_exp6, 
                                colData=table7,design= ~ condition)
# row.names(dds_1) <- dfExp
dds_6 <- DESeq(dds_6)

#normalized reads count
norCounts_6 <- counts(dds_6, normalized=TRUE)

#DEseq results
res_6 <- results(dds_6, alpha=0.01, contrast = c("condition", "treatment", "control"))
head(res_6)

#combining all the DEseq results together
res_combined <- merge(as.data.frame(res_1), as.data.frame(res_2), as.data.frame(res_3), as.data.frame(res_4), as.data.frame(res_5), as.data.frame(res_6), by = "row.names")

#which(!is.na(res_1<0.1)) ??

# Extract log2FoldChange and padj columns from each results object and # Combine data frames into a single data frame
res1_df <- results(res_1)[, c("log2FoldChange", "padj")]
res2_df <- results(res_2)[, c("log2FoldChange", "padj")]
res3_df <- results(res_3)[, c("log2FoldChange", "padj")]
res4_df <- results(res_4)[, c("log2FoldChange", "padj")]
res5_df <- results(res_5)[, c("log2FoldChange", "padj")]
res6_df <- results(res_6)[, c("log2FoldChange", "padj")]
