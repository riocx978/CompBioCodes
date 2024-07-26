setwd('C:\\Users/riocx/Documents/Masters new/Spring 2023/Applied computational genomics/Applied comp project/')
input <- 'C:\\Users/riocx/Documents/Masters new/Spring 2023/Applied computational genomics/Applied comp project/all_samples_exons.txt'
df <- read.table(input,
                 header = T,
                 sep = '\t',
                 row.names = 1)

df <- df[,6:11]
### Check the distribution of counts per sample ----

par(mfrow = c(2,3), mar = c(4,4,1,1))
print( apply(df, 2, function(x) quantile(as.numeric(x))) )
hist.plots <- apply(df, 2, function(x) {hist(log2(x), breaks = 100)})

### Keep genes with expression bigger than 20 ----
expressed.ids <- apply(df, 1, function(x) any(x > 20))
dfExp <- df[expressed.ids, ]

## check distribution again after filtering
par(mfrow = c(2,3), mar = c(4,4,1,1))
print( apply(dfExp, 2, function(x) quantile(as.numeric(x))) )
hist.plots <- apply(dfExp, 2, function(x) {hist(log2(x), breaks = 100)})

#Instal new packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomeInfoDb")
BiocManager::install("DESeq2")
BiocManager::install("RCurl")
##
packages <- c("BiocManager","gplots","knitr","RColorBrewer")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos='https://cloud.r-project.org') 
}

# install our primary analysis-packages from bioconductor if they aren't already installed
bioconductor.packages <- c("GO.db", "DESeq2")
if (length(setdiff(bioconductor.packages, rownames(installed.packages()))) > 0){
  BiocManager::install(bioconductor.packages,
                       ask = FALSE,
                       quiet = TRUE,
                       verbose = FALSE)
}

# clean up leftover variables from above
rm("packages", "bioconductor.packages")

##
#building using DESeq2
library(DESeq2)
table2 <- data.frame(name = colnames(dfExp),
                     condition = c('control','treatment','treatment','treatment','treatment','treatment'))
dds <- DESeqDataSetFromMatrix(dfExp, 
                              colData=table2, design= ~ condition)
dds <- DESeq(dds)
#normalized reads count
norCounts <- counts(dds, normalized=TRUE)
res <- results(dds)
####Getting differentially expressed genes ----
# take a look at our results object ("res")
res
# res is a dataframe object, so we can check out metadata for what the columns mean
def <- mcols(res, use.names=TRUE)
def@listData[["description"]]
## extract the genes with adjusted P-value < 0.01
resSig <- res[ which(res$padj < 0.01), ]
dim(resSig)
#Generating the boxplots 
box_plot_data <- resSig
colnames(box_plot_data) <- c('Wildtype','pqm(rhd90)','pqm(ok485)','ceh60','ceh60_pqm(rhd90)','ceh60_pqm(ok485)')

boxplot(log2(box_plot_data), outline = F, col = c('purple','orange',
                                           'red','navy', 'blue', 'green'),
        border = c('purple','orange','red','navy', 'blue', 'green'),
        medcol = rep('white', 6), ylab = 'log2(FPKM)')

#contrasting the 5 genes against the wildtype
# Extract expression values for the 6 genes
gene_names <- c('Wildtype','pqm(rhd90)','pqm(ok485)','ceh60','ceh60_pqm(rhd90)','ceh60_pqm(ok485)')
gene_expression <- norCounts[gene_names,]

# Create a data frame with gene expression values
df_gene_expression <- data.frame(expression=gene_expression, 
                                 gene=rep(gene_names, each=ncol(gene_expression)))

# Create a boxplot of the differential gene expression (log2) of Wildtype against the other 5 genes
library(ggplot2)
ggplot(df_gene_expression, aes(x=gene, y=log2(Wildtype+1))) + 
  geom_boxplot() + 
  labs(x="Genes", y="Differential gene expression (log2)") + 
  theme_bw()
#####

#Calculating pearsons correlation
fpkm <- read.table(file = 'C:\\Users/riocx/Documents/Masters new/Spring 2023/Applied computational genomics/Applied comp project/aa_samples_exons.fpkm.txt', header = T, sep = '\t', row.names = 1)
corF <- cor(fpkm)
colnames(corF) <- c('Wildtype','pqm(rhd90)','pqm(ok485)','ceh60','ceh60_pqm(rhd90)','ceh60_pqm(ok485)')
rownames(corF) <- c('Wildtype','pqm(rhd90)','pqm(ok485)','ceh60','ceh60_pqm(rhd90)','ceh60_pqm(ok485)')
library(reshape2)
melted_corF <- melt(corF) # change the correlation matrix to the pair 
length(corF)

get_upper_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

upper_tri        <- get_upper_tri(corF)
melted_upper_tri <- melt( upper_tri )

ggplot(data = melted_upper_tri, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "navy", high = "red", mid = "white", 
                       midpoint = 0.9, limit = c(0.8,1), 
                       #space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

#Hierarchical clustering using pearsons
data(USArrests)
# Compute distances and hierarchical clustering
cor_matrix <- as.matrix(fpkm)
dd <- dist(scale(cor_matrix), method = "euclidean")
set.seed(100)
hc <- hclust(dd)
#edited
fpkm_matrix <- as.matrix(fpkm)
cor_matrix <- cor(fpkm_matrix)
dd <- as.dist(1 - cor_matrix)
set.seed(100)
hc <- hclust(dd, method = "ward.D2")
plot(hc)