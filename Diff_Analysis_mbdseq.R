# DS project 2024
# Differenatial Methylation Analysis
# Caroline Forsythe,Astghik Sarukhanyan, Julia Lapucha
# 11.06.2024

# Import libraries
library(edgeR)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(EnhancedVolcano)

setwd("D:/GIT/Master/DataScience/Project")

# Read the count matrix from a text file
count_matrix <- read.table("./output/Annotation/77Samples/annotedMethylationCounts_allsamples.csv", sep=",", header = T, row.names = 1)

# Check format
head(count_matrix)


# Read sample data 
sample_data <- read.table("./Data/all_Methylation/allMethylationSample_info.txt", header = TRUE, sep=",")
head(sample_data)

sample_data$Sample
dim(sample_data)

# Select only patient IDs that also have RNA-seq data from sample txt file
count_matrix <- count_matrix %>%
  dplyr::select(all_of(sample_data$Sample))

colnames(count_matrix) %>% length() # 77 

# check if order is the same 
sample_data <- sample_data[match(colnames(count_matrix), sample_data$Sample), ] 
all(sample_data$Sample == colnames(count_matrix))

dim(sample_data)

col_data <- data.frame(
  sample = colnames(count_matrix),
  conditions = sample_data$Group
)

######### Methylation Analysis using EdgeR


# Creat DGEList object
group <- col_data$conditions
y <- DGEList(counts = count_matrix, group = group)

# Filtering
keep <- rowSums(cpm(y) > 5) >= 2
y <- y[keep,]

dim(y)
head(cpm(y))


# Normalize 
y <- calcNormFactors(y, method = "TMM")


# Estimate dispersion
design <- model.matrix(~group)
y <- estimateDisp(y, design)
print("Dispersion estimates:")
print(y$common.dispersion)

# Fit the model and perform the statistical test
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)

# Results
results <- topTags(lrt, n = Inf, adjust.method = "BH")
print("Differential expression results:")
print(results)

# Results table
results_table <- as.data.frame(results$table)
print("Full results table:")
print(head(results_table))



# write.csv(results_table, "./output/DMR_allsamples.csv")

# map ensemblids to external gene name
geneList <- read.csv( "./output/Annotation/77Samples/biomart_output.csv")

# result table to geneList
results_table$ensembl_gene_id <- rownames(results_table)
rownames(results_table) <- 1:nrow(results_table)

externamNamesResults <- results_table%>%
  left_join(
    geneList,
    by = join_by(ensembl_gene_id),
    multiple = "all"
  )

externamNamesResults <-  externamNamesResults[!duplicated(externamNamesResults), ]

# Filter results with FDR <= 0.05
sig <- externamNamesResults %>%
  filter(FDR <= 0.05)

cat("Significant results (FDR <= 0.05):\n")
print(sig)

# save significant DMRs
saveRDS(sig, file = "./output/Annotation/77Samples/sig_DMR.rds")

####################################### volcano plot

alpha = 0.05
lfc_cutoff = 1

# significant DMR
sig_volc <- EnhancedVolcano(externamNamesResults,
                lab = externamNamesResults$external_gene_name,
                x = 'logFC',
                y = 'FDR',
                FCcutoff = lfc_cutoff, # cut-off for absolute log2 fold-change
                title = "Volcano Plot of significant DMRs",
                caption = paste("p < ", alpha, "and LFC Cutoff:", lfc_cutoff),
                pCutoff = alpha,
                drawConnectors = TRUE,
                widthConnectors = 1.5,
                pointSize=3)
# save plot
ggsave(plot=sig_volc,
       filename = paste0(getwd(), "output/sig_volc.png"),
       height = 5, width=9,dpi=320)


################################################ Boxplot gen expression

#filter sig genes for the 10 with biggest logfold 2 change
logfc<-sig[order(sig$logFC , decreasing = TRUE), ]
rownames(logfc) <- NULL
top10_logfc<- logfc$ensembl_gene_id[1:10]

count_matrix$ensembl_gene_id <- rownames(count_matrix)
filtered_counts<-as.data.frame( t(filter(count_matrix,
                        count_matrix$ensembl_gene_id %in% top10_logfc)))

filtered_counts$Sample<-row.names(filtered_counts)
top10_counts <- merge(filtered_counts, sample_data, by="Sample")

top10_counts[, 2:11] <- sapply(top10_counts[, 2:11], as.numeric)
top10_counts[, 2:11] <- log(top10_counts[, 2:11]) # try again

long_sig <- top10_counts %>% 
  pivot_longer(
    cols = "ENSG00000100083":"ENSG00000289901", 
    names_to = "ensembl_gene_id",
    values_to = "Count"
  )

boxplot_gene <- ggplot(long_sig, 
                       aes(x = reorder(ensembl_gene_id, Group), 
                           y = Count, 
                           color = Group)) +
  geom_boxplot() +
  geom_point(aes(color  = Group))+
  labs(title = "Methylation Distribution of log transformed counts (Top 10 genes; p.adj < 0.05)",
       x = "Ensemblids",
       y = "Log transformed counts") + 
  theme_gray()+
  theme(axis.text.x = element_text(face = "bold", 
                                   size = 10, angle = 60, hjust = 1)) +
  theme(axis.text.y = element_text(face = "bold", size =10))

# saving plot
ggsave(plot=boxplot_gene,
       filename = paste0(getwd(), "output/boxplot_top10.png"),
       height = 5, width=9,dpi=320)


############################################# preparing matrix for ML
count_matrix$ensembl_gene_id <- rownames(count_matrix)

count_matrix_sig <- filter(count_matrix,
                       ensembl_gene_id %in% sig$ensembl_gene_id)
count_matrix_sig<-count_matrix_sig  %>% dplyr::select(-ensembl_gene_id)


# transpose counts
transposed_counts <- t(count_matrix_sig)
write.csv(transposed_counts, "./output/Annotation/77Samples/transposedMethylationCounts.csv")


