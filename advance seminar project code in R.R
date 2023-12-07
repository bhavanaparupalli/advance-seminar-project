# Read the table from the specified file
file_path <- 'C:/Users/bhava/Downloads/GSE246995_counts (1).txt.gz'
counts1.a <- read.table(file_path, sep="\t", header=FALSE, stringsAsFactors=FALSE)
# Display the data using the View function
View(counts1.a)
# Extract the data from the read table, excluding the first row and first column
counts1 <- as.matrix(counts1.a[-1, -1])
# Convert the matrix to integer mode
mode(counts1) <- "integer"
# Set rownames to geneIDs (excluding the first column)
rownames(counts1) <- counts1.a[-1, 1]
# Set colnames to the values from the first row (excluding the first element)
colnames(counts1) <- counts1.a[1, -1]
#Display the number of columns in the counts1 matrix
num_cols <- ncol(counts1)
cat("Number of columns:", num_cols, "\n")
# Display the first few rows of the counts1 matrix
head(counts1)

# Read the CSV file with column metadata
colData_file_path <- 'C:/Users/bhava/Downloads/Book 2 (2).csv'
colData_original <- read.csv(colData_file_path)
# Display the original data using the View function
View(colData_original)
# Extract the data from the read CSV, excluding the first row and first column
colData_matrix <- as.matrix(colData_original[-1, -1])
# Uncomment the line below if you want to set the mode to "integer"
# mode(colData_matrix) <- "integer"
# Set rownames to the values in the first column (excluding the first element)
rownames(colData_matrix) <- colData_original[-1, 1]
# Set colnames to the values from the first row (excluding the first element)
colnames(colData_matrix) <- colData_original[1, -1]
# Display the resulting matrix using the View function
View(colData_matrix)

# Make sure the column names in 'counts2' and 'colData2' are identical and in the correct format
all(colnames(counts1) %in% rownames(colData1))
all(colnames(counts1) == rownames(colData1)) 

# Create a DESeqDataSet object
dds1.0 <- DESeqDataSetFromMatrix(countData = counts1, 
                                 colData = colData1, 
                                 design = ~ Condition)

# Set the filtering criteria: Keep genes with at least 10 counts
keep <- rowSums(counts(dds1.0)) >= 10
# Filter the original dataset based on the criteria
dds1 <- dds1.0[keep,]
# Display the filtered dataset
dds1
# Set the factor level: Change the reference level of the 'Condition' variable to "control"
# Note: Assuming 'Condition' is a factor variable in dds1.0
dds1$condition <- relevel(dds1.0$Condition, ref = "control")
# Display the levels of the 'Condition' variable in the original dataset (dds1.0)
levels(dds1.0$Condition)
# At this point, 'dds1' contains the filtered and re-leveled data
# Assuming 'dds_filtered' is your filtered and re-leveled dataset
dds_filtered <- DESeq(dds1)
# Get DESeq results
res <- results(dds_filtered)
res
summary(res)

# Additional analysis with a significance level of 0.05
res1 <- results(dds_filtered, alpha = 0.05)
res1
summary(res1)

# Plotting MA plot for visualization
plotMA(res1)
# Decrease margins and adjust plotting area
par(mar = c(1, 3, 1, 1) + 0.1)  # Adjust the margins as needed
plotDispEsts(dds_filtered, cex = 0.5)
# Assuming 'deg_results' is the data frame containing DEGs from your analysis
deg_results <- as.data.frame(subset(res1, (abs(res1$log2FoldChange) >= 1) & (res1$padj <= 0.05)))
# Display the number of DEGs
num_DEGs <- nrow(deg_results)
print(num_DEGs)

# Write the DEG results to a CSV file
write.csv(deg_results, file = "DEG_results.csv")

# View the CSV file (open it in a separate window)
View(read.csv("DEG_results.csv"))

# Set threshold values
padj.cutoff <- 0.05
lfc.cutoff <- 2
# Create a threshold for significant results
threshold <- res1$padj < padj.cutoff & abs(res1$log2FoldChange) > lfc.cutoff
# Display the number of significant results based on the threshold
num_significant <- length(which(threshold))
print(num_significant)
# Add the threshold information to the results dataframe
res1$threshold <- threshold
# Create a dataframe with only significant results
sig_results <- data.frame(subset(res1, threshold = TRUE))
# Order significant results by adjusted p-value
sig_ordered <- sig_results[order(sig_results$padj), ]
# Select the top 20 significant results
top20_sig <- rownames(sig_ordered[1:20, ])
# Get normalized counts for the top 20 significant genes
normalized_counts <- counts(dds_filtered, normalized = TRUE)
top_20_sig_norm <- as.data.frame(normalized_counts[top20_sig, ])
# Add gene names as a separate column
top_20_sig_norm$gene <- rownames(top_20_sig_norm)
# View the top 20 significant genes with normalized counts
View(top_20_sig_norm)

# Load required libraries (if not already installed)
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("DESeq2")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(DESeq2)
library(pheatmap)
# Assuming 'dds1' is your DESeqDataSet
# Replace 'dds1' with your actual DESeqDataSet
dds1 <- DESeqDataSetFromMatrix(countData = counts1, colData = colData1, design = ~ Condition)
# Assuming 'dds_filtered' is your filtered and re-leveled dataset
dds_filtered <- DESeq(dds1)
# Get DESeq results
res <- results(dds_filtered)
# Additional analysis with a significance level of 0.05
res1 <- results(dds_filtered, alpha = 0.05)
# Extract normalized counts for significant genes
sig_genes <- rownames(subset(res1, padj < 0.05 & abs(log2FoldChange) > 2))
normalized_counts_sig <- counts(dds_filtered, normalized = TRUE)[sig_genes, ]

# Create a heatmap
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(log2(normalized_counts_sig + 1), cluster_cols = FALSE, col = heatmap_colors,
         main = 'Heatmap of Normalized Counts for DEGs', fontsize = 8)
# Assuming 'res1' is the results dataframe containing log2FoldChange and padj
# Convert the 'padj' column to numeric
res1$padj <- as.numeric(res1$padj)

# Assuming 'res1' is the results dataframe containing log2FoldChange and padj
# Convert the 'padj' column to numeric
res1$padj <- as.numeric(res1$padj)

# Load the ggplot2 library
library(ggplot2)

# Create a Volcano Plot
ggplot(res1, aes(x = log2FoldChange, y = -log10(padj), color = ifelse(padj < 0.05, 'red', 'violet'))) +
  geom_point() +
  ggtitle('Volcano Plot') +
  xlab('Log2 Fold Change') +
  ylab('-log10(Adjusted P-value)') +
  theme_minimal()
