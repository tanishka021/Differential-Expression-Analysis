if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("limma")
library(GEOquery)
library(limma)
install.packages("ggplot2")
install.packages("ggrepel")
gse <- getGEO("GSE16561", GSEMatrix = TRUE)
gse <- gse[[1]]
fdata <- fData(gse)
head(fdata)
gene_symbols<- fdata$ILMN_Gene
# Expression matrix
expr_raw <- exprs(gse)
# Check range
summary(expr_raw[,1:5])
sum(is.na(expr_raw))
expr1 <- expr_raw
expr1[is.na(expr1)] <- rowMeans(expr1 , na.rm = TRUE)[row(expr1)[is.na(expr1)]]
sum(is.na(expr1))
exprSet<- expr1
# Step 2: KNN imputation
library(impute)
exprSet <- impute.knn(exprSet)$data
sum(is.na(exprSet))
exprSet <- cbind(GeneSymbol = gene_symbols, exprSet)
exprSet <- as.data.frame(exprSet)
exprSet <- exprSet[!is.na(exprSet$GeneSymbol),]
sum(is.na(exprSet))
exprSet[,-1] <- lapply(exprSet[,-1], as.numeric)
str(exprSet)
library(dplyr)
exprSet <- exprSet %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), mean))
exprSet <- as.data.frame(exprSet)
rownames(exprSet) <- exprSet$GeneSymbol
exprSet <- exprSet[, -1]
str(exprSet)
dim(exprSet)
head(rownames(exprSet))
genes_to_check <- c("CD163", "LAMP2", "PICALM", "RGS2","PIN1")
genes_to_check %in% rownames(exprSet)
# Metadata
meta <- pData(gse)
head(meta)
# Extract group from title
group_list <- ifelse(grepl("Stroke", meta$title), "Stroke", "Control")
# Convert to factor
group_list <- factor(group_list, levels = c("Control", "Stroke"))
# Check counts
save.image("stroke-script.RData")
table(group_list)
str(group_list)
levels(group_list)
design <- model.matrix(~0 + group_list)
colnames(design) <- levels(group_list)
design
library(limma)
fit <- lmFit(exprSet, design)
contrast.matrix <- makeContrasts(Stroke - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, adjust = "fdr", number = Inf)
head(deg)
genes_to_check %in% rownames(deg)
deg[genes_to_check,]
deg_filtered <- deg[abs(deg$logFC) > 0.5 & deg$adj.P.Val < 0.05, ]
#acc to ppr - 482 , mine: 525
genes_to_check %in% rownames(deg_filtered)
library(ggplot2)
library(ggrepel)
# Add significance column
deg$Significance <- "Not Significant"
deg$Significance[deg$logFC > 0.5 & deg$adj.P.Val < 0.05] <- "Upregulated"
deg$Significance[deg$logFC < -0.5 & deg$adj.P.Val < 0.05] <- "Downregulated"
ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significance), alpha = 0.6) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  geom_point(data = deg[rownames(deg) %in% genes_to_check, ],
             color = "black", size = 2) +
  geom_text_repel(
    data = deg[rownames(deg) %in% genes_to_check, ],
    aes(label = rownames(deg)[rownames(deg) %in% genes_to_check]),
    size = 4,
    max.overlaps = 10
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal()
sum(deg$Significance %in% c("Upregulated", "Downregulated"))
up <- sum(deg$Significance == "Upregulated")
down <- sum(deg$Significance == "Downregulated")
print(up)
print(down)