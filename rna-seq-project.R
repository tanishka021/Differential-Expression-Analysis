#getting data
install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
install.packages("tidyverse")
BiocManager::install("Biobase")
BiocManager::install("limma")
library(Biobase)
library(limma)
library(tidyverse)
library(dplyr)
gset1 <- GEOquery::getGEO("GSE16561", AnnotGPL =TRUE)
eset1 <-gset1[[1]]
expr1 <- exprs(eset1)
pdata1<- pData(eset1)
fdata1 <- fData(eset1)
expr1[is.na(expr1)] <- rowMeans(expr1 , na.rm = TRUE)[row(expr1)[is.na(expr1)]]
#normalize the technical noices or samples should have same distribution
expr_norm <- normalizeBetweenArrays(expr1 ,  method = "quantile")
#group labels
group <- ifelse(grepl("Control", pdata1$description),"Control","Stroke")
group<- factor(group)
table(group)

#design matrix
design1 <- model.matrix(~0 + group)
design1
colnames(design1) <- levels(group)

#fit linear model and contrast
fit1 <- lmFit(expr_norm , design1)
contrast.matrix <- makeContrasts(Stroke_vs_Control = Stroke - Control, levels = design1)
fit2 <- contrasts.fit(fit1 , contrast.matrix)
fit2 <- eBayes(fit2)
fit2

#DE genes
deg_genes1 <- topTable(fit2 , coef = 1 , number = Inf , adjust.method = "BH")
head(deg_genes1)
nrow(deg_genes1)

#saving results
write.csv(deg_genes1, "results_of_limma_of_allgenes.csv" , row.names = TRUE)

str(deg_genes1)
#visualization
install.packages("ggplot2")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(deg_gene_level,
                lab= deg_gene_level$`Gene symbol`,
                x= "logFC",
                y= "adj.P.Val",
                pCutoff = 0.05,
                FCcutoff = 1,
                title="Volcano plot = stroke vs control")

library(EnhancedVolcano)
is.na(deg_genes1)
BiocManager::install("pheatmap")
library(pheatmap)
top50 <- rownames(deg_genes1)[1:50]
sum(is.na(expr_norm[top50,]))
annotation <- data.frame(Group = group)
rownames(annotation) <- colnames(expr_norm[top50,])   # must match sample names!

pheatmap(expr_50,
         scale = "row",
         annotation_col = annotation,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE)

save.image("rna-seq-project.RData")
load("rna-seq-project.RData")

top50 <- deg_gene_level %>%
  arrange(adj.P.Val)%>%
  head(50)
top50
expr_50 <- expr1[rownames(expr1) %in% top50$ID ,]
expr_50
annotation <- data.frame(Group = group)
rownames(annotation) <- colnames(expr_50)
#change probe ids with gene symbol

deg_genes1$ID <- rownames(deg_genes1)
fdata_new <- fdata1[,c("ID","Gene symbol")]
deg_annot <- merge(deg_genes1, fdata_new, by="ID")
head(deg_genes1)
deg_genes_clean <- deg_annot[! (is.na(deg_annot$`Gene symbol`) | deg_annot$`Gene symbol` ==""),]

#eliminate same gene symbols and keep only relevant one
deg_gene_level <- deg_genes_clean %>%
  group_by(`Gene symbol`) %>%
  slice_min(order_by = adj.P.Val, n=1) %>%
  ungroup()
nrow(deg_gene_level)

vars <- load("rna-seq-project.RData")
print(vars)
rm(annotated_expr,annotated_expr_clean,annotation,dat,dat_symbol)
rm(dat.long,deg_filtered,design,expr_df,expr_norm,expr_raw,expr_raw_clean,expr_symbol,feature,fit,gpl_annot)
rm(clean_names,cols_to_keep,group,group_list,ifc_cut,metadata_ids_clean,p_cut,sample_names,sample_names_clean,symbol_col,top_5000_genes)
rm(deg_results)
head(fdata1)
sum(is.na(deg_annot))
sum(deg_annot$`Gene symbol`=="")

#GO/KEGG - functional enrichment analysis
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library("clusterProfiler")
library("org.Hs.eg.db")

deg_strict <- subset(deg_gene_level ,P.Value < 0.05 & abs(logFC) > 0.5)
write.csv(deg_strict , "DEG_strict_654.csv")
gene_symbol <- deg_strict$`Gene symbol`
entrez_id <- bitr(gene_symbol,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
#go enrichment(enriched biological processes)
ego <- enrichGO(gene = entrez_id$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = "TRUE")
#kegg enrichment(enriched pathways)
ekegg <- enrichKEGG(gene = entrez_id$ENTREZID,
                    organism = "hsa",
                    pvalueCutoff = 0.05)
#visualization
dotplot(ego, showCategory=15 , split="ONTOLOGY") + facet_grid(ONTOLOGY ~ .,scale="free")
dotplot(ekegg, showCategory = 20)

#PPI NETWORK
deg_symbol <- deg_strict$`Gene symbol`

write.table(deg_symbol , "DEG_symbols_654.txt",
            quote=FALSE , row.names = FALSE , col.names = FALSE)

install.packages("pROC")
library("pROC")
top_genes <- deg_gene_level %>% arrange(adj.P.Val) %>% slice_head(n=10)
probe_ids <- top_genes$ID
probe_ids<- probe_ids[probe_ids %in% rownames(expr_norm)]
#function to compute auc per gene- deg can be checked statsistically , are they able to really predict - a person is stroked or not ?

gene_aucs <- lapply(seq_along(probe_ids), function(i){
  prb <- probe_ids[i]
  gene <- top_genes$`Gene symbol`[match(prb, top_genes$ID )]
  x= as.numeric(expr_norm[prb,])
  y=as.numeric(group =="Stroke")
  roc_obj <- roc(y,x, quiet=TRUE , direction = "<")
  data.frame(Gene = gene, ProbeID = prb, AUC = as.numeric(auc(roc_obj)))

})
gene_aurs <- bind_rows(gene_aucs) %>% arrange(desc(AUC))
gene_aurs
#build a simple 5 gene model
k<- 5
sel_probes <- gene_aurs$ProbeID[1:k]
x <- t(expr_norm[sel_probes,])
df <- data.frame(x,label=group)

#logistic regression
fit_new <- glm(label ~ .,data= df , family = binomial())
pred <- predict(fit_new,type = "response")
pred
roc_multi <- roc(df$label == "Stroke" , pred)
auc(roc_multi)
write.csv(gene_aurs,"roc_single_gene_auc.csv" , row.names=FALSE)

BiocManager::install("STRINGdb")
library("STRINGdb")

string_db <- STRINGdb$new(version = "12" , species = 9606, score_threshold = 700, input_directory = "")
mappped <-string_db$map(data.frame(gene = deg_symbol), "gene", removeUnmappedRows = TRUE)
hits <- string_db$get_interactions(mapped$string_db)




#roc with best genes
best_probeid <- gene_aurs$ProbeID[1:3]
best_gene <- gene_aurs$ProbeID[1]
best_name <- gene_aurs$Gene[1:3]
cols <- c("red","blue","green")

#empty plot
plot(NULL, xlim = c(1,0) , ylim =c(0,1),
     xlab = "False positive rate (1- specificity)",
     ylab = "true positive rate (Specificity)",
     main = "roc curves top 3 gene ")

#loop over top 3
for(i in 1:3){
  prb <- best_probeid[i]
  gene <- best_name[i]
  expr_vals <- as.numeric(expr_norm[best_gene, ])
  labels <- as.numeric(group == "Stroke")
  roc_obj <- roc(labels, expr_vals, quiet = TRUE, direction = "<")
  lines(roc_obj , col = cols[i], lwd = 2)
  auc_val <- round(auc(roc_obj), 3)
  best_name[i] <- paste0(gene, " (AUC=", auc_val, ")")
}

expr_vals <- as.numeric(expr1[best_gene, ])
labels <- as.numeric(group == "Stroke")

roc_obj <- roc(labels, expr_vals, quiet = TRUE, direction = "<")
plot.roc(roc_obj, main = paste("ROC curve -", best_name), col = "blue", lwd = 2)
par(mfrow=c(1,2))
for (prb in best_gene){
  gene <- gene_aucs$Gene[gene_aucs$ProbeID == prb]
  expr_vals <- as.numeric(expr1[prb,])
  labels <- as.numeric(group == "Stroke")
  
  roc_obj <- roc(labels , expr_vals , quiet = TRUE, direction = "<")
  plot(roc_obj , main = paste("ROC curve -", gene))
  }

length(group)
ncol(expr_norm)
length(expr_vals)
