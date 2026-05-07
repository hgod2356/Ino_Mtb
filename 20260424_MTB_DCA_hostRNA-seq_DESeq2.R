if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")

#BiocManager::install(c("tximport", "GenomicFeatures", "tximeta" , force = TRUE))
#BiocManager::install("BiocParallel", force = TRUE )
#BiocManager::install("txdbmaker" , force = TRUE)

setwd("D:/1.Experiment/충남대/MTBDCA_Deseq2")

library(tximport)
library(GenomicFeatures)
library(txdbmaker)

txdb2 <- makeTxDbFromGFF(file="gencode.vM38.annotation.gff3.gz")
transcripts(txdb2, columns=c("tx_id", "tx_name"))
genes(txdb2)
txdb2
columns(txdb2)

#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(tximeta)
columns(org.Mm.eg.db)

### extract gene_name in '.gff3'
library(rtracklayer)
#aa <- import.gff3("Mus_musculus.GRCm39.111.gff3.gz") # 25s
aa2 <- import.gff3("gencode.vM38.annotation.gff3.gz")
gene.DF <- data.frame(aa2$gene_name,aa2$gene_id)
u.gene.DF<- unique(gene.DF)
head(u.gene.DF)

u.gene.DF[,2] <- gsub("\\..*", "", u.gene.DF[,2])
head(u.gene.DF)

### 
k <- keys(txdb2, keytype = "TXNAME")
tx2gene <- select(txdb2, k, "GENEID", "TXNAME")
head(tx2gene)

#tx2gene에서 dot remove하기
tx2gene[,1] <- gsub("\\..*", "", tx2gene[,1])
tx2gene[,2] <- gsub("\\..*", "", tx2gene[,2])

head(tx2gene)

#write.csv(as.data.frame(u.gene.DF), file = "gene_name_20251106.csv")
#write.csv(as.data.frame(tx2gene), file = "gene_name_tx2gene_20251106.csv")
#write.csv(as.data.frame(txi.salmon$counts), file = "gene_count_20251106.csv")

###
setwd("D:/1.Experiment/충남대/MTBDCA_Deseq2/salmon")
dir <- getwd()
files = paste0(dir, "/",c("Con1","Con2","Con3","Exp1","Exp2","Exp3"), "/quant.sf")

names(files) <- paste0( c("Con1","Con2","Con3","Exp1","Exp2","Exp3"))
file.exists(files)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

head(txi.salmon)

###deseq2 - using counts
#BiocManager::install("DESeq2")
setwd("D:/1.Experiment/충남대/MTBDCA_Deseq2")
library(DESeq2)
sampleTable <- data.frame(
  condition = factor(c(rep("Mtb",3), rep("Mtb_DCA",3)))
)
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)
deseqresult <- DESeq(dds)
deseqresult

#condition1 Mtb vs Mtb_DCA
res01 <- results(deseqresult, contrast = c("condition","Mtb_DCA","Mtb"), alpha = 0.05)
res01
summary(res01)
plotMA(res01, ylim=c(-10,10))
#write.csv(as.data.frame(res01), file = "1.MtbvsMtbDCA_results_count_20251106.csv")

# Genename mapping
res01.df <- as.data.frame(res01)
res01.df$gene_id <- rownames(res01.df)
gene_map <- u.gene.DF
colnames(gene_map) <- c("genename", "gene_id")
merged <- merge(res01.df, gene_map, by = "gene_id", all.x = TRUE)
rownames(merged) <- merged$gene_id
merged2<-merged[,-1]
dups <- merged$gene_name[!is.na(merged$gene_name)]
dups <- dups[duplicated(dups)]
sort(unique(dups))
#write.csv(as.data.frame(merged2), file = "1.MtbvsMtbDCA_results_count_genename_20251106.csv")


##ggplot -volcano plot
library(ggplot2)
library(ggrepel)
res01genename <- read.csv("1.MtbvsMtbDCA_results_count_genename_20251106.csv", stringsAsFactors = F)
topT.res01 <- as.data.frame(res01genename)
topT.res01 <- topT.res01[!is.na(topT.res01$padj), ]

cut_lfc <- 0.5
cut_pvalue <- 0.05

topT.res01$group <- "Not significant"
topT.res01$group[topT.res01$padj < cut_pvalue & topT.res01$log2FoldChange >  cut_lfc] <- "Mtb"
topT.res01$group[topT.res01$padj < cut_pvalue & topT.res01$log2FoldChange < -cut_lfc] <- "Mtb+DCA"

volcano_colors <- c("Mtb"="#E63943", 
                    "Mtb+DCA"="#2A9D8F", 
                    "Not significant"="grey80")
exclude_patterns <- "^Gm[0-9]|Rik$|-ps|^ENSMUSG"
label_data <- subset(topT.res01, 
                     group != "Not significant" & 
                       (abs(log2FoldChange) >= 5 | -log10(padj) >= 25) &
                       !grepl(exclude_patterns, genename))

# ggplot
ggplot(topT.res01, aes(x=log2FoldChange, y=-log10(padj), color=group)) +
  geom_point(alpha=0.8, size=1.8) +
  scale_color_manual(values=volcano_colors) +
  geom_text_repel(data = label_data, 
                  aes(label = genename),
                  size = 3.5,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  force = 2,
                  max.overlaps = Inf,
                  show.legend = FALSE) +
  geom_vline(xintercept = c(-cut_lfc, 0, cut_lfc), 
             linetype = c("dashed","dotted","dashed"), color="black") +
  geom_hline(yintercept = -log10(cut_pvalue), linetype="dashed", color="black") +
  labs(title="Volcano plot: Mtb vs Mtb+DCA",
       x=expression(Log[2]~fold~change),
       y=expression(-Log[10]~Q~value),
       color="Group") +
  theme_bw(base_size=14) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold"),
    legend.position = "right",
    legend.box.background = element_rect(color="black", fill="white", size=1.0),
    legend.key = element_rect(fill="white", color="white")
  )

#interest gene labeling
target_genes <- c("Pparg", "Gdf15", "Fabp3", "Acsl1", "Acsl5", 
                  "Ppard", "Auh", "Cpt2", "Nr1h3", "Ldlr", "Sik1", "Il6", "Cxcl1", "Cxcl3", "Ifng", "Prkaa2", "Prkab2", "Abhd8", "Muc2", "Muc19")

label_data <- subset(topT.res01, genename %in% target_genes)

volcano_colors <- c("Mtb"="#E63943", 
                    "Mtb+DCA"="#2A9D8F", 
                    "Not significant"="grey80")

ggplot(topT.res01, aes(x=log2FoldChange, y=-log10(padj), color=group)) +
  geom_point(alpha=0.6, size=1.8) +
  scale_color_manual(values=volcano_colors) +
  geom_text_repel(data = label_data, 
                  aes(label = genename),
                  size = 4,
                  fontface = "plain",
                  color = "black",
                  box.padding = 0.8,
                  point.padding = 0.5,
                  force = 3,
                  segment.color = "grey30",
                  max.overlaps = Inf,
                  show.legend = FALSE) +
  geom_vline(xintercept = c(-cut_lfc, 0, cut_lfc), 
             linetype = c("dashed","dotted","dashed"), color="black") +
  geom_hline(yintercept = -log10(cut_pvalue), linetype="dashed", color="black") +
  labs(title="Volcano plot: Mtb (Left) vs Mtb+DCA (Right)",
       subtitle="Labeled: Selected Lipid Metabolism Genes",
       x=expression(Log[2]~fold~change),
       y=expression(-Log[10]~Q~value),
       color="Group") +
  theme_bw(base_size=14) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size=11, color="grey30"),
    legend.position = "right",
    legend.box.background = element_rect(color="black", fill="white", size=0.8)
  )

#Extracting transformed values - vst, rlog, ntd
vsd <- vst(deseqresult, blind = "TRUE")
rld <- rlog(deseqresult, blind = TRUE)
ntd <- normTransform(deseqresult)
head(assay(vsd), 3)
head(assay(rld), 3)
head(assay(ntd), 3)

matched_symbols <- u.gene.DF$aa2.gene_name[match(rownames(vsd), u.gene.DF$aa2.gene_id)]
matched_symbols[is.na(matched_symbols)] <- rownames(vsd)[is.na(matched_symbols)]
final_names <- make.unique(matched_symbols)
rownames(vsd) <- final_names
head(rownames(vsd))

#ggplot - PCA
vsdpcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
vsdpcaData$group <- factor(vsdpcaData$group, levels = c("Mtb","Mtb_DCA"))
ggplot(vsdpcaData, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", round(attr(vsdpcaData, "percentVar")[1] * 100), "% variance"),
       y = paste0("PC2: ", round(attr(vsdpcaData, "percentVar")[2] * 100), "% variance"),
       color = "Group") +
  scale_color_manual(values = c("Mtb"="#E63943", "Mtb_DCA"="#2A9D8F")) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


#Pathway analysis ---------------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(stringr)
library(enrichplot)
library(dplyr)

up_genes <- topT.res01$genename[topT.res01$log2FoldChange > 0.25 & topT.res01$padj < 0.05]
up_genes <- unique(na.omit(up_genes))

down_genes <- topT.res01$genename[topT.res01$log2FoldChange < -0.25 & topT.res01$padj < 0.05] # 오타 주의: padj < 0.05
down_genes <- unique(na.omit(down_genes))

gene_list <- list(Mtb_DCA = up_genes, Mtb = down_genes)

# GO Enrichment (Biological Process)
ck_go <- compareCluster(geneCluster = gene_list, 
                        fun         = "enrichGO", 
                        OrgDb       = org.Mm.eg.db, 
                        keyType     = "SYMBOL", 
                        ont         = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)

ck_go_df <- as.data.frame(ck_go)
pathways_mtb <- ck_go_df %>% filter(Cluster == "Mtb") %>% pull(Description) %>% unique()
pathways_dca <- ck_go_df %>% filter(Cluster == "Mtb_DCA") %>% pull(Description) %>% unique()

# 3.unique pathway extraction
only_mtb_list <- setdiff(pathways_mtb, pathways_dca)
only_dca_list <- setdiff(pathways_dca, pathways_mtb)

only_mtb_df <- ck_go_df %>% 
  filter(Cluster == "Mtb" & Description %in% only_mtb_list) %>%
  arrange(p.adjust)

only_dca_df <- ck_go_df %>% 
  filter(Cluster == "Mtb_DCA" & Description %in% only_dca_list) %>%
  arrange(p.adjust)

# 5. CSV save
#write.csv(only_mtb_df, "Pathways_Only_in_MTB.csv", row.names = FALSE)
#write.csv(only_dca_df, "Pathways_Only_in_MTB_DCA.csv", row.names = FALSE)
#write.csv(ck_go_df, file = "GO_comparison_Mtb_vs_MtbDCA_20251106.csv", row.names = FALSE)

# 1. Shared Pathway
common_pathways_list <- intersect(pathways_mtb, pathways_dca)
cat("Mtb Mtb_DCA shared Pathway:", length(common_pathways_list), "\n")

common_mtb_info <- ck_go_df %>%
  filter(Cluster == "Mtb" & Description %in% common_pathways_list) %>%
  select(Description, GeneRatio, p.adjust, Count, zScore, FoldEnrichment)

common_dca_info <- ck_go_df %>%
  filter(Cluster == "Mtb_DCA" & Description %in% common_pathways_list) %>%
  select(Description, GeneRatio, p.adjust, Count, zScore, FoldEnrichment)

common_comparison_df <- full_join(common_mtb_info, common_dca_info, 
                                  by = "Description", 
                                  suffix = c("_Mtb", "_DCA")) %>%
  mutate(zScore_Diff = zScore_DCA - zScore_Mtb) %>%
  arrange(p.adjust_DCA)

# write.csv(common_comparison_df, "Common_Pathways_Comparison.csv", row.names = FALSE)

# target pathway visualization
target_pathways <- c(
  "lipid storage",
  "lipid droplet formation",
  "positive regulation of lipid biosynthetic process",
  "glycerolipid biosynthetic process",
  "cholesterol biosynthetic process",
  "cholesterol storage",
  "triglyceride biosynthetic process",
  "foam cell differentiation",
  "macrophage derived foam cell differentiation",
  "long-chain fatty acid transport",
  "lipid oxidation"
)

# dotplot
p_lipid <- dotplot(ck_go, showCategory = target_pathways) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
  labs(title = "Lipid Metabolism & Foam Cell Differentiation",
       subtitle = "Mtb vs Mtb+DCA",
       x = "Group",
       y = "GO Term Description") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_lipid

#bar plot
plot_data <- as.data.frame(ck_go) %>%
  filter(Description %in% target_pathways) %>%
  mutate(logP = -log10(p.adjust)) %>%
  mutate(Description = str_wrap(Description, width = 35)) %>%
  arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description))

max_logp <- max(plot_data$logP)
max_count <- max(plot_data$Count)
scale_factor <- max_count / max_logp

ggplot(plot_data, aes(x = Description)) +
  geom_bar(aes(y = logP), stat = "identity", fill = "black", width = 0.7) +
  geom_line(aes(y = Count / scale_factor, group = 1, color = "No. of genes"), 
            size = 1) +
  geom_point(aes(y = Count / scale_factor, color = "No. of genes"), 
             size = 2, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_manual(name = NULL, values = c("No. of genes" = "#F4A261")) +
  scale_y_continuous(
    name = expression(bold(-log[10](p.adjust))), 
    expand = expansion(mult = c(0, 0.1)),
    sec.axis = sec_axis(~ . * scale_factor
                        #,
                        #name = "No. of genes"
                        )
  ) +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(x = NULL, 
       title = "Lipid Metabolism & Foam Cell Pathway",
       subtitle = "(H37Rv only enriched pathway, p.adjust < 0.05)", 
       ) +
  theme(
    legend.position = c(0.8, 0.1),      
    legend.background = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_text(color = "black", size = 10, face = "bold"),
    axis.line.x.top = element_line(color = "black"),
    axis.text.x.top = element_text(color = "black", face = "bold"),
    axis.title.x.top = element_text(color = "black", face = "bold")
  )

#write.csv(plot_data, "Lipid_Metabolism_Barplot_Data.csv", row.names = FALSE)

# Inflammation related pathways
immune_target_pathways <- c(
  "leukocyte homeostasis",
  "chemotaxis",
  "leukocyte migration",
  "interleukin-6 production",
  "T cell homeostasis",
  "cell chemotaxis",
  "regulation of T cell activation",
  "positive regulation of innate immune response",
  "myeloid leukocyte activation",
  "pattern recognition receptor signaling pathway",
  "tumor necrosis factor production",
  "leukocyte chemotaxis",
  "macrophage chemotaxis"
)

ck_go@compareClusterResult$Cluster <- factor(ck_go@compareClusterResult$Cluster, 
                                             levels = c("Mtb", "Mtb_DCA"))

p_immune <- dotplot(ck_go, showCategory = immune_target_pathways) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
  labs(title = "Inflammation related pathways",
       subtitle = "Mtb vs Mtb+DCA",
       x = "Group",
       y = "GO Term Description") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_immune

immune_dotplot_data <- ck_go@compareClusterResult %>%
  filter(Description %in% immune_target_pathways)
head(immune_dotplot_data)
#write.csv(immune_dotplot_data, "Inflammation_related_pathways_dotplot_raw.csv", row.names = FALSE)

# heatmap------------
library(ComplexHeatmap)
library(RColorBrewer)

target_pathways <- c(
  "lipid storage", 
  "positive regulation of lipid biosynthetic process",
  "glycerolipid biosynthetic process", 
  "cholesterol biosynthetic process",
  "triglyceride biosynthetic process", 
  "macrophage derived foam cell differentiation",
  "long-chain fatty acid transport",
  "lipid droplet formation", 
  "cholesterol storage", 
  "foam cell differentiation",
  "lipid oxidation"
)

ck_go_df <- as.data.frame(ck_go)
target_pathway_data <- ck_go_df[ck_go_df$Description %in% target_pathways, ]
all_genes <- unique(unlist(strsplit(target_pathway_data$geneID, "/")))

gene_anno_mat <- matrix("0", nrow = length(all_genes), ncol = length(target_pathways))
rownames(gene_anno_mat) <- all_genes
colnames(gene_anno_mat) <- target_pathways

for(i in 1:nrow(target_pathway_data)) {
  pw_desc <- target_pathway_data$Description[i]
  genes <- unlist(strsplit(target_pathway_data$geneID[i], "/"))
  gene_anno_mat[genes, pw_desc] <- "1"
}

mat <- assay(vsd)[rownames(vsd) %in% all_genes, ]
mat_scaled <- t(scale(t(mat))) 
mat_scaled_horiz <- t(mat_scaled) 

gene_anno_mat_ordered <- gene_anno_mat[colnames(mat_scaled_horiz), , drop=FALSE]

anno_colors_vector <- brewer.pal(min(12, length(target_pathways)), "Set3")
anno_col_list <- list()
for(i in seq_along(target_pathways)) {
  pw <- target_pathways[i]
  anno_col_list[[pw]] <- c("0" = "white", "1" = anno_colors_vector[i])
}

bottom_ha_desc = columnAnnotation(
  df = as.data.frame(gene_anno_mat_ordered),
  col = anno_col_list,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_name_rot = 0,              
  annotation_name_side = "right",
  show_legend = FALSE,
  simple_anno_size = unit(3, "mm")
)

col_ann_df <- as.data.frame(colData(vsd)[, "condition", drop=FALSE])
row_ha_samples = rowAnnotation(
  condition = col_ann_df$condition,
  col = list(condition = c("Mtb" = "#E63943", "Mtb_DCA" = "#2A9D8F")),
  show_annotation_name = FALSE
)

new_sample_names <- c("MTB_1", "MTB_2", "MTB_3", "MTB+DCA_1", "MTB+DCA_2", "MTB+DCA_3")

hp_horiz <- Heatmap(mat_scaled_horiz, 
                    name = "Z-score", 
                    col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                    column_title = "Lipid accumulation & Foam cell associated gene",
                    row_labels = new_sample_names,
                    left_annotation = row_ha_samples, 
                    bottom_annotation = bottom_ha_desc,
                    show_column_names = TRUE, 
                    column_names_gp = gpar(fontsize = 8), 
                    row_names_gp = gpar(fontsize = 10, face = "bold"),
                    column_names_rot = 45,
                    cluster_rows = TRUE, 
                    cluster_columns = TRUE,
                    row_dend_side = "left",
                    column_dend_side = "top")

draw(hp_horiz, merge_legend = TRUE)


#lipid accumulation, foam cell differentiation gene volcano plot
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

target_pathways <- c(
  "lipid storage", 
  "glycerolipid biosynthetic process", "cholesterol biosynthetic process",
  "triglyceride biosynthetic process", "macrophage derived foam cell differentiation",
  "long-chain fatty acid transport", "lipid droplet formation", 
  "cholesterol storage", "foam cell differentiation"
)

target_genes_from_go <- as.data.frame(ck_go) %>%
  filter(Description %in% target_pathways) %>%
  pull(geneID) %>%
  strsplit("/") %>% 
  unlist() %>% 
  unique()

label_data <- subset(topT.res01, 
                     genename %in% target_genes_from_go & 
                       group != "Not significant")

ggplot(topT.res01, aes(x=log2FoldChange, y=-log10(padj), color=group)) +
  geom_point(alpha=0.4, size=6) + 
  geom_point(data = label_data, aes(x=log2FoldChange, y=-log10(padj)), 
             alpha=1, size=6, shape=21, fill=NA, stroke=0.8, color="black") +
  scale_color_manual(values=volcano_colors) +
  geom_text_repel(data = label_data, 
                  aes(label = genename),
                  size = 4.5,            
                  fontface = "bold",    
                  color = "black",
                  box.padding = 0.6, 
                  point.padding = 0.4,
                  force = 5,             
                  segment.color = "grey20",
                  max.overlaps = 50,    
                  show.legend = FALSE) +
  geom_vline(xintercept = c(-cut_lfc, 0, cut_lfc), 
             linetype = c("dashed","dotted","dashed"), color="black") +
  geom_hline(yintercept = -log10(cut_pvalue), linetype="dashed", color="black") +
  labs(title="Volcano plot: Lipid & Foam Cell Related Genes",
       #subtitle=paste0("Labeled: Genes in ", length(target_pathways), " Lipid-related GO Terms"),
       x=expression(Log[2]~fold~change),
       y=expression(-Log[10]~Q~value),
       color="Group") +
  theme_bw(base_size=15) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold", size=18),
    plot.subtitle = element_text(hjust=0.5, size=12, color="grey30"),
    legend.position = "right",
    legend.title = element_text(face="bold"),
    legend.box.background = element_rect(color="black", fill="white", size=0.8)
  )
