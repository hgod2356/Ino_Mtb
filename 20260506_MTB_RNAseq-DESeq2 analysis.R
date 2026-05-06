setwd("D:/1.Experiment/충남대/MTB RNA-seq")

#package load
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(clusterProfiler)
library(tidyr)
library(pheatmap)
library(grid)

gtf <- import("Reference/GCF_000195955.2_ASM19595v2_genomic.gtf")

gtf_df <- as.data.frame(gtf)
#write.csv(gtf.df, file="MTB_gtf_table.csv")

# 1. gene type - extract locus_tag, gene symbol
genes <- gtf_df %>%
  filter(type == "gene") %>%
  select(locus_tag, gene) %>%
  distinct()

# 2. CDS type - extract locus_tag, product 
products <- gtf_df %>%
  filter(type == "CDS") %>%
  select(locus_tag, product) %>%
  distinct()

# 3. merge table (locus_tag)
gene_map <- left_join(genes, products, by = "locus_tag")

# results check
head(gene_map)
#write.csv(gene_map, file="gene_map.csv")

# loading STAR output files
files <- c("Aligned_STAR/Bacteria_1_ReadsPerGene.out.tab",
           "Aligned_STAR/Bacteria_2_ReadsPerGene.out.tab",
           "Aligned_STAR/Combi_1_ReadsPerGene.out.tab",
           "Aligned_STAR/Combi_2_ReadsPerGene.out.tab")

read_counts <- function(file_path) {
  dat <- read.table(file_path, skip = 4)
  return(dat[, c(1, 3)]) # 3열(Forward) 선택
}

counts_list <- lapply(files, read_counts)

count_matrix <- do.call(cbind, lapply(counts_list, function(x) x[, 2]))
rownames(count_matrix) <- counts_list[[1]][, 1]
colnames(count_matrix) <- c("Bacteria_1", "Bacteria_2", "Combi_1", "Combi_2")

head(count_matrix)

# metadata fixing and DESeq2
col_data <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(c("MTB", "MTB", "MTB_DCA", "MTB_DCA"), 
                     levels = c("MTB", "MTB_DCA")) 
)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ condition)

# 3. Low count filtering
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 4. DESeq2
dds <- DESeq(dds)

# 5. Results check
res <- results(dds, contrast = c("condition", "MTB_DCA", "MTB"), alpha = 0.05)
summary(res)
res_df <- as.data.frame(res[order(res$padj), ])

# gene name mapping
res_df$locus_tag <- rownames(res_df)
res_final <- left_join(res_df, gene_map, by = "locus_tag")

res_final <- res_final %>%
  mutate(gene = ifelse(is.na(gene), locus_tag, gene))%>%
  # 4. locus_tag 컬럼을 맨 앞으로 이동
  relocate(locus_tag, .before = 1)

#write.csv(res_final, "DESeq2_Results_MTB_vs_MTB_DCA.csv", row.names = FALSE)

#value transformation
vsd <- vst(dds, blind = "TRUE")
rld <- rlog(dds, blind = TRUE)
ntd <- normTransform(dds)

# PCA - ggplot2 visualization
vsdpcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
vsdpcaData$group <- factor(vsdpcaData$group, levels = c("MTB","MTB_DCA"))

ggplot(vsdpcaData, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", round(attr(vsdpcaData, "percentVar")[1] * 100), "% variance"),
       y = paste0("PC2: ", round(attr(vsdpcaData, "percentVar")[2] * 100), "% variance"),
       color = "Group") +
  scale_color_manual(values = c("MTB"="#E63943", "MTB_DCA"="#2A9D8F")) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# volcano plot
# 1. grouping column creation
res_final <- res_final %>%
  mutate(change = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Mtb+DCA",
    padj < 0.05 & log2FoldChange < -1 ~ "Mtb",
    TRUE ~ "NS"
  ))

# 2. data filtering 
label_data <- res_final %>% 
  filter(change != "NS") %>%                       # 유의미한 유전자 중
  filter(!is.na(gene) & gene != locus_tag)        # 진짜 이름(Symbol)이 있는 경우만

# 3. Volcano Plot
ggplot(res_final, aes(x = log2FoldChange, y = -log10(padj))) +
  # 배경 점 (유의미하지 않은 유전자)
  geom_point(data = filter(res_final, change == "NS"), 
             color = "grey80", alpha = 0.4, size = 1.5) +
  geom_point(data = filter(res_final, change == "Mtb"), 
             aes(color = "Mtb"), alpha = 0.8, size = 2) +
  geom_point(data = filter(res_final, change == "Mtb+DCA"), 
             aes(color = "Mtb+DCA"), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Mtb" = "#E63943", "Mtb+DCA" = "#2A9D8F")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,    
    box.padding = 0.5,     
    point.padding = 0.3,   
    segment.color = "grey50", 
    force = 2            
  ) +
  theme_bw() +
  labs(title = "Differential Expression: MTB vs MTB+DCA",
       x = "log2 Fold Change", 
       y = "-log10 Adjusted P-value",
       color = "Group") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


# pathway analysis
# 1. gene list extraction
up_genes <- res_final %>%
  filter(padj < 0.05 & log2FoldChange > 0.5) %>%
  pull(locus_tag)

down_genes <- res_final %>%
  filter(padj < 0.05 & log2FoldChange < -0.5) %>%
  pull(locus_tag)

# 2. KEGG Enrichment (Mtb code: mtu)
kegg_enrich <- enrichKEGG(gene = up_genes,
                          organism = 'mtu',
                          keyType = 'kegg',
                          pvalueCutoff = 0.05)

kegg_res_enr_df <- as.data.frame(kegg_enrich)

kegg_deplete <- enrichKEGG(gene = down_genes,
                           organism = 'mtu',
                           keyType = 'kegg', 
                           pvalueCutoff = 0.05)

kegg_res_dep_df <- as.data.frame(kegg_deplete)

# 3. Dot plot visualization
dotplot(kegg_enrich, showCategory = 10) +
  ggtitle("Enriched KEGG Pathways (Up-regulated in MTB+DCA)")

dotplot(kegg_deplete, showCategory = 10) +
  ggtitle("Enriched KEGG Pathways (Up-regulated in MTB)")

# Bi-derectional plot
kegg_up_df <- as.data.frame(kegg_enrich) %>%
  top_n(10, wt = -p.adjust) %>%
  mutate(Score = -log10(p.adjust),
         Group = "Mtb+DCA")

kegg_down_df <- as.data.frame(kegg_deplete) %>%
  top_n(10, wt = -p.adjust) %>%
  mutate(Score = -log10(p.adjust) * -1, 
         Group = "Mtb")

kegg_tornado_data <- rbind(kegg_up_df, kegg_down_df)

ggplot(kegg_tornado_data, aes(x = reorder(Description, Score), y = Score, fill = Group)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Mtb" = "#E63943", "Mtb+DCA" = "#2A9D8F")) +
  scale_y_continuous(labels = abs) +
  theme_bw() +
  labs(title = "KEGG Pathway Enrichment: Mtb vs Mtb+DCA",
       subtitle = "Left: Higher in Mtb | Right: Higher in Mtb+DCA",
       x = "KEGG Pathways",
       y = "-log10(Adjusted P-value)") +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

# GO pathway analysis
# 1. Mycobrowser tap-separated files loading
mycobrowser_data <- read.delim("Mycobacterium_tuberculosis_H37Rv_txt_v5.txt", sep="\t", stringsAsFactors = FALSE)

# 2. GO mapping table creation
# Locus(Rv number) Gene.Ontology column
go_anno <- mycobrowser_data %>%
  select(Locus, Gene.Ontology) %>%
  filter(Gene.Ontology != "" & !is.na(Gene.Ontology)) %>%
  separate_rows(Gene.Ontology, sep = "[,; ]+") %>%
  mutate(Gene.Ontology = trimws(Gene.Ontology)) %>%
  filter(grepl("GO:", Gene.Ontology)) %>% 
  select(Gene.Ontology, Locus) 

# GO analysis
go_enrich <- enricher(gene = up_genes,
                      TERM2GENE = go_anno,
                      pvalueCutoff = 0.05)

go_deplet <- enricher(gene = down_genes,
                      TERM2GENE = go_anno,
                      pvalueCutoff = 0.05)

go_enrich_df <- as.data.frame(go_enrich)
go_deplet_df <- as.data.frame(go_deplet)

write.csv(go_enrich_df, file="go_enrich.csv")
write.csv(go_deplet_df, file="go_deplet.csv")

dotplot(go_enrich, showCategory = 20) + ggtitle("GO Enrichment: Up-regulated in MTB+DCA")
dotplot(go_deplet, showCategory = 20) + ggtitle("GO Enrichment: Up-regulated in MTB")

#bi-directional bar chart
# MTB+DCA (Up)
up_go <- go_enrich_df %>%
  top_n(10, wt = -p.adjust) %>% 
  mutate(Score = -log10(p.adjust),
         Group = "MTB+DCA")

# MTB (Down)
down_go <- go_deplet_df %>%
  top_n(10, wt = -p.adjust) %>%
  mutate(Score = -log10(p.adjust) * -1, 
         Group = "MTB")

tornado_data <- rbind(up_go, down_go)

ggplot(tornado_data, aes(x = reorder(Description, Score), y = Score, fill = Group)) +
  geom_col(width = 0.8) +
  coord_flip() + 
  scale_fill_manual(values = c("MTB" = "#E63943", "MTB+DCA" = "#2A9D8F")) +
  scale_y_continuous(labels = abs) +
  theme_bw() +
  labs(title = "Top Enriched GO Terms (Biological Process)",
       subtitle = "Left: Enriched in MTB | Right: Enriched in MTB+DCA",
       x = "GO Terms",
       y = "-log10(Adjusted P-value)") +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank())

# specific GO visualization
go_names <- c(
  "GO:0031177" = "Phosphopantetheine binding",
  "GO:0009058" = "Biosynthetic process",
  "GO:0006633" = "Fatty acid biosynthetic process",
  "GO:0006950" = "Response to stress"
)

up_go_clean <- go_enrich_df %>%
  mutate(Description = ifelse(ID %in% names(go_names), go_names[ID], ID),
         Score = -log10(p.adjust),
         Group = "Mtb+DCA")

down_go_clean <- go_deplet_df %>%
  mutate(Description = ifelse(ID %in% names(go_names), go_names[ID], ID),
         Score = -log10(p.adjust) * -1,
         Group = "Mtb")

tornado_data <- bind_rows(up_go_clean, down_go_clean) %>%
  filter(!is.na(Description)) # 매핑된 것만 남기기

ggplot(tornado_data, aes(x = reorder(Description, Score), y = Score, fill = Group)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Mtb" = "#E63943", "Mtb+DCA" = "#2A9D8F")) +
  scale_y_continuous(labels = abs) +
  theme_bw() +
  labs(title = "Top Enriched GO Terms",
       subtitle = "Left: Enriched in MTB | Right: Enriched in MTB+DCA",
       x = "Biological Pathway",
       y = "-log10(Adjusted P-value)",
       fill = "Enriched Group") +
  theme(axis.text.y = element_text(size = 11, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

# Heatmap
sig_genes_df <- res_final %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  filter(!is.na(gene) & gene != locus_tag) %>%     # Rv 번호인 경우 제외
  arrange(desc(log2FoldChange)) # 발현 차이 순서로 정렬

new_row_names <- paste0(sig_genes_df$gene, " (", sig_genes_df$locus_tag, ")")
sig_locus_tags <- sig_genes_df$locus_tag
plot_matrix <- assay(vsd)[sig_locus_tags, ]
rownames(plot_matrix) <- new_row_names

annotation_col <- col_data %>% select(condition)
ann_colors <- list(
  condition = c(MTB = "#E63943", MTB_DCA = "#2A9D8F")
)

pheatmap(plot_matrix, 
         scale = "row",                      
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         show_colnames = TRUE,
         show_rownames = TRUE, 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         main = "Heatmap of DEGs: MTB vs MTB+DCA",
         border_color = NA)

# manual heatmap - each categories visualization
manual_list <- list(
  "Virulence & Secretion" = c("Rv3875", "Rv3874", "Rv3615c", "Rv1038c", "Rv1197", "Rv0288", "Rv0287", "Rv1168c", "Rv0304c", "Rv3892c", "Rv0240"),
  "Dormancy & Stress Response" = c("Rv3133c", "Rv3132c", "Rv2031c", "Rv3130c", "Rv1737c", "Rv1736c", "Rv2626c", "Rv1720"),
  "Lipid Metabolism & Uptake" = c("Rv0165c", "Rv0166", "Rv2482c", "Rv1639c", "Rv0503c", "Rv3229c", "Rv1569", "Rv0468", "Rv0243"),
  "Cell Wall & Drug Resistance" = c("Rv0342", "Rv0341", "Rv0678", "Rv0353", "Rv2461", "Rv3625c", "Rv3791")
)

col_data_final <- data.frame(
  condition = factor(c("Mtb", "Mtb", "Mtb+DCA", "Mtb+DCA"), levels = c("Mtb", "Mtb+DCA")),
  row.names = c("Mtb_1", "Mtb_2", "Mtb+DCA_1", "Mtb+DCA_2")
)

ann_colors <- list(
  condition = c("Mtb" = "#E63943", "Mtb+DCA" = "#2A9D8F")
)

heatmap_results <- list()

for (cat_name in names(manual_list)) {
  
  tags <- manual_list[[cat_name]]
  target_info <- res_final %>% filter(locus_tag %in% tags)
  valid_tags <- target_info$locus_tag
  
  if(length(valid_tags) < 1) next
  
  plot_mat <- assay(vsd)[valid_tags, ]
  rownames(plot_mat) <- paste0(target_info$gene, " (", target_info$locus_tag, ")")
  colnames(plot_mat) <- c("Mtb_1", "Mtb_2", "Mtb+DCA_1", "Mtb+DCA_2")
  
  p <- pheatmap(plot_mat, 
                scale = "row",
                cluster_cols = TRUE,
                annotation_col = col_data_final,
                annotation_colors = ann_colors,
                #color = colorRampPalette(c("#E63943", "white", "#2A9D8F"))(100),
                main = paste("Category:", cat_name),
                border_color = "grey90",
                fontsize_row = 10,
                silent = TRUE)
  
  heatmap_results[[cat_name]] <- p
}

# 1. Virulence Heatmap 
grid.newpage()
grid.draw(heatmap_results[["Virulence & Secretion"]]$gtable)

# 2. Dormancy Heatmap 
grid.newpage()
grid.draw(heatmap_results[["Dormancy & Stress Response"]]$gtable)

# 3. Lipid Metabolism Heatmap 
grid.newpage()
grid.draw(heatmap_results[["Lipid Metabolism & Uptake"]]$gtable)

# 4. Cell Wall & Drug Resistance Heatmap
grid.newpage()
grid.draw(heatmap_results[["Cell Wall & Drug Resistance"]]$gtable)
