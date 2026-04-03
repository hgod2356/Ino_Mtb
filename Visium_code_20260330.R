#start
setwd("G:/CNU/visium")

# ==============================================================================
# Project: Spatial Transcriptomic Analysis of MDR-Mtb Infection and Cs Treatment
# Script: Processing, Deconvolution Integration, and Niche Analysis
# Author: June-Young Lee (KAIST)
# ==============================================================================

# [Setup] Load Required Libraries
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)
library(reshape2)
library(Matrix)
library(ggrepel)

# ==============================================================================
# STEP 1. Data Loading and Quality Control
# ==============================================================================

# 1.1. Load 10X Genomics Spatial Datasets
con <- Load10X_Spatial(data.dir = "MDR_Con/", filename = "Spatial_matrix/filtered_feature_bc_matrix.h5", slice = "MDR_Con")
cs  <- Load10X_Spatial(data.dir = "MDR_Cs/",  filename = "Spatial_matrix/filtered_feature_bc_matrix.h5", slice = "MDR_Cs")

con$orig.ident <- "MDR_Con"
cs$orig.ident  <- "MDR_Cs"

con$barcodes <- rownames(con@meta.data)
cs$barcodes <- rownames(cs@meta.data)

# 1.2. Quality Filtering (Remove low-quality spots)
summary(con$nCount_Spatial)
summary(cs$nCount_Spatial)

summary(con$nFeature_Spatial)
summary(cs$nFeature_Spatial)

MDR_combined <- merge(con, cs, add.cell.ids = c("CON", "CS"))

VlnPlot(MDR_combined, features = c("nCount_Spatial", "nFeature_Spatial"), cols = c("#0072b2", "#EE6677"),
        group.by = "orig.ident", pt.size = 0.1)

SpatialFeaturePlot(MDR_combined, features = "nCount_Spatial") + theme(legend.position = "right")
SpatialFeaturePlot(MDR_combined, features = "nFeature_Spatial") + theme(legend.position = "right")

con$lowQC <- ifelse(con$nFeature_Spatial >300, "PASS", "FAIL")
cs$lowQC <- ifelse(cs$nFeature_Spatial >300, "PASS", "FAIL")

plot1 <- FeatureScatter(con, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", group.by = "lowQC", cols = c("PASS" = "black", "FAIL" = "red")) + theme(legend.position = "none")
print(plot1)
plot2 <- FeatureScatter(cs, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", group.by = "lowQC", cols = c("PASS" = "black", "FAIL" = "red"))
print(plot2)

con <- subset(con, subset = lowQC == "PASS") 
cs <- subset(cs, subset = lowQC == "PASS") 

MDR_combined <- merge(con, cs, add.cell.ids = c("MDR_Con", "MDR_CS"))

# 1.3. Merge Objects
MDR_combined <- JoinLayers(MDR_combined) # Required for Seurat v5

# ==============================================================================
# STEP 2. Normalization and Batch Correction (Harmony)
# ==============================================================================

MDR_combined <- NormalizeData(MDR_combined, normalization.method = "LogNormalize")
MDR_combined <- FindVariableFeatures(MDR_combined, assay = "Spatial", selection.method = "vst", nfeatures = 2000)
MDR_combined <- ScaleData(MDR_combined, features = rownames(MDR_combined))

MDR_combined <- RunPCA(MDR_combined, features = VariableFeatures(object = MDR_combined))
ElbowPlot(MDR_combined, ndims =50)

MDR_combined <- MDR_combined %>% RunUMAP(dims = 1:10) %>%
  FindNeighbors(dims = 1:10) %>% FindClusters(verbose = FALSE, resolution = 1.0)

DimPlot(MDR_combined, split.by = "orig.ident", label = T)

# Run Harmony to mitigate batch effects between samples
MDR_combined <- RunHarmony(MDR_combined, group.by.vars = "orig.ident")

# UMAP and Clustering based on Harmony embeddings
MDR_combined <- MDR_combined %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8)

DimPlot(MDR_combined, split.by = "orig.ident", label = T)

SpatialDimPlot(MDR_combined, images = "MDR_Con") + 
  SpatialDimPlot(MDR_combined, images = "MDR_Cs")

# Nodal clusters group (cluster tree) ------------------------------------------------------------
MDR_combined <- BuildClusterTree(
  MDR_combined,
  dims = 1:20,
  reduction = "harmony",
  reorder = FALSE,
  verbose = TRUE
)

a <- Tool(object = MDR_combined, slot = "BuildClusterTree")
ape::plot.phylo(x = a, direction = "rightwards", label.offset=TRUE)

# Tree visualization
PlotClusterTree(MDR_combined)

# Cluster check
Idents(MDR_combined) <- MDR_combined$seurat_clusters

# Cluster -> grouping
new_ids <- c("0" = "Group1", "1" = "Group1", 
             "2" = "Group2", 
             "3" = "Group3", 
             "4" = "Group4", 
             "5" = "Group5", "7" = "Group5",
             "6" = "Group6", 
             "8" = "Group7",
             "9" = "Group8", "10" = "Group8"
)

# grouping apply
MDR_combined <- RenameIdents(MDR_combined, new_ids)
MDR_combined$groups <- Idents(MDR_combined)

# Grouping check 
SpatialDimPlot(MDR_combined, images = "MDR_Con") + 
  SpatialDimPlot(MDR_combined, images = "MDR_Cs")

# group array ->factoring
group_levels <- c("Group1", "Group2", "Group3", "Group4", "Group5", "Group6", "Group7", "Group8")
MDR_combined$groups <- factor(MDR_combined$groups, levels = group_levels)
Idents(MDR_combined) <- MDR_combined$groups

# group coloring
group_cols = c(
  "Group1" = "#7bc07d", # Green
  "Group2" = "#3467a6", # Blue
  "Group3" = "#7BD3EA", # Cyan
  "Group4" = "#f4bb84", # Orange
  "Group5" = "#de1779", # Pink
  "Group6" = "#A084E8", # Purple (추가)
  "Group7" = "#F9D923", # Yellow (추가)
  "Group8" = "#607274"  # Grey-Blue (추가)
)

DimPlot(MDR_combined, label = T, split.by = "orig.ident", cols = group_cols, label.size = 3, repel = T)

# Control visualization
p1 <- SpatialDimPlot(MDR_combined, pt.size.factor = 2, stroke = NA, label = FALSE, images = "MDR_Con") + 
  scale_fill_manual(values = group_cols) + 
  theme(legend.position = "right") + 
  ggtitle("Control (MDR_Con)")

# MDR+Exp 샘플 시각화
p2 <- SpatialDimPlot(MDR_combined, pt.size.factor = 2, stroke = NA, label = FALSE, images = "MDR_Cs") + 
  scale_fill_manual(values = group_cols) + 
  theme(legend.position = "right") + 
  ggtitle("MDR+Exp (MDR_Cs)")

p1 + p2

meta_data <- MDR_combined@meta.data

# compare number of spots
plot_count <- ggplot(meta_data, aes(x = orig.ident, fill = groups)) +
  geom_bar(color = "black", width = 0.5) +
  theme_bw() +
  labs(title = "Cluster Count", x = "Condition", y = "Number of Spots") +
  scale_fill_manual(values = group_cols)

# compare spot percentage
plot_percent <- ggplot(meta_data, aes(x = orig.ident, fill = groups)) +
  geom_bar(position = "fill", color = "black", width = 0.5) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Cluster Proportion", x = "Condition", y = "Percentage") +
  scale_fill_manual(values = group_cols)

plot_count + plot_percent

#all group marker gene finding
group_markers <- FindAllMarkers(MDR_combined, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)

#write.csv(group_markers, file = "MDR_Group_Markers_Result_20260224.csv")

# ==============================================================================
# STEP 3. Integration of Cell2location Results
# ==============================================================================

# Count file extraction for Cell2location
counts_vis <- GetAssayData(MDR_combined, assay = "Spatial", layer = "counts")
writeMM(counts_vis, "vis_counts.mtx")

# Save metadata
#write.csv(MDR_combined@meta.data, "vis_metadata.csv", row.names = TRUE)

# Save genename
#write.csv(rownames(MDR_combined), "vis_genes.csv", row.names = FALSE)

#-> after cell2location

# 3.1. Import Cell2location output
c2l_res <- read.csv("cell2location_MDR_results.csv", row.names = 1)
MDR_combined <- AddMetaData(MDR_combined, metadata = c2l_res)

MDR_combined@meta.data


# full name mapping
name_map <- c(
  "AT1" = "Alveolar Type 1 cell",
  "AT2.1" = "Alveolar Type 2 cell (1)",
  "AT2.2" = "Alveolar Type 2 cell (2)",
  "Alv.Mf" = "Alveolar Macrophage",
  "Art" = "Arterial Endothelial cell",
  "B.cell.1" = "B cell (1)",
  "B.cell.2" = "B cell (2)",
  "CD4.T.cell.1" = "CD4+ T cell (1)",
  "CD4.T.cell.2" = "CD4+ T cell (2)",
  "CD8.T.cell.1" = "CD8+ T cell (1)",
  "CD8.T.cell.2" = "CD8+ T cell (2)",
  "Cap" = "General Capillary Endothelial",
  "Cap.a" = "Aerocyte (Alveolar Capillary)",
  "Ciliated" = "Ciliated Epithelial cell",
  "Club" = "Club cell",
  "Col13a1..fibroblast" = "Col13a1+ Adventitial Fibroblast",
  "Col14a1..fibroblast" = "Col14a1+ Alveolar Fibroblast",
  "DC1" = "Dendritic Cell (Type 1)",
  "DC2" = "Dendritic Cell (Type 2)",
  "ILC2" = "Innate Lymphoid Cell Type 2",
  "Int.Mf" = "Interstitial Macrophage",
  "Lymph" = "Lymphatic Endothelial cell",
  "Mast.Ba2" = "Mast cell / Basophil",
  "Mesothelial" = "Mesothelial cell",
  "Mono" = "Monocyte",
  "Myofibroblast" = "Myofibroblast",
  "NK.cell" = "Natural Killer cell",
  "Neut.1" = "Neutrophil (Type 1)",
  "Neut.2" = "Neutrophil (Type 2)",
  "Pericyte.1" = "Pericyte (1)",
  "Pericyte.2" = "Pericyte (2)",
  "SMC" = "Smooth Muscle Cell",
  "Vein" = "Venous Endothelial cell",
  "gd.T.cell" = "Gamma-Delta T cell"
)

# spot - cell type matching
max_cell_type <- colnames(c2l_res)[apply(c2l_res, 1, which.max)]
MDR_combined$predicted_celltype <- max_cell_type
table(MDR_combined$predicted_celltype)

DimPlot(MDR_combined, 
        group.by = "predicted_celltype", 
        split.by = "orig.ident", 
        label = T, 
        label.size = 3, 
        repel = T) + 
  ggtitle("Cell Type Distribution")+
    theme(legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

# label->full name change
MDR_combined$predicted_celltype_full <- ifelse(
  MDR_combined$predicted_celltype %in% names(name_map),
  name_map[MDR_combined$predicted_celltype],
  MDR_combined$predicted_celltype
)

table(MDR_combined$predicted_celltype_full)
table(MDR_combined$predicted_celltype)

DimPlot(MDR_combined, 
        group.by = "predicted_celltype_full", 
        split.by = "orig.ident", 
        label = T, 
        label.size = 3, 
        repel = T) + 
  ggtitle("Cell Type Distribution (Predicted by Cell2location)") + 
  theme(legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

# Merge Cell Lineages (16 Major Groups)
merge_list <- list(
  "Neutrophils"               = c("Neut.1", "Neut.2"),
  "Alveolar macrophages"      = "Alv.Mf",
  "Interstitial macrophages"  = "Int.Mf",
  "Monocytes"                 = "Mono",
  "Dendritic cells"           = c("DC1", "DC2"),
  "Alveolar epithelial cells" = c("AT1","AT2.1", "AT2.2"),
  "Endothelial cells"         = c("Cap", "Cap.a", "Art", "Vein", "Lymph"),
  "Fibroblasts"               = c("Col13a1..fibroblast", "Col14a1..fibroblast", "Myofibroblast"),
  "CD4+ T cells"              = c("CD4.T.cell.1", "CD4.T.cell.2"),
  "CD8+ T cells"              = c("CD8.T.cell.1", "CD8.T.cell.2"),
  "B cells"                   = c("B.cell.1", "B.cell.2"),
  "Natural Killer cells"      = "NK.cell",
  "Mast cells / Basophils"    = "Mast.Ba2",
  "Other T cells"             = c("gd.T.cell", "ILC2"),
  "Airway epithelial cells"   = c("Ciliated", "Club"),
  "Others"                    = c("SMC", "Pericyte.1", "Pericyte.2", "Mesothelial")
)

# Calculate aggregated abundance for merged groups
abundance_final <- data.frame(matrix(NA, nrow = ncol(MDR_combined), ncol = length(merge_list)))
colnames(abundance_final) <- names(merge_list)
for(group in names(merge_list)){
  abundance_final[, group] <- rowSums(c2l_res[, merge_list[[group]], drop = FALSE])
}
MDR_combined <- AddMetaData(MDR_combined, metadata = abundance_final)

# Assign Dominant Cell Type to each spot
MDR_combined$main_cellgroup <- factor(names(merge_list)[apply(abundance_final, 1, which.max)], levels = names(merge_list))

table(MDR_combined$main_cellgroup)

DimPlot(MDR_combined, 
        group.by = "main_cellgroup", 
        split.by = "orig.ident", 
        label = T, 
        label.size = 3, 
        repel = T) + 
  ggtitle("Cell Type Distribution (Predicted by Cell2location)") + 
  theme(legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

#Neutrophil, Alv.Mf abundance spatial plot
MDR_combined$Neutrophils_abundance <- abundance_final$Neutrophils
MDR_combined$Alvmf_abundance <- abundance_final$`Alveolar macrophages`

FeaturePlot(MDR_combined, 
            features = "Neutrophils_abundance",
            split.by = "orig.ident", 
            cols = c("lightgrey", "red"))

p123 <- SpatialFeaturePlot(MDR_combined, 
                           features = "Neutrophils_abundance", 
                           images = c("MDR_Con", "MDR_Cs"), 
                           pt.size.factor = 1.6,
                           min.cutoff = "q0",
                           max.cutoff = "q99"
)

p234 <- SpatialFeaturePlot(MDR_combined, 
                           features = "Alvmf_abundance", 
                           images = c("MDR_Con", "MDR_Cs"),
                           pt.size.factor = 1.6,
                           min.cutoff = "q0",
                           max.cutoff = "q99"
)

# Color and title
p123 & scale_fill_gradientn(colors = c("lightgrey", "yellow", "red")) &
  patchwork::plot_annotation(title = "Neutrophils Abundance")

p234 & scale_fill_gradientn(colors = c("lightgrey", "yellow", "red")) &
  patchwork::plot_annotation(title = "Alveolar Macrophages Abundance")

# ==============================================================================
# STEP 4. Spatial Niche Stratification
# ==============================================================================

# 1. Neutrophil marker list (just select)
neut_features <- list(c("S100a8", "S100a9", "Ly6g"))

# 2. Neutrophil Score (Module Score)
MDR_combined <- AddModuleScore(
  object = MDR_combined,
  features = neut_features,
  name = "Neutrophil_Score"
)

Neut.scoreplot<-SpatialFeaturePlot(MDR_combined, features = "Neutrophil_Score1", 
                   images = c("MDR_Con", "MDR_Cs"), pt.size.factor = 1.6)

Neut.scoreplot + 
  plot_annotation(title = "Neutrophils score plot",
                  theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))

# Neutrophil Score compare Vlnplot (Control vs Cs)
VlnPlot(MDR_combined, features = "Neutrophil_Score1", group.by = "orig.ident", 
        cols = c("#0072b2", "#EE6677"), pt.size = 0) + 
  geom_boxplot(width=0.1, fill="white") +
  labs(title = "Neutrophil Score Comparison", y = "Score (S100a8, S100a9, Ly6g)")

# 2. 통계적 유의성 확인 (Wilcoxon rank sum test)
score_data <- MDR_combined@meta.data
wilcox_res <- wilcox.test(Neutrophil_Score1 ~ orig.ident, data = score_data)
print(wilcox_res$p.value) # p-value가 0.05보다 작으면 유의미한 변화입니다.

# 그룹별 Neutrophil Score 확인
VlnPlot(MDR_combined, features = "Neutrophil_Score1", group.by = "groups", 
        cols = group_cols) + 
  ggtitle("Neutrophil Infiltration by Spatial Groups")

# ==============================================================================
# STEP 5. Marker gene identification
# ==============================================================================

# Reference data load - GSE151974
ref_data <- read.csv("GSE151974_reference_signatures.csv", row.names = 1)

# specific marker gene expression in cell
check_gene <- function(gene_name) {
  if(gene_name %in% rownames(ref_data)) {
    # 행 데이터를 추출하고 벡터로 변환 (unlist)
    gene_values <- unlist(ref_data[gene_name, ])
    # 내림차순 정렬
    return(sort(gene_values, decreasing = TRUE))
  } else {
    return("Gene not found")
  }
}

check_gene("Cd36")
check_gene("S100a9")
check_gene("Fth1")
check_gene("Tmsb4x")

# Reference fixing <- merge list apply
merged_ref <- matrix(0, nrow = nrow(ref_data), ncol = length(merge_list))
colnames(merged_ref) <- names(merge_list)
rownames(merged_ref) <- rownames(ref_data)

for (category in names(merge_list)) {
  subtypes <- merge_list[[category]]
  subtypes_fixed <- gsub("[+ /]", ".", subtypes)
  
  valid_cols <- subtypes_fixed[subtypes_fixed %in% colnames(ref_data)]
  
  if (length(valid_cols) > 1) {
    merged_ref[, category] <- rowMeans(ref_data[, valid_cols]) # 평균값 사용
  } else if (length(valid_cols) == 1) {
    merged_ref[, category] <- ref_data[, valid_cols]
  }
}
merged_ref <- as.data.frame(merged_ref)

get_top_genes_no_filter <- function(df, target_cat, top_n = 50) {
  target_vals <- df[[target_cat]]
  names(target_vals) <- rownames(df)
  return(sort(target_vals, decreasing = TRUE)[1:top_n])
}

top_neut_list <- get_top_genes_no_filter(merged_ref, "Neutrophils")
top_am_list <- get_top_genes_no_filter(merged_ref, "Alveolar macrophages")

print("--- Neutrophils Top Genes (No Filter) ---")
print(names(top_neut_list))

print("--- Alveolar Macrophages Top Genes (No Filter) ---")
print(names(top_am_list))

intersect_genes <- intersect(names(top_neut_list), names(top_am_list))
print(intersect_genes)

unique_neut_genes <- setdiff(names(top_neut_list), names(top_am_list))
unique_am_genes <- setdiff(names(top_am_list), names(top_neut_list))

print(unique_neut_genes)
print(unique_am_genes)


# ==============================================================================
# STEP 6. module scoring - neutrophils, macrophages
# ==============================================================================

#ADD module score for specific cell type
am_list <- list(c("Fth1", "Tmsb4x", "Ftl1", "Ccl6"))
neut_list <- list(c("Ly6g", "S100a8", "S100a9"))

# AddModuleScore
MDR_combined <- AddModuleScore(
  object = MDR_combined,
  features = am_list,
  name = "AM_Specific_Score"
)

MDR_combined <- AddModuleScore(
  object = MDR_combined,
  features = neut_list,
  name = "Neut_Specific_Score"
)

Am_scorespatial <- SpatialFeaturePlot(MDR_combined, features = "AM_Specific_Score1")
print(Am_scorespatial)

Neut_scorespatial <- SpatialFeaturePlot(MDR_combined, features = "Neut_Specific_Score1")
print(Neut_scorespatial)

# score comparison
my_comparisons <- list( c("MDR_Con", "MDR_Cs") )

# AM specific score Vln plot
max_score_Am <- max(MDR_combined$AM_Specific_Score1, na.rm = TRUE)
plot_top_Am <- ceiling(max_score_Am * 1.1 / 0.5) * 0.5 # 최대값보다 높으면서 0.5의 배수가 되도록 설정

VlnPlot(MDR_combined, 
        features = "AM_Specific_Score1", 
        group.by = "orig.ident", 
        cols = c("#0072b2", "#EE6677"), 
        pt.size = 0,
        y.max = plot_top_Am) + 
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) + 
  labs(title = "Alveolar Macrophage Score Comparison", y = "Activity Score") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = max_score_Am * 1.05) +
  scale_y_continuous(breaks = seq(0, plot_top_Am, by = 0.5),
                     limits = c(NA, plot_top_Am)) 

max_score_neu <- max(MDR_combined$Neut_Specific_Score1, na.rm = TRUE)
plot_top_net <- ceiling(max_score_neu * 1.1 / 0.5) * 0.5 

VlnPlot(MDR_combined, 
        features = "Neut_Specific_Score1", 
        group.by = "orig.ident", 
        cols = c("#0072b2", "#EE6677"), 
        pt.size = 0,
        y.max = plot_top_net) + 
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) + 
  labs(title = "Neutrophils Score Comparison", y = "Activity Score") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = max_score_neu * 1.1) +
  scale_y_continuous(breaks = seq(0, plot_top_net, by = 0.5),
                     limits = c(NA, plot_top_net))           

# ==============================================================================
# STEP 7. Gene spatial plotting
# ==============================================================================

genes_to_plot <- c("Foxo1", "Foxo3", "Acadl", "Hadha")
target_clusters <- c("Neutrophils", "Alveolar macrophages")

p_genes <- SpatialFeaturePlot(MDR_combined, 
                              features = genes_to_plot, 
                              images = c("MDR_Con", 
                                         "MDR_Cs"), 
                              pt.size.factor = 2,
                              ncol = 4,
                              min.cutoff = "q5", 
                              max.cutoff = "q95")

p_genes + 
  plot_annotation(title = "Expression of Foxo and FAO-related Genes",
                  subtitle = "MDR_Con and MDR_Cs",
                  theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))

# dotplot
Idents(MDR_combined) <- "main_cellgroup"

DotPlot(MDR_combined, 
        features = genes_to_plot, 
        group.by = "main_cellgroup",
        cols = c("dodgerblue3", "firebrick1"), 
        dot.scale = 8) + 
  RotatedAxis() + 
  ggtitle("Gene Expression across Cell Populations")


# ==============================================================================
# STEP 8. Pathway analysis
# ==============================================================================

library(clusterProfiler)
library(org.Mm.eg.db)

# Neutrophils
neut_median <- median(MDR_combined$Neutrophils)
neut_mean <- mean(MDR_combined$Neutrophils)
print(paste("Neutrophil Median Value:", neut_median))
print(paste("Neutrophil Mean Value:", neut_mean))

MDR_combined$neut_group_median <- ifelse(MDR_combined$Neutrophils > neut_median, 
                                         "Neut_High", "Neut_Low")
MDR_combined$neut_group_mean <- ifelse(MDR_combined$Neutrophils > neut_mean, 
                                       "Neut_High", "Neut_Low")

table(MDR_combined$orig.ident, MDR_combined$neut_group_median)
table(MDR_combined$orig.ident, MDR_combined$neut_group_mean)

# Tagging
MDR_combined$target_compare <- "Others"
MDR_combined$target_compare[MDR_combined$orig.ident == "MDR_Con" & MDR_combined$neut_group_median == "Neut_High"] <- "High_Con"
MDR_combined$target_compare[MDR_combined$orig.ident == "MDR_Cs" & MDR_combined$neut_group_median == "Neut_Low"] <- "Low_Cs"

MDR_combined$target_compare_mean <- "Others" #기존데이터 덮어쓰기로 중복 제거
MDR_combined$target_compare_mean[MDR_combined$orig.ident == "MDR_Con" & MDR_combined$neut_group_mean == "Neut_High"] <- "High_Con"
MDR_combined$target_compare_mean[MDR_combined$orig.ident == "MDR_Cs"  & MDR_combined$neut_group_mean == "Neut_Low"] <- "Low_Cs"

table(MDR_combined$target_compare)
table(MDR_combined$target_compare_mean)

# DEG (High_Con vs Low_Cs)
Idents(MDR_combined) <- "target_compare"
final_degs <- FindMarkers(MDR_combined, 
                          ident.1 = "High_Con", 
                          ident.2 = "Low_Cs",
                          logfc.threshold = 0.25)

# MDR_Con (Neut High)
ego_high <- enrichGO(gene = rownames(final_degs[final_degs$avg_log2FC > 0.25 & final_degs$p_val_adj < 0.05, ]),
                     OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', ont = "BP")

# MDR_Cs (Neut Low)
ego_low <- enrichGO(gene = rownames(final_degs[final_degs$avg_log2FC < -0.25 & final_degs$p_val_adj < 0.05, ]),
                    OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', ont = "BP")

df_high <- as.data.frame(ego_high)
if(nrow(df_high) > 0) df_high$Group <- "MDR_Con (Neut-High)"

df_low <- as.data.frame(ego_low)
if(nrow(df_low) > 0) df_low$Group <- "MDR_Cs (Neut-Low)"

#write.csv(as.data.frame(ego_high), file = "GO_Pathways_MDR_Con_High.csv", row.names = FALSE)
#write.csv(as.data.frame(ego_low), file = "GO_Pathways_MDR_Cs_Low.csv", row.names = FALSE)

plot_df <- rbind(head(df_high, 10), head(df_low, 10))

# visualization
ggplot(plot_df, aes(x = Group, y = Description, size = GeneRatio, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(title = "Comparison of Biological Processes", x = "", y = "GO Pathway") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Alv.Mf
am_median <- median(MDR_combined$Alv.Mf)
am_mean <- mean(MDR_combined$Alv.Mf)
print(paste("Alveolar Macrophage Median Value:", am_median))
print(paste("Alveolar Macrophage Mean Value:", am_mean))

# tagging
MDR_combined$am_group_median <- ifelse(MDR_combined$Alv.Mf > am_median, "AM_High", "AM_Low")
MDR_combined$am_group_mean <- ifelse(MDR_combined$Alv.Mf > am_mean, "AM_High", "AM_Low")

table(MDR_combined$orig.ident, MDR_combined$am_group_median)
table(MDR_combined$orig.ident, MDR_combined$am_group_mean)

# AmHigh_Con vs AmLow_Cs)
MDR_combined$am_target_compare <- "Others"
MDR_combined$am_target_compare[MDR_combined$orig.ident == "MDR_Con" & MDR_combined$am_group_median == "AM_High"] <- "AM_High_Con"
MDR_combined$am_target_compare[MDR_combined$orig.ident == "MDR_Cs" & MDR_combined$am_group_median == "AM_Low"] <- "AM_Low_Cs"

MDR_combined$am_target_compare <- "Others"
MDR_combined$am_target_compare[MDR_combined$orig.ident == "MDR_Con" & MDR_combined$am_group_mean == "AM_High"] <- "AM_High_Con"
MDR_combined$am_target_compare[MDR_combined$orig.ident == "MDR_Cs"  & MDR_combined$am_group_mean == "AM_Low"] <- "AM_Low_Cs"

table(MDR_combined$am_target_compare)

# DEG 
Idents(MDR_combined) <- "am_target_compare"
am_final_degs <- FindMarkers(MDR_combined, 
                             ident.1 = "AM_High_Con", 
                             ident.2 = "AM_Low_Cs",
                             logfc.threshold = 0.25)

#write.csv(am_final_degs, file = "DEG_AlvMf_HighCon_vs_LowCs.csv")

# AmHigh_Con
ego_am_high <- enrichGO(gene = rownames(am_final_degs[am_final_degs$avg_log2FC > 0.25 & am_final_degs$p_val_adj < 0.05, ]),
                        OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', ont = "BP")

# AMLow_Cs
ego_am_low <- enrichGO(gene = rownames(am_final_degs[am_final_degs$avg_log2FC < -0.25 & am_final_degs$p_val_adj < 0.05, ]),
                       OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', ont = "BP")

#write.csv(as.data.frame(ego_am_high), "GO_AlvMf_High_Con.csv")
#write.csv(as.data.frame(ego_am_low), "GO_AlvMf_Low_Cs.csv")

#visualization
library(dplyr)

# GO results merging
df_neut_high <- as.data.frame(ego_high) %>% mutate(CellType = "Neutrophil", Group = "MDR_Con (High)")
df_neut_low  <- as.data.frame(ego_low)  %>% mutate(CellType = "Neutrophil", Group = "MDR_Cs (Low)")
df_am_high   <- as.data.frame(ego_am_high) %>% mutate(CellType = "Alv.Mf", Group = "MDR_Con (High)")
df_am_low    <- as.data.frame(ego_am_low)  %>% mutate(CellType = "Alv.Mf", Group = "MDR_Cs (Low)")

combined_df <- rbind(df_neut_high, df_neut_low, df_am_high, df_am_low)

target_pathways <- c(
  "lipid localization", 
  "lipid transport", 
  "lipid storage", 
  "neutrophil migration", 
  "lipid homeostasis", 
  "lipid oxidation", 
  "macrophage activation", 
  "myeloid leukocyte activation",
  "leukocyte chemotaxis",
  "fatty acid metabolic process",
  "epithelial cell development",
  "cilium organization"
)

# 3. 정확히 일치하는 항목만 필터링
filtered_df <- combined_df %>%
  filter(tolower(Description) %in% tolower(target_pathways))

# 4. 시각화 (Bubble Plot)
ggplot(filtered_df, aes(x = Group, y = Description, size = FoldEnrichment, color = p.adjust)) +
  geom_point() +
  facet_wrap(~CellType, scales = "free_x") +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  labs(
    title = "GO Pathways analysis",
    subtitle = "Neutrophils and Alveolar Macrophages Comparison",
    x = "Comparison Group",
    y = "Gene Ontology Terms",
    size = "Fold enrichment",
    color = "Adj. P-value"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10),
    strip.background = element_rect(fill = "gray20"),
    strip.text = element_text(color = "white", face = "bold")
  )

#IFN-beta target pathway
target_pathways2 <- c(
  "cellular response to interferon-beta", 
  "interferon-beta production", 
  "response to interferon-beta", 
  "interferon-mediated signaling pathway", 
  "regulation of interferon-beta production", 
  "innate immune response-activating signaling pathway"
)

# 3. 정확히 일치하는 항목만 필터링
filtered_df2 <- combined_df %>%
  filter(tolower(Description) %in% tolower(target_pathways2))

# 4. 시각화 (Bubble Plot)
ggplot(filtered_df2, aes(x = Group, y = Description, size = FoldEnrichment, color = p.adjust)) +
  geom_point() +
  facet_wrap(~CellType, scales = "free_x") +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  labs(
    title = "GO Pathways analysis",
    subtitle = "Neutrophils and Alveolar Macrophages Comparison",
    x = "Comparison Group",
    y = "Gene Ontology Terms",
    size = "Fold enrichment",
    color = "Adj. P-value"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10),
    strip.background = element_rect(fill = "gray20"),
    strip.text = element_text(color = "white", face = "bold")
  )


# ==============================================================================
# STEP 8. Volcano plot
# ==============================================================================

# DEG (MDR_Con vs MDR_Cs)
cs_markers <- FindMarkers(MDR_combined, ident.1 = "MDR_Con", 
                          group.by ="orig.ident", test.use = "MAST")

write.csv(cs_markers, file = "MDR_Cluster_convscs_markers_20260225.csv")

res_data <- read.csv("MDR_Cluster_convscs_markers_20260225.csv", stringsAsFactors = F)
colnames(res_data)[1] <- "genename"

res_data$p_val_adj[res_data$p_val_adj == 0] <- 1e-300
res_data <- res_data[!is.na(res_data$p_val_adj), ]

# Left: MDR_Con, Right: MDR_Cs
res_data$log2FoldChange <- -(res_data$avg_log2FC)

cut_lfc <- 0.5
cut_pvalue <- 0.05

#LFC negative (MDR_Con), LFC positive (MDR_Cs)
res_data$group <- "Not significant"
res_data$group[res_data$p_val_adj < cut_pvalue & res_data$log2FoldChange < -cut_lfc] <- "MDR_Con"
res_data$group[res_data$p_val_adj < cut_pvalue & res_data$log2FoldChange >  cut_lfc] <- "MDR_Cs"

# legend order
res_data$group <- factor(res_data$group, levels = c("MDR_Con", "MDR_Cs", "Not significant"))

# gene list
target_genes <- c(
  "Srebf1", "Fasn", "Acly", "Acaca", "Acacb", "Gpat1", "Gpat2", "Gpat3", "Gpat4", 
  "Acsl1", "Acsl5", "Lpl", "Slc27a1", "Lpin2", "Lpgat1", "Dgke", "Elovl5", "Abhd8",
  "Sqstm1", "Nbr1", "Optn", "Spart", "Osbpl8", "Ddhd2", "Vps4a", "Atg14", "Tp53inp2",
  "Tfeb", "Tfe3", "Foxo1", "Ppara", "Pparg", "Nr1h3", "Sik1", "Prkaa2", "Prkab2",
  "Gdf15", "Fabp3", "Ppard", "Auh", "Cpt2", "Ldlr", "Igf1r", "Egr1", "Muc2", "Muc19"
)

# gene extraction
label_data <- subset(res_data, genename %in% target_genes)

# visualization
volcano_colors <- c("MDR_Con"="#E63943", "MDR_Cs"="#2A9D8F", "Not significant"="grey80")

ggplot(res_data, aes(x=log2FoldChange, y=-log10(p_val_adj), color=group)) +
  geom_point(alpha=0.6, size=1.8) +
  scale_color_manual(values=volcano_colors) +
  geom_text_repel(data = label_data, 
                  aes(label = genename),
                  size = 3.5,
                  fontface = "plain",
                  color = "black",
                  box.padding = 0.6,
                  point.padding = 0.4,
                  force = 4,
                  segment.color = "grey30",
                  max.overlaps = Inf,
                  show.legend = FALSE) +
  geom_vline(xintercept = c(-cut_lfc, 0, cut_lfc), 
             linetype = c("dashed","dotted","dashed"), color="black") +
  geom_hline(yintercept = -log10(cut_pvalue), linetype="dashed", color="black") +
    labs(title="Volcano plot: MDR_Con (Left) vs MDR_Cs (Right)",
       subtitle="Labeled: Selected Lipid Metabolism & Host Protection Genes",
       x=expression(Log[2]~fold~change),
       y=expression(-Log[10]~Q~value),
       color="Group") +
  theme_bw(base_size=14) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size=10, color="grey30"),
    legend.position = "right",
    legend.box.background = element_rect(color="black", fill="white", size=0.8)
  )


# ==============================================================================
# STEP 9. Pathway scoring spatial plot
# ==============================================================================

# Each gene spatial-------------------------------------------------------------
# pathway related gene spatial feature plot
ifnb_production <- "interferon-beta production"
ifnb_response <- "response to interferon-beta"

ifnb_production_from_go <- ego_high@result %>%
  filter(Description == ifnb_production) %>%
  pull(geneID) %>%
  strsplit("/") %>% # 보통 유전자들이 "Gene1/Gene2/Gene3" 형태로 들어있습니다
  unlist()

ifnb_response_from_go <- ego_high@result %>%
  filter(Description == ifnb_response) %>%
  pull(geneID) %>%
  strsplit("/") %>% # 보통 유전자들이 "Gene1/Gene2/Gene3" 형태로 들어있습니다
  unlist()

SpatialFeaturePlot(MDR_combined, 
                   features = ifnb_production_from_go[1:8],
                   images = c("MDR_Con", "MDR_Cs"), 
                   pt.size.factor = 1.6,
                   ncol = 4, 
                   min.cutoff = "q5", 
                   max.cutoff = "q95") +
  plot_annotation(title = paste("Genes involved in:", ifnb_production),
                  subtitle = "MDR_Con vs MDR_Cs",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

SpatialFeaturePlot(MDR_combined, 
                   features = ifnb_production_from_go[1:8],
                   images = c("MDR_Con", "MDR_Cs"), 
                   pt.size.factor = 1.6,
                   ncol = 4, 
                   min.cutoff = "q5", 
                   max.cutoff = "q95") +
  plot_annotation(title = paste("Genes involved in:", ifnb_response),
                  subtitle = "MDR_Con vs MDR_Cs",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))



# addmodulescore apply----------------------------------------------------------
#GSEA
# 1. TNF-alpha Signaling Gene List (HALLMARK_TNFA_SIGNALING_VIA_NFKB)
tnf_nfkb_genes <- c(
  "Mcl1","Cd80","Traf1","F2rl1","Dusp2","Tnc","Fosl2","Stat5a","Vegfa","Efna1","Relb","Rela",
  "Cebpd","Ptger4","Cdkn1a","Ptx3","Il15ra","Atp2b1","Nfkbia","Tnf","Ier3","Ier2","Hes1",
  "Tnfaip2","Dusp1","Eif1","Fosl1","Bcl6","Ifngr2","Tank","Gadd45b","Gadd45a","Tubb2a",
  "Sqstm1","Il18","Rhob","Cxcl1","Btg2","Nfe2l2","Tsc22d1","Irs2","Atf3","Nfil3","Btg3",
  "Jag1","Ackr3","Cxcl5","Phlda1","Bhlhe40","Per1","Plk2","Nfkb2","Tnfsf9","Tnfrsf9",
  "Klf10","Tgif1","Nfkbie","Tnfaip6","Tnfaip3","Ninj1","Birc3","Birc2","Smad3","Dennd5a",
  "Socs3","Phlda2","Olr1","Snn","Bcl2a1d","Egr3","Sphk1","G0s2","Cd83","Ccl20","Klf9",
  "Msc","Cflar","Ier5","Gfpt2","Csf2","Csf1","Sgk1","Cxcl2","Ehd1","Fjx1","Klf4","Klf2",
  "Tlr2","Klf6","Map2k3","Map3k8","Sdc4","Cxcl10","Nr4a1","Nr4a2","Nr4a3","Icosl","Nfat5",
  "Panx1","Cxcl11","Plek","Rcan1","Ripk2","Dnajb4","Plpp3","Pnrc1","Kynu","Ifih1","Dram1",
  "Ccrl2","Spsb1","Rnf19b","Ccnl1","Tnip1","Ppp1r15a","B4galt5","Pdlim5","Litaf","Pmepa1",
  "Nampt","Clcf1","Il23a","Zbtb10","Slc16a6","Trip10","Tnfaip8","Tiparp","Pfkfb3",
  "Zc3h12a","Tnip2","Yrdc","Gpr183","Dusp4","Rigi","Slc2a6","Trib1","Klf10","Dusp5",
  "Areg","Bcl3","Bmp2","Btg1","Ccnd1","Cd44","Cd69","Cebpb","F3","Ccn1","Serpinb8",
  "Edn1","Egr1","Egr2","Ets2","Fos","Fosb","Fut4","Gch1","B4galt1","Slc2a3","Hbegf",
  "Icam1","Id2","Il12b","Il1a","Il1b","Il6","Il6st","Il7r","Inhba","Irf1","Jun","Junb",
  "Ldlr","Lif","Marcks","Mxd1","Maff","Myc","Nfkb1","Serpine1","Serpinb2","Plau","Plaur",
  "Ptgs2","Ptpre","Rel","Sat1","Ccl5","Sod2","Tap1","Zfp36","Ifit2","Pde4b","Abca1","Gem","Lamb3"
)

# Check gene availability in the object
# Only use genes that actually exist in your spatial data
available_genes <- tnf_nfkb_genes[tnf_nfkb_genes %in% rownames(MDR_combined)]

# Add Module Score
# The score will be saved in the metadata as 'TNF_NFkB_Score1'
MDR_combined <- AddModuleScore(
  object = MDR_combined,
  features = list(available_genes),
  name = "TNF_NFkB_Score"
)

# Spatial Feature Plot Visualization
p_score <- SpatialFeaturePlot(MDR_combined, 
                              features = "TNF_NFkB_Score1", 
                              images = c("MDR_Con", "MDR_Cs"), 
                              pt.size.factor = 2,
                              min.cutoff = "q5", 
                              max.cutoff = "q95") 

# Add Title and Labels
p_score_final <- p_score + 
  plot_annotation(title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

p_score_final

# 2. HALLMARK_IL6_JAK_STAT3_SIGNALING
il6_jak_stat3_genes <- c(
  "Stat3","Stat2","Stat1","Il12rb1","Ccr1","Pla2g2a","Il15ra","Ltb","Tnf","Ltbr",
  "Lepr","Il13ra1","Il4ra","Il18r1","Il17ra","Cd38","Irf9","Ifngr2","Ifngr1","Ifnar1",
  "Cd36","Myd88","Inhbe","Il10rb","Bak1","Socs3","Tnfrsf1b","Tnfrsf1a","Osmr","Acvr1b",
  "Acvrl1","Csf2","Csf1","Csf2ra","Csf3r","Cxcl2","Tlr2","Hax1","Map3k8","Tnfrsf12a",
  "Cxcl9","Cxcl10","Ebi3","Socs1","Il17rb","Pdgfc","Cxcl11","Cxcl13","Pf4","Crlf2",
  "Stam2","Tyk2","Tnfrsf21","Pik3r5","A2m","Cbl","Cd14","Cd44","Cd9","Fas","Grb2",
  "Hmox1","Il1b","Il1r1","Il1r2","Il2ra","Il2rg","Il3ra","Il6","Il6st","Il7","Il9r",
  "Irf1","Itga4","Itgb3","Jun","Pim1","Ptpn1","Ptpn2","Reg1","Dntt","Tgfb1","Ptpn11",
  "Ccl7","Cntfr"
)

# Filter for genes present in the dataset
available_jak_genes <- il6_jak_stat3_genes[il6_jak_stat3_genes %in% rownames(MDR_combined)]

# Calculate Module Score
# The score will be saved as 'JAK_STAT3_Score1'
MDR_combined <- AddModuleScore(
  object = MDR_combined,
  features = list(available_jak_genes),
  name = "JAK_STAT3_Score"
)

# Spatial Feature Plot Visualization
p_jak_score <- SpatialFeaturePlot(MDR_combined, 
                              features = "JAK_STAT3_Score1", 
                              images = c("MDR_Con", "MDR_Cs"), 
                              pt.size.factor = 2,
                              min.cutoff = "q5", 
                              max.cutoff = "q95") 

# Add Title and Labels
p_jak_score_final <- p_jak_score + 
  plot_annotation(title = "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

p_jak_score_final

# 3. HALLMARK_IL2_STAT5_SIGNALING
il2_stat5_genes <- c(
  "Cd86","Traf1","Pou2f1","Ecm1","Serpinb6a","Cish","Fgl2","Tiam1","Tnfrsf4","Cdkn1c",
  "Penk","Rora","Gpx4","Ltb","Ahr","Slc1a5","Il4ra","Ptch1","Adam19","Il18r1","Nrp1",
  "Myo1c","Myo1e","Tnfsf10","Capn3","Scn9a","Ifngr1","Casp3","Gadd45b","Ccr4","Nop2",
  "Lrig1","Emp1","Rhob","Gpr65","Rgs16","Plpp1","Ncs1","Mapkapk2","Nfil3","Bmpr2",
  "Cd81","Irf4","Phlda1","Bhlhe40","Capg","Amacr","Rragd","Hycc2","She","Tnfsf11",
  "Plagl1","Tnfrsf9","Rnh1","Dennd5a","Eomes","Map6","Socs2","Ncoa3","Plec","Coch",
  "F2rl2","Cst7","Itgae","Umps","Swap70","Alcam","Hipk2","Tnfrsf1b","Hk2","Dhrs3",
  "Ahnak","St3gal4","Cd83","Syngr2","Cyfip1","P2rx4","S100a1","Csf2","Csf1","Ndrg1",
  "Gsto1","Ikzf2","Ikzf4","Spry4","Cdc6","Slc29a2","Klf6","Map3k8","Cxcl10","Rabgap1l",
  "Socs1","Icos","Batf","Irf6","Syt11","Praf2","Gucy1b1","Ctsz","Ifitm3","Snx9",
  "Slc39a8","Gabarapl1","Eef1akmt1","Pdcd2l","Sh3bgrl2","Wls","Phtf2","Dcps","Hopx",
  "Ttc39b","Glipr2","Cdc42se2","Rhoh","Batf3","Itih5","Gbp3","Huwe1","Pus1","Smpdl3a",
  "Nfkbiz","Uck2","Twsg1","Lrrc8c","Spred2","Tnfrsf21","Snx14","Tlr7","Cdcp1","Galm",
  "Etfbkmt","Ptrh2","Ckap4","Lclat1","Drc1","Ahcyl","Plin2","Anxa4","Aplp1","Serpinc1",
  "Bcl2","Bcl2l1","Bmp2","Car2","Ccnd2","Ccnd3","Ccne1","Cd44","Cd48","Col6a1",
  "Ctla4","Plscr1","Ager","Tnfrsf18","Eno3","Fah","Flt3l","Gata1","Gpr83","Slc2a3",
  "Irf8","Cd79b","Igf1r","Igf2r","Il10","Il10ra","Il13","Il1r2","Il2ra","Il2rb",
  "Il3ra","Itga6","Itgav","Lif","Mxd1","Maff","Muc1","Myc","Pnp","Enpp1","Odc1",
  "P4ha1","Furin","Abcb1a","Pim1","Prkch","Prnp","Ptger2","Pth1r","Sell","Selp",
  "Spp1","Il1rl1","Tgm2","Xbp1","Etv4","Arl4a","Nt5e","Tnfrsf8"
)

# Filter available genes
available_il2_genes <- il2_stat5_genes[il2_stat5_genes %in% rownames(MDR_combined)]

# 2. Calculate Module Score
MDR_combined <- AddModuleScore(
  object = MDR_combined,
  features = list(available_il2_genes),
  name = "IL2_STAT5_Score"
)

p_il2_score <- SpatialFeaturePlot(MDR_combined, 
                                  features = "IL2_STAT5_Score1", 
                                  images = c("MDR_Con", "MDR_Cs"), 
                                  pt.size.factor = 2,
                                  min.cutoff = "q5", 
                                  max.cutoff = "q95") 

# Add Title and Labels
p_il2_score_final <- p_il2_score + 
  plot_annotation(title = "HALLMARK_IL2_STAT5_SIGNALING",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

p_il2_score_final

# 4. HALLMARK_BILE_ACID_METABOLISM
bilemetabolism_genes <- c(
  "Acsl1","Fdxr","Atxn1","Cyp7b1","Gclm","Hsd17b4","Cyp7a1","Pxmp2","Abca4","Npc1",
  "Amacr","Pipox","Gnmt","Pex7","Agxt","Soat2","Ch25h","Pex19","Slc22a18","Nr1i2",
  "Cyp8b1","Pex11a","Pex16","Dio2","Cyp46a1","Slc23a1","Klf1","Gnpat","Nr0b2","Slc27a2",
  "Slc27a5","Abcd1","Abcd3","Abcd2","Abca3","Abca8b","Hsd17b6","Slco1a4","Nr1h4","Bcar3",
  "Aldh1a1","Slc23a2","Prdx5","Aldh9a1","Abcg4","Aqp9","Bbox1","Isoc1","Lonp2","Retsat",
  "Pnpla8","Abcg8","Nudt12","Paox","Pex1","Optn","Efhc1","Acsl5","Pex13","Pex11g",
  "Slc35b2","Pex26","Crot","Dhcr24","Abca6","Gstk1","Fads1","Sult2b1","Slc29a1","Cyp39a1",
  "Mlycd","Hacl1","Fads2","Sult1b1","Hsd3b7","Pex12","Pecr","Hsd17b11","Akr1d1","Pex6",
  "Abca5","Abca9","Idi1","Tfcp2l1","Aldh8a1","Apoa1","Ar","Bmp6","Cat","Serpina6",
  "Cyp27a1","Phyh","Dio1","Gc","Hao1","Idh1","Idh2","Lck","Lipe","Nedd4","Pfkm",
  "Rbp1","Rxra","Rxrg","Scp2","Sod1","Ttr","Nr3c2","Ephx2","Abca2","Abca1"
)

# Filter available genes
available_ba_genes <- bilemetabolism_genes[bilemetabolism_genes %in% rownames(MDR_combined)]

# 2. Calculate Module Score
MDR_combined <- AddModuleScore(
  object = MDR_combined,
  features = list(available_ba_genes),
  name = "BA_metabolism_Score"
)

p_ba_score <- SpatialFeaturePlot(MDR_combined, 
                                  features = "BA_metabolism_Score1", 
                                  images = c("MDR_Con", "MDR_Cs"), 
                                  pt.size.factor = 2,
                                  min.cutoff = "q5", 
                                  max.cutoff = "q95") 

# Add Title and Labels
p_ba_score_final <- p_ba_score + 
  plot_annotation(title = "HALLMARK_BILE_ACID_METABOLISM",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

p_ba_score_final



# ==============================================================================
#GO
# 1. Extract Gene Lists from GO Results
# Ensure the Description strings match your ego_high@result exactly
ifnb_prod_genes <- ego_high@result %>%
  filter(Description == "interferon-beta production") %>%
  pull(geneID) %>%
  strsplit("/") %>% unlist()

ifnb_resp_genes <- ego_high@result %>%
  filter(Description == "response to interferon-beta") %>%
  pull(geneID) %>%
  strsplit("/") %>% unlist()

lipidstorage_genes <- ego_high@result %>%
  filter(Description == "lipid storage") %>%
  pull(geneID) %>%
  strsplit("/") %>% unlist()

lipidoxidation_genes <- ego_low@result %>%
  filter(Description == "lipid oxidation") %>%
  pull(geneID) %>%
  strsplit("/") %>% unlist()

# Filter for genes actually present in the dataset
ifnb_prod_genes <- ifnb_prod_genes[ifnb_prod_genes %in% rownames(MDR_combined)]
ifnb_resp_genes <- ifnb_resp_genes[ifnb_resp_genes %in% rownames(MDR_combined)]
lipid_storage_genes <- lipidstorage_genes[lipidstorage_genes %in% rownames(MDR_combined)]
lipid_oxidation_genes <- lipidoxidation_genes[lipidoxidation_genes %in% rownames(MDR_combined)]

# 2. Calculate Module Scores
# Scoring for 'Production' pathway
MDR_combined <- AddModuleScore(
  object = MDR_combined,
  features = list(ifnb_prod_genes),
  name = "IFNB_Production_Score"
)

# Scoring for 'Response' pathway
MDR_combined <- AddModuleScore(
  object = MDR_combined,
  features = list(ifnb_resp_genes),
  name = "IFNB_Response_Score"
)

# Scoring for lipid associated pathway
MDR_combined <- AddModuleScore(
  object = MDR_combined,
  features = list(lipid_storage_genes),
  name = "lipid_storage_Score"
)

MDR_combined <- AddModuleScore(
  object = MDR_combined,
  features = list(lipid_oxidation_genes),
  name = "lipid_oxidation_Score"
)

# 3. Spatial Visualization
# Note: Seurat appends '1' to the name provided in AddModuleScore
p_ifn_scores1 <- SpatialFeaturePlot(MDR_combined, 
                              features = "IFNB_Production_Score1", 
                              images = c("MDR_Con", "MDR_Cs"), 
                              pt.size.factor = 2,
                              min.cutoff = "q5", 
                              max.cutoff = "q95") 

p_ifn_scores2 <- SpatialFeaturePlot(MDR_combined, 
                                    features = "IFNB_Response_Score1", 
                                    images = c("MDR_Con", "MDR_Cs"), 
                                    pt.size.factor = 2,
                                    min.cutoff = "q5", 
                                    max.cutoff = "q95") 

p_lipid_scores1 <- SpatialFeaturePlot(MDR_combined, 
                                    features = "lipid_storage_Score1", 
                                    images = c("MDR_Con", "MDR_Cs"), 
                                    pt.size.factor = 2,
                                    min.cutoff = "q5", 
                                    max.cutoff = "q95") 

p_lipid_scores2 <- SpatialFeaturePlot(MDR_combined, 
                                      features = "lipid_oxidation_Score1", 
                                      images = c("MDR_Con", "MDR_Cs"), 
                                      pt.size.factor = 2,
                                      min.cutoff = "q5", 
                                      max.cutoff = "q95") 

p_ifn_scores1_plot <- p_ifn_scores1 + 
  plot_annotation(title = "GO:0032608 (IFN-beta production)",
                  theme = theme(plot.title = element_text(size = 18, 
                                                          #face = "bold", 
                                                          hjust = 0.5)))

p_ifn_scores2_plot <- p_ifn_scores2 + 
  plot_annotation(title = "GO:0035456 (Response to IFN-beta)",
                  theme = theme(plot.title = element_text(size = 18, 
                                                          #face = "bold", 
                                                          hjust = 0.5)))

p_lipid_scores1_plot <- p_lipid_scores1 + 
  plot_annotation(title = "GO:0019915 (Lipid storage)",
                  theme = theme(plot.title = element_text(size = 18, 
                                                          #face = "bold", 
                                                          hjust = 0.5)))

p_lipid_scores2_plot <- p_lipid_scores2 + 
  plot_annotation(title = "GO:0034440 (Lipid oxidation)",
                  theme = theme(plot.title = element_text(size = 18, 
                                                          #face = "bold", 
                                                          hjust = 0.5)))

p_ifn_scores1_plot
p_ifn_scores2_plot
p_lipid_scores1_plot
p_lipid_scores2_plot
