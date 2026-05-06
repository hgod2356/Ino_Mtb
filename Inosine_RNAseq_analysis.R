################################################################################
#
# RNA-seq Analysis Pipeline: Effect of Inosine on Mtb (H37Rv) Infection
# ------------------------------------------------------------------------------
# Description : Single-end RNA-seq pipeline from raw FASTQ to differential
#               expression and pathway enrichment analysis.
# Workflow    : (1) Genome indexing      -> Rsubread::buildindex
#               (2) Read alignment        -> Rsubread::align
#               (3) Gene-level counting   -> Rsubread::featureCounts
#               (4) Quality control       -> edgeR::plotMDS
#               (5) DEG analysis          -> edgeR (TMM + GLM/LRT)
#               (6) Heatmap visualization -> ComplexHeatmap
#               (7) Pathway enrichment    -> clusterProfiler::gseKEGG
# Organism    : Mus musculus (GRCm39, Ensembl release 112)
# Design      : 3 vs 3 (Mtb vs Mtb + Inosine), filtered to 2 vs 2 after QC
#
################################################################################


# ============================================================================
# 0. Dependencies
# ============================================================================

# Install (run once):
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("Rsubread", "edgeR", "DESeq2", "tximport", "ShortRead",
#                        "org.Mm.eg.db", "AnnotationDbi", "ComplexHeatmap",
#                        "circlize", "clusterProfiler", "enrichplot", "msigdbr"))
# install.packages(c("ggplot2"))

suppressPackageStartupMessages({
  library(Rsubread)
  library(edgeR)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ComplexHeatmap)
  library(circlize)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(ggplot2)
})


# ============================================================================
# 1. Configuration  (edit BASE_DIR only)
# ============================================================================

BASE_DIR <- "/Users/darrenjung/Documents/3 Projects/Inosine"

# Sub-directories (auto-derived from BASE_DIR)
fastq_dir   <- file.path(BASE_DIR, "fastq")
output_dir  <- file.path(BASE_DIR, "output_bam")
ref_dir     <- file.path(BASE_DIR, "ref")
index_dir   <- file.path(BASE_DIR, "genome_index")
result_dir  <- file.path(BASE_DIR, "results")

# Reference files
genome_fa  <- file.path(ref_dir, "Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz")
gtf_file   <- file.path(ref_dir, "Mus_musculus.GRCm39.112.gtf.gz")
index_path <- file.path(index_dir, "subread_index")

# Create output folders if missing
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(index_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

# Sample names (must match FASTQ filenames without ".fastq.gz")
sample_names <- c("control-1_2",   "control-2_2",   "control-3_2",
                  "Experient-1_2", "Experient-2_2", "Experient-3_2")

# Group labels (length must equal length(sample_names))
group_labels <- c(rep("H37Rv",     3),
                  rep("H37Rv_INO", 3))


# ============================================================================
# 2. Build genome index  (run once)
# ============================================================================

buildindex(
  basename  = index_path,
  reference = genome_fa
)


# ============================================================================
# 3. Read alignment (single-end)
# ============================================================================

for (sample in sample_names) {
  readfile1  <- file.path(fastq_dir,  paste0(sample, ".fastq.gz"))
  output_bam <- file.path(output_dir, paste0(sample, ".bam"))

  align(
    index                  = index_path,
    readfile1              = readfile1,
    output_file            = output_bam,
    type                   = "rna",
    nthreads               = 4,
    unique                 = TRUE,
    sortReadsByCoordinates = TRUE
  )

  cat("Finished aligning:", sample, "\n")
}


# ============================================================================
# 4. Gene-level read counting (featureCounts)
# ============================================================================

bam_files <- file.path(output_dir, paste0(sample_names, ".bam"))

fc <- featureCounts(
  files               = bam_files,
  annot.ext           = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType     = "exon",
  GTF.attrType        = "gene_id",
  isPairedEnd         = FALSE,
  useMetaFeatures     = TRUE,
  nthreads            = 4
)

count_matrix <- fc$counts
colnames(count_matrix) <- sample_names

write.csv(count_matrix,
          file = file.path(result_dir, "count_matrix_raw.csv"))


# ============================================================================
# 5. Sample metadata
# ============================================================================

metadata <- data.frame(
  row.names = sample_names,
  condition = group_labels
)


# ============================================================================
# 6. Differential expression analysis (edgeR, full 3 vs 3)
# ============================================================================

group <- factor(metadata$condition)

# Build DGEList and filter low-expression genes
dge  <- DGEList(counts = count_matrix, group = group)
keep <- filterByExpr(dge)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes retained after filtering:", nrow(dge$counts), "\n")

# TMM normalization and dispersion estimation
dge    <- calcNormFactors(dge)
design <- model.matrix(~ group)
dge    <- estimateDisp(dge, design)

# GLM fit + likelihood ratio test
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

# Extract full results table
results   <- topTags(lrt, n = Inf)
deg_table <- results$table
deg_table$ensembl_id <- rownames(deg_table)

# Map Ensembl IDs to gene symbols
deg_table$gene_name <- mapIds(
  org.Mm.eg.db,
  keys      = deg_table$ensembl_id,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

write.csv(deg_table,
          file      = file.path(result_dir, "DEG_results_3vs3.csv"),
          row.names = FALSE)


# ============================================================================
# 7. Quality control: MDS plot
# ============================================================================

mds <- plotMDS(dge, plot = FALSE)

mds_df <- data.frame(
  x      = mds$x,
  y      = mds$y,
  sample = colnames(dge),
  group  = group_labels
)

p_mds <- ggplot(mds_df, aes(x = x, y = y, color = group, label = sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3.5) +
  scale_color_manual(values = c("H37Rv"     = "#4472C4",
                                "H37Rv_INO" = "#ED7D31")) +
  labs(title = "MDS Plot",
       x     = "Leading logFC dim 1",
       y     = "Leading logFC dim 2",
       color = "Group") +
  theme_bw()

ggsave(file.path(result_dir, "MDS_plot.pdf"),
       p_mds, width = 6, height = 5)


# ============================================================================
# 8. Outlier removal and re-analysis (2 vs 2)
# ============================================================================
# Based on MDS, control-1_2 and Experient-3_2 were excluded as outliers.

set.seed(100)

keep_samples <- c("control-2_2", "control-3_2",
                  "Experient-1_2", "Experient-2_2")

dge_filtered <- dge[, keep_samples]
dge_filtered$samples$group <- factor(c("H37Rv", "H37Rv",
                                       "H37Rv_INO", "H37Rv_INO"))

# Re-filter and re-normalize on the smaller design
keep         <- filterByExpr(dge_filtered)
dge_filtered <- dge_filtered[keep, , keep.lib.sizes = FALSE]
cat("Genes retained after refiltering (2v2):",
    nrow(dge_filtered$counts), "\n")

dge_filtered <- calcNormFactors(dge_filtered)

design_f     <- model.matrix(~ dge_filtered$samples$group)
dge_filtered <- estimateDisp(dge_filtered, design_f)

fit_f <- glmFit(dge_filtered, design_f)
lrt_f <- glmLRT(fit_f)

results_filtered <- topTags(lrt_f, n = Inf)
deg_table_f <- results_filtered$table
deg_table_f$ensembl_id <- rownames(deg_table_f)
deg_table_f$gene_name  <- mapIds(
  org.Mm.eg.db,
  keys      = deg_table_f$ensembl_id,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

write.csv(deg_table_f,
          file      = file.path(result_dir, "DEG_results_2vs2.csv"),
          row.names = FALSE)


# ============================================================================
# 9. Heatmap of cytokine / chemokine genes
# ============================================================================

# Compute log-CPM and convert row names to gene symbols
log_cpm <- cpm(dge, log = TRUE)
rownames(log_cpm) <- deg_table$gene_name[
  match(rownames(log_cpm), deg_table$ensembl_id)
]

# Target gene panel (cytokines / chemokines of interest)
target_genes <- c("Cxcl1", "Cxcl12", "Cxcl17", "Cxcl15",
                  "Ccl1",  "Cx3cl1", "Tnf",    "Ccl2",
                  "Ccl7",  "Cxcl14", "Il21",   "Ccl20",
                  "Tnfsf15")

# Subset matrix to target genes (keep input order) and 2v2 samples
mat1 <- log_cpm[rownames(log_cpm) %in% target_genes, ]
mat1 <- mat1[match(target_genes[target_genes %in% rownames(mat1)],
                   rownames(mat1)), ]
mat1 <- mat1[, keep_samples]

# Row-wise Z-score
mat1_z <- t(scale(t(mat1)))

col_anno <- HeatmapAnnotation(
  Group = c(rep("Mtb", 2), rep("Mtb + INO", 2)),
  col   = list(Group = c("Mtb"       = "#4472C4",
                         "Mtb + INO" = "#ED7D31")),
  show_annotation_name = FALSE
)

ht_cytokine <- Heatmap(
  mat1_z,
  name              = "Z",
  col               = colorRamp2(c(-2, 0, 2),
                                 c("blue", "white", "red")),
  top_annotation    = col_anno,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  show_column_names = FALSE,
  row_names_side    = "right",
  row_names_gp      = gpar(fontsize = 9),
  column_title      = "Cytokine & Chemokine",
  column_title_gp   = gpar(fontface = "bold")
)

pdf(file.path(result_dir, "Heatmap_cytokine_chemokine.pdf"),
    width = 5, height = 5)
draw(ht_cytokine)
dev.off()


# ============================================================================
# 10. KEGG pathway enrichment (GSEA, on 2 vs 2 results)
# ============================================================================

# Map Ensembl -> Entrez (required by gseKEGG)
entrez_map <- mapIds(
  org.Mm.eg.db,
  keys      = deg_table_f$ensembl_id,
  column    = "ENTREZID",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# Build ranked gene list (logFC), drop NAs/duplicates, sort descending
gene_list <- deg_table_f$logFC
names(gene_list) <- entrez_map
gene_list <- gene_list[!is.na(names(gene_list))]
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Ranked gene list size:", length(gene_list), "\n")

# Run KEGG GSEA
gsea_kegg <- gseKEGG(
  geneList     = gene_list,
  organism     = "mmu",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE,
  seed         = 42
)

gsea_df <- as.data.frame(gsea_kegg)
write.csv(gsea_df,
          file      = file.path(result_dir, "GSEA_KEGG_results.csv"),
          row.names = FALSE)


# ============================================================================
# 11. GSEA bubble-bar plot for selected pathways
# ============================================================================

target_pathways <- c(
  "mmu00190",   # Oxidative phosphorylation
  "mmu04660",   # T cell receptor signaling pathway
  "mmu04060",   # Cytokine-cytokine receptor interaction
  "mmu04976",   # Bile secretion
  "mmu04151",   # PI3K-Akt signaling pathway
  "mmu04020"    # Calcium signaling pathway
)

gsea_plot_df <- gsea_df[gsea_df$ID %in% target_pathways,
                        c("ID", "Description", "NES", "p.adjust")]
gsea_plot_df <- gsea_plot_df[order(gsea_plot_df$NES, decreasing = TRUE), ]
gsea_plot_df$Description  <- factor(gsea_plot_df$Description,
                                    levels = rev(gsea_plot_df$Description))
gsea_plot_df$neg_log_padj <- -log10(gsea_plot_df$p.adjust)
gsea_plot_df$group        <- ifelse(gsea_plot_df$NES > 0, "pos", "neg")

# Build per-pathway colour map: red gradient (positive NES),
# blue gradient (negative NES); intensity scales with |NES|
pos_df <- gsea_plot_df[gsea_plot_df$NES > 0, ]
neg_df <- gsea_plot_df[gsea_plot_df$NES < 0, ]
pos_df <- pos_df[order(pos_df$NES), ]
neg_df <- neg_df[order(neg_df$NES, decreasing = TRUE), ]

pos_colors <- colorRampPalette(c("#FFAAAA", "#CC0000"))(nrow(pos_df))
neg_colors <- colorRampPalette(c("#AACCFF", "#003399"))(nrow(neg_df))

color_map <- c(
  setNames(pos_colors, pos_df$ID),
  setNames(neg_colors, neg_df$ID)
)

p_gsea <- ggplot(gsea_plot_df, aes(x = NES, y = Description)) +

  geom_col(aes(fill = ID), width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = color_map) +

  # KEGG ID inside the bar
  geom_text(aes(label = paste0("KEGG:", ID), x = NES / 2),
            color = "white", size = 3.2, fontface = "bold") +

  # NES value outside the bubble
  geom_text(aes(label = round(NES, 2),
                x = ifelse(NES > 0, NES + 0.55, NES - 0.55)),
            hjust = 0.5, size = 3.5, fontface = "bold") +

  # Pathway description on the side opposite to the bar
  geom_text(aes(label = Description,
                x = ifelse(NES > 0, -0.08, 0.08)),
            hjust = ifelse(gsea_plot_df$NES > 0, 1.05, -0.05),
            size = 3.3, color = "black") +

  # -log10(p.adjust) bubble at the bar tip
  geom_point(aes(x = ifelse(NES > 0, NES + 0.3, NES - 0.3),
                 size = neg_log_padj),
             color = "#CCCC00", shape = 16, alpha = 0.9) +

  scale_size_continuous(name   = "-log(p.adjust)",
                        range  = c(2, 10),
                        breaks = c(0, 2, 5, 10),
                        labels = c("0", "2", "5", "10")) +
  scale_x_continuous(limits = c(-3.2, 3.2),
                     breaks = seq(-3, 3, 1)) +

  geom_vline(xintercept = 0, color = "black", linewidth = 0.7) +

  annotate("text", x = -3.0,
           y = length(levels(gsea_plot_df$Description)) + 0.3,
           label = "H37Rv",
           color = "#003399", fontface = "bold", size = 4.5) +
  annotate("text", x = 3.0,
           y = length(levels(gsea_plot_df$Description)) + 0.3,
           label = "H37Rv + INO",
           color = "#CC0000", fontface = "bold", size = 4.5) +

  labs(title = "GSEA Analysis",
       x     = "Normalized enrichment score",
       y     = NULL) +

  theme_minimal() +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.x        = element_text(size = 10),
    axis.title.x       = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.title       = element_text(size = 9),
    legend.key.size    = unit(0.5, "cm")
  )

ggsave(file.path(result_dir, "GSEA_KEGG_bubble.pdf"),
       p_gsea, width = 9, height = 5)


# ============================================================================
# 12. Export raw count matrix for GEO submission
# ============================================================================
# After QC the four samples used in downstream analysis are exported as a
# tab-delimited text file in GEO's preferred format.

count_matrix_geo <- count_matrix[, keep_samples]

geo_dir <- file.path(BASE_DIR, "GEO_submission")
dir.create(geo_dir, showWarnings = FALSE, recursive = TRUE)

write.table(
  count_matrix_geo,
  file      = file.path(geo_dir, "raw_count_matrix.txt"),
  sep       = "\t",
  quote     = FALSE,
  row.names = TRUE,
  col.names = NA
)

cat("Pipeline finished. Results written to:", result_dir, "\n")

