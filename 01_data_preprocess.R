# 01_data_preprocess.R
# Childhood asthma BALF/PBMC scRNA-seq full pipeline
# Seurat v5 | Harmony | scDblFinder | SingleR | Monocyte subclustering
# All raw sequencing data & large *.rds files are LOCAL ONLY, NOT uploaded to GitHub
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggforce)
library(ggsci)
library(harmony)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(celldex)
library(SingleR)
library(SummarizedExperiment)

# ====================== 1. Sample meta & path definition ======================
# !!! base_dir: raw fastq/mtx data path, local only, do NOT upload raw data to GitHub
base_dir <- "path/to/scRNAdata"

sample_ids <- c("sample1", "sample2", "sample3",   # Batch 1
                "sample4", "sample5", "sample6")           # Batch 2
group_info <- c("Healthy", "Healthy", "Asthma",
                "Asthma", "Asthma", "Healthy")
batch_info <- c("Batch1", "Batch1", "Batch1",
                "Batch2", "Batch2", "Batch2")
sample_names <- c("Control1", "Control2", "AS1", "AS2", "AS3", "Control3")

# ====================== 2. Read individual samples & create Seurat list ======================
sc_list <- list()
for (i in 1:length(sample_ids)) {
  sample <- sample_ids[i]
  path <- file.path(base_dir, sample)
  message("Loading: ", sample, " ...")
  
  mtx_path <- file.path(path, "matrix.mtx.gz")
  cells_path <- file.path(path, "barcodes.tsv.gz")
  features_path <- file.path(path, "features.tsv.gz")
  
  counts <- ReadMtx(
    mtx = mtx_path,
    cells = cells_path,
    features = features_path,
    feature.column = 1,
    cell.column = 1,
    skip.cell = 0,
    skip.feature = 0
  )
  
  obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)
  obj$Group <- group_info[i]
  obj$Batch <- batch_info[i]
  obj$SampleID <- sample_names[i]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  sc_list[[i]] <- obj
}

# Merge all samples
sc_combined <- merge(sc_list[[1]], y = sc_list[-1], add.cell.ids = sample_ids, project = "DNBC4")

# ====================== 3. Pre-filter QC visualization ======================
p_vln_pre <- VlnPlot(sc_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     ncol = 3, pt.size = 0, group.by = "SampleID") +
  theme_bw() + theme(panel.grid = element_blank()) + scale_fill_jco()

p_scatter1_pre <- FeatureScatter(sc_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                                 group.by = "Batch") +
  geom_hline(yintercept = c(500, 6000), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = c(1000, 40000), linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(title = "Before filtering: nFeature_RNA vs nCount_RNA") +
  theme_bw() + theme(panel.grid = element_blank())

p_scatter2_pre <- FeatureScatter(sc_combined, feature1 = "nCount_RNA", feature2 = "percent.mt",
                                 group.by = "Batch") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(title = "Before filtering: percent.mt vs nCount_RNA") +
  theme_bw() + theme(panel.grid = element_blank())

# ====================== 4. Cell filtering ======================
sc_combined_filtered <- subset(sc_combined,
                               subset = nFeature_RNA > 500 &
                                 nFeature_RNA < 6000 &
                                 nCount_RNA > 1000 &
                                 nCount_RNA < 40000 &
                                 percent.mt < 10)

# Post-filter QC plots
p_vln_post <- VlnPlot(sc_combined_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      ncol = 3, pt.size = 0, group.by = "SampleID") +
  theme_bw() + theme(panel.grid = element_blank()) + scale_fill_jco()

p_scatter1_post <- FeatureScatter(sc_combined_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                                  group.by = "Batch") +
  geom_hline(yintercept = c(500, 6000), linetype = "dashed") +
  geom_vline(xintercept = c(1000, 40000), linetype = "dashed") +
  theme_bw() + theme(panel.grid = element_blank())

p_scatter2_post <- FeatureScatter(sc_combined_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt",
                                  group.by = "Batch") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  theme_bw() + theme(panel.grid = element_blank())

# ====================== 5. Normalization & PCA + Harmony batch correction ======================
sc_combined <- NormalizeData(sc_combined)
sc_combined <- FindVariableFeatures(sc_combined, selection.method = "vst", nfeatures = 2000)
sc_combined <- ScaleData(sc_combined, vars.to.regress = "percent.mt")
sc_combined <- RunPCA(sc_combined)
sc_combined <- RunHarmony(sc_combined, group.by.vars = "Batch", plot_convergence = TRUE)

sc_combined <- RunUMAP(sc_combined, reduction = "harmony", dims = 1:30)
sc_combined <- FindNeighbors(sc_combined, reduction = "harmony", dims = 1:30)
sc_combined <- FindClusters(sc_combined, resolution = 0.5)

# ====================== 6. Find all cluster markers ======================
DefaultAssay(sc_combined) <- "RNA"
sc_combined <- JoinLayers(sc_combined, layers = "data", new.layer = "data")
Idents(sc_combined) <- "seurat_clusters"

markers <- FindAllMarkers(sc_combined,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

# Output marker tables to result folder
write.csv(markers, "../results_tables/all_markers.csv", row.names = FALSE)
write.csv(top10, "../results_tables/top10_markers.csv", row.names = FALSE)

marker_genes <- c("CD3D", "CD14", "MS4A1", "NKG7", "EPCAM")
p_dot <- DotPlot(sc_combined, features = marker_genes, group.by = "seurat_clusters") +
  RotatedAxis() + theme_bw() + theme(panel.grid = element_blank())

# ====================== 7. scDblFinder: doublet detection & remove doublets ======================
# Merge counts layer for doublet detection
sc_combined <- JoinLayers(sc_combined, layers = "counts", new.layer = "counts")
bp <- MulticoreParam(4)
doublet_labels_v5 <- rep(NA, ncol(sc_combined))
names(doublet_labels_v5) <- colnames(sc_combined)

for (s in unique(sc_combined$SampleID)) {
  message("Processing sample: ", s)
  cells_in_sample <- WhichCells(sc_combined, expression = SampleID == s)
  counts_mat <- GetAssayData(sc_combined, assay = "RNA", layer = "counts")[, cells_in_sample, drop = FALSE]
  sce <- SingleCellExperiment(list(counts = counts_mat))
  sce <- scDblFinder(sce,
                     samples = rep(s, ncol(sce)),
                     dbr = 0.03,
                     BPPARAM = bp)
  doublet_labels_v5[cells_in_sample] <- sce$scDblFinder.class
}

sc_combined$doublet_v5 <- doublet_labels_v5
sc_clean <- subset(sc_combined, subset = doublet_v5 == 1)

# ====================== 8. Re-run PCA / Harmony / UMAP after doublet removal ======================
DefaultAssay(sc_clean) <- "RNA"
sc_clean <- NormalizeData(sc_clean)
sc_clean <- FindVariableFeatures(sc_clean, selection.method = "vst", nfeatures = 2000)
sc_clean <- ScaleData(sc_clean, vars.to.regress = "percent.mt")
sc_clean <- RunPCA(sc_clean, npcs = 50)
ElbowPlot(sc_clean, ndims = 50)

sc_clean <- RunHarmony(sc_clean,
                       group.by.vars = "Batch",
                       reduction.use = "pca",
                       dims.use = 1:30,
                       plot_convergence = TRUE)

sc_clean <- RunUMAP(sc_clean, reduction = "harmony", dims = 1:30)
sc_clean <- FindNeighbors(sc_clean, reduction = "harmony", dims = 1:30)
sc_clean <- FindClusters(sc_clean, resolution = 0.5)

# Batch assessment plots
p_batch <- DimPlot(sc_clean, reduction = "umap", group.by = "Batch",
                   cols = c("Batch1" = "#E41A1C", "Batch2" = "#377EB8")) +
  ggtitle("UMAP after Harmony (by Batch)") + scale_y_reverse() +
  theme_bw() + theme(panel.grid = element_blank())

p_sample <- DimPlot(sc_clean, reduction = "umap", group.by = "SampleID",
                    label = FALSE, repel = TRUE) +
  ggtitle("UMAP after Harmony (by Sample)") + scale_y_reverse() +
  theme_bw() + theme(panel.grid = element_blank())

p_cluster <- DimPlot(sc_clean, reduction = "umap", group.by = "seurat_clusters",
                     label = TRUE, repel = TRUE) +
  ggtitle("UMAP after Harmony (Clusters)") + scale_y_reverse() +
  theme_bw() + theme(panel.grid = element_blank())

p_umap_cluster_label <- DimPlot(sc_clean, reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend() + scale_y_reverse() + theme_bw() + theme(panel.grid = element_blank())

# ====================== 9. Re-calculate markers for cleaned object & manual cell annotation ======================
Idents(sc_clean) <- "seurat_clusters"
markers_clean <- FindAllMarkers(sc_clean,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)

top10_clean <- markers_clean %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(markers_clean, "../results_tables/all_markers_clean.csv", row.names = F)
write.csv(top10_clean, "../results_tables/top10_markers_clean.csv", row.names = F)

# Manual cell type annotation
cluster_annotations <- c(
  "0" = "Naive CD4+ T cells",
  "1" = "NK cells",
  "2" = "Th cells",
  "3" = "CD8+ T cells",
  "4" = "Classical Monocytes",
  "5" = "Naive B cells",
  "6" = "Effector memory T cells",
  "7" = "Cytotoxic CD8+ T cells",
  "8" = "Regulatory T cell",
  "9" = "Th17 cells",
  "10" = "CD4+ T cells",
  "11" = "Non classical Monocytes",
  "12" = "cDC2",
  "13" = "Naive CD4+ T cells",
  "14" = "Megakaryocytes",
  "15" = "CD8+ T cells",
  "16" = "Proliferating cells",
  "17" = "Plasmacytoid dendritic cells",
  "18" = "Plasma cells",
  "19" = "Mast cells"
)

ann_vec <- cluster_annotations[as.character(sc_clean$seurat_clusters)]
names(ann_vec) <- NULL
sc_clean$celltype_clean <- ann_vec
sc_clean$celltype_clean <- factor(sc_clean$celltype_clean, levels = unique(cluster_annotations))

my_colors <- c(
  "#F9D9B5", "#C7E9C0", "#F4C6C6", "#B3CDE3", "#CCEBC5",
  "#E6F5C9", "#FDD0A2", "#FCCDE5", "#D9D9D9", "#E0F3F8",
  "#FFF2AE", "#D1E5F0", "#FDC086", "#FEF0D9", "#BEAED4",
  "#F4A582", "#ABD9E9", "#FEE0D2", "#B2DF8A", "#FFD92F"
)
my_colors <- my_colors[1:20]

p_umap_manual_anno <- DimPlot(sc_clean, group.by = "celltype_clean", label = TRUE, repel = TRUE,
                              label.size = 3) + NoLegend() +
  scale_color_manual(values = my_colors) +
  ggtitle("Single-Cell Landscape of PBMC (Asthma vs Healthy)") +
  theme_bw() + theme(panel.grid = element_blank())

p_feature_markers <- FeaturePlot(sc_clean, features = c("CD3D", "CD4", "CD8A", "NKG7", "MS4A1",
                                                        "CD14", "FCGR3A", "CLEC4C", "PPBP"),
                                 ncol = 3) & theme_bw() & theme(panel.grid = element_blank())

# ====================== 10. Monocyte subclustering ======================
Idents(sc_clean) <- "seurat_clusters"
mono_subset <- subset(sc_clean, idents = c("4", "11"))
DefaultAssay(mono_subset) <- "RNA"

mono_subset <- SCTransform(mono_subset, vars.to.regress = "percent.mt", verbose = TRUE)
mono_subset <- RunPCA(mono_subset, npcs = 30, verbose = TRUE)
mono_subset <- RunHarmony(mono_subset, group.by.vars = "Batch", dims.use = 1:30, plot_convergence = TRUE)
mono_subset <- RunUMAP(mono_subset, reduction = "harmony", dims = 1:20, return.model = TRUE)
mono_subset <- FindNeighbors(mono_subset, reduction = "harmony", dims = 1:20)
mono_subset <- FindClusters(mono_subset, resolution = 0.2)

Idents(mono_subset) <- "seurat_clusters"
mono_markers <- FindAllMarkers(mono_subset,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
top5_mono <- mono_markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
write.csv(mono_markers, "../results_tables/monocyte_subcluster_markers.csv", row.names = F)

p_mono_raw <- DimPlot(mono_subset, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Monocyte subclusters (resolution=0.2)") +
  theme_bw() + theme(panel.grid = element_blank())

# Monocyte subtype annotation
mono_subset$celltype_rough <- case_when(
  mono_subset$seurat_clusters == "0" ~ "Classical Monocytes",
  mono_subset$seurat_clusters == "1" ~ "Classical Monocytes",
  mono_subset$seurat_clusters == "2" ~ "Nonclassical Monocytes",
  mono_subset$seurat_clusters == "3" ~ "Intermediate Monocytes",
  mono_subset$seurat_clusters == "4" ~ "B cell contaminant",
  TRUE ~ "Unknown"
)

p_mono_anno_all <- DimPlot(mono_subset, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Monocyte subtypes") +
  scale_color_manual(
    values = c(
      "Classical Monocytes" = "#1f77b4",
      "Nonclassical Monocytes" = "#ff7f0e",
      "Intermediate Monocytes" = "#2ca02c",
      "B cell contaminant" = "#d62728",
      "Unknown" = "#949494"
    )
  ) + theme_bw() + theme(panel.grid = element_blank())

# Remove B cell contamination
Idents(mono_subset) <- "seurat_clusters"
mono_subset_filtered <- subset(mono_subset, idents = c("0", "1", "2", "3"))
mono_subset_filtered$celltype_rough_filtered <- case_when(
  mono_subset_filtered$seurat_clusters == "0" ~ "Classical Monocytes",
  mono_subset_filtered$seurat_clusters == "1" ~ "Classical Monocytes",
  mono_subset_filtered$seurat_clusters == "2" ~ "Nonclassical Monocytes",
  mono_subset_filtered$seurat_clusters == "3" ~ "Intermediate Monocytes",
  TRUE ~ "Unknown"
)

p_mono_anno_filtered <- DimPlot(mono_subset_filtered, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Monocyte subtypes (B cells removed)") +
  scale_color_manual(values = c("Classical Monocytes" = "#1f77b4",
                                "Nonclassical Monocytes" = "#ff7f0e",
                                "Intermediate Monocytes" = "#2ca02c",
                                "Unknown" = "gray70")) +
  theme_bw() + theme(panel.grid = element_blank())

# Map monocyte subtype back to main seurat object
sc_clean$mono_subtype <- NA
cells_to_update <- colnames(mono_subset_filtered)
annotations <- mono_subset_filtered$celltype_rough_filtered
sc_clean$mono_subtype[cells_to_update] <- annotations

# ====================== 11. SingleR supervised annotation (MonacoImmune reference) ======================
ref_monaco <- MonacoImmuneData()
query_expr <- GetAssayData(sc_clean, layer = "data")
query_expr <- as.matrix(query_expr)

pred_monaco_fine <- SingleR(test = query_expr,
                            ref = ref_monaco,
                            labels = ref_monaco$label.fine,
                            assay.type.test = 1,
                            assay.type.ref = 1)

sc_clean$SingleR_monaco_fine <- pred_monaco_fine$pruned.labels

p_umap_singler <- DimPlot(sc_clean, group.by = "SingleR_monaco_fine", label = TRUE, repel = TRUE,
                          label.size = 3) + NoLegend() +
  ggtitle("MonacoImmuneData Fine Labels") +
  theme_bw() + theme(panel.grid = element_blank())

# Extract DC cells summary
dc_cells <- subset(sc_clean, SingleR_monaco_fine == "Myeloid dendritic cells")
table_dc_group <- table(dc_cells$Group)
write.csv(as.data.frame(table_dc_group), "../results_tables/DC_cell_count_by_group.csv", row.names = F)

# ====================== NOTE: Save RDS locally, NOT upload to GitHub ======================
# saveRDS(sc_clean, file = "sc_clean_with_mono_subtypes.rds")
# saveRDS(mono_subset_filtered, file = "mono_subset_filtered.rds")

cat("01 Complete single-cell preprocessing, QC, doublet removal, cell annotation and monocyte subclustering script finished.\n")
cat("Result tables saved to ../results_tables/\n")
cat("Large rds objects are for local storage only.\n")
