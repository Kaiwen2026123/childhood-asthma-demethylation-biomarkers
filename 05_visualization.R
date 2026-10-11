# 05_visualization.R
# 论文可视化：UMAP、FeaturePlot、DE热图、火山图
library(Seurat)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(EnhancedVolcano)

# seurat_integrated <- readRDS("your_local_path/seurat_integrated.rds")
# de_res <- read.csv("../results_tables/DE_genes_results.csv")

# UMAP图
# DimPlot(seurat_integrated, reduction = "umap", group.by = "cell_type", label=T) +
#   theme_bw() + theme(panel.grid = element_blank()) + scale_color_npg()
# ggsave("../figures/Figure2_UMAP.pdf", width=8, height=6)

# 火山图
# EnhancedVolcano(de_res, lab = de_res$gene, x = "avg_log2FC", y = "p_val_adj",
#                 pCutoff = 0.05, FCcutoff = 0.25,
#                 col = c("grey30", "grey50", "blue", "red3"),
#                 labSize = 3) + theme_bw() + theme(panel.grid = element_blank())
# ggsave("../figures/Figure_volcano.pdf", width=9, height=7)

# 差异基因热图
# top_de_genes <- de_res %>% slice_min(p_val_adj, n=50) %>% pull(gene)
# DoHeatmap(seurat_integrated, features = top_de_genes) +
#   theme_bw() + theme(panel.grid = element_blank())
# ggsave("../figures/Figure_heatmap.pdf", width=10, height=7)

cat("05 可视化绘图脚本，输出pdf矢量图到figures文件夹\n")
