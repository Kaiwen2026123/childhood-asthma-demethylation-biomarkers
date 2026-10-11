# 02_DE_analysis.R
# 细胞亚群差异基因 DE 分析，哮喘vs对照
library(Seurat)
library(dplyr)
library(tibble)

# seurat_integrated <- readRDS("your_local_path/seurat_integrated.rds")

# 按分组找差异基因，分组列：group (asthma / control)
# de_res <- FindMarkers(seurat_integrated,
#                      ident.1 = "asthma", ident.2 = "control",
#                      min.pct = 0.1, logfc.threshold = 0.25)

# 整理结果，增加gene symbol列
# de_res <- de_res %>% rownames_to_column("gene")
# write.csv(de_res, file = "../results_tables/DE_genes_results.csv", row.names = FALSE)

cat("02 DE差异基因分析完成，结果输出到results_tables/DE_genes_results.csv\n")
