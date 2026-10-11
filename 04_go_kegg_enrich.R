# 04_go_kegg_enrich.R
# GO KEGG富集分析，基于DE基因/DMR关联基因
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

# de_res <- read.csv("../results_tables/DE_genes_results.csv")
# sig_genes <- de_res %>% filter(p_val_adj <0.05) %>% pull(gene)

# GO富集
# go_enrich <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db,
#                       keyType = "SYMBOL", ont = "ALL", pAdjustMethod = "fdr", qvalueCutoff = 0.2)
# write.csv(go_enrich@result, "../results_tables/GO_enrich.csv", row.names = F)

# KEGG富集
# gene_entrez <- bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
# kegg_enrich <- enrichKEGG(gene = gene_entrez, organism = "hsa", pvalueCutoff = 0.2)
# write.csv(kegg_enrich@result, "../results_tables/KEGG_enrich.csv", row.names = F)

cat("04 GO/KEGG富集脚本，富集结果输出csv\n")
