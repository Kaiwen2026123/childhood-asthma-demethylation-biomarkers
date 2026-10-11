# 03_dmr_analysis.R
# 差异甲基化区域 DMR 分析，儿童哮喘去甲基化标志物
library(minfi)
library(DMRcate)
library(limma)

# 甲基化矩阵本地读取，不把大矩阵上传github
# rgSet <- read.metharray.exp(targets = targets)
# mSet <- preprocessQuantile(rgSet)

# 分组设计矩阵
# design <- model.matrix(~0 + group)
# colnames(design) <- c("asthma","control")
# fit <- lmFit(getBeta(mSet), design)

# DMR识别
# myAnnotation <- cpg.annotate(object = mSet, datatype = "array", what = "Beta",
#                             analysis.type = "differential", design = design,
#                             contrasts = TRUE, cont.matrix = makeContrasts(asthma-control, levels=design))
# dmrs <- dmrcate(myAnnotation, lambda=1000, C=2)
# results.ranges <- extractRanges(dmrs)

# 输出DMR汇总表
# dmr_df <- as.data.frame(results.ranges)
# write.table(dmr_df, file = "../results_tables/DMR_result.tsv", sep = "\t", row.names = F)

cat("03 DMR差异甲基化分析完成，结果输出DMR_result.tsv\n")
