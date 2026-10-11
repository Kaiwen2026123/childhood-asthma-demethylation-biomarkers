# 06_machine_learning.R
# 7-gene signature 诊断预测模型，儿童哮喘生物标志物
library(caret)
library(randomForest)
library(pROC)
library(dplyr)

# 表达矩阵+临床样本信息本地读取
# expr_matrix <- read.csv("your_local_path/sample_gene_expression.csv")
# meta <- read.csv("../results_tables/metadata.csv")

# 7基因签名
# sig7_genes <- c("HLA-DPB1","CIITA","IRF7","CXCL10","ANO1","SYK","PRKAA2")
# dat <- expr_matrix[,colnames(expr_matrix) %in% c("group", sig7_genes)]

# 训练集测试集划分
# set.seed(123)
# trainIndex <- createDataPartition(dat$group, p=0.7, list=F)
# train <- dat[trainIndex,]
# test <- dat[-trainIndex,]

# 随机森林模型
# rf_ctrl <- trainControl(method="cv", number=5, classProbs = T, summaryFunction = twoClassSummary)
# rf_model <- train(group ~ ., data=train, method="rf", trControl=rf_ctrl)

# 预测+ROC
# pred_prob <- predict(rf_model, newdata=test, type="prob")[,2]
# roc_obj <- roc(response = test$group, predictor = pred_prob)
# auc(roc_obj)

# 保存模型结果表
# res_df <- data.frame(gene=sig7_genes, importance=varImp(rf_model)$importance[,1])
# write.csv(res_df, "../results_tables/signature_7gene_result.csv", row.names=F)

cat("06 7基因signature机器学习模型脚本完成\n")
