# 07_molecular_docking.R
# SYK抑制剂OXSI-2与PRKAA2分子对接准备 & 结果分析
library(bio3d)
library(ggplot2)
library(ggsci)

# 蛋白结构、配体文件本地路径，不上传pdb/sdf大文件
# protein <- read.pdb("your_local_path/PRKAA2.pdb")
# ligand <- read.pdb("your_local_path/OXSI-2.pdb")

# 对接后结果解析
# dock_result <- read.table("your_local_path/dock_out.txt", header=T)
# ggplot(dock_result, aes(x=binding_energy)) +
#   geom_histogram(fill="#E64B35CC", bins=15) +
#   theme_bw() + theme(panel.grid = element_blank()) +
#   xlab("Binding energy (kcal/mol)")
# ggsave("../figures/docking_energy.pdf", width=6, height=4)

cat("07 分子对接结果解析脚本\n")
