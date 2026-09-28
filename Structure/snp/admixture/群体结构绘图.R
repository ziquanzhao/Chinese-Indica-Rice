rm(list = ls())
library(ggplot2)
library(reshape2)
library(patchwork)
library(openxlsx)

#------------------------------
# ① 定义绘图函数
#------------------------------
StructurePlot <- function(mydata, Q){
  myplot <- ggplot(data = mydata, aes(x = Sample, y = Proportion, fill = Group, colour = Group)) +
    geom_bar(stat = "identity", position = "stack", width = 1, linewidth = 0) +
    scale_color_manual(values = c("#C5000B","#004586","#579D1C","#FF950E","#7E0021","#83CAFF","#314004","#AECF00","#4B1F6F","#DC3912","#0084D1","#FFFFCC","#FF5500","#654CFF","#FF99BF"))+
    scale_fill_manual(values = c("#C5000B","#004586","#579D1C","#FF950E","#7E0021","#83CAFF","#314004","#AECF00","#4B1F6F","#DC3912","#0084D1","#FFFFCC","#FF5500","#654CFF","#FF99BF"))+
    theme_bw() +
    theme(legend.position = "none",
          plot.margin=unit(c(0,0,0,0),"mm"),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
  return(myplot)
}

# -----------------------------
# ② 自动获取样本ID（来自fam）
# -----------------------------
prefix <- "667XianSample.PASS.missing0.9.maf0.01.Biallelic.OnlyGT.WildOutgroup.OnlySNP.LDfilter"
fam_file <- paste0(prefix, ".fam")
fam_data <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
samples_fam <- fam_data$V2

# -----------------------------
# ③ 读取外部样本顺序文件（可选）
# -----------------------------
sample_order_file <- "667XianSample.SampleOrder.NoWild.Reroot.xlsx"  # ← 若存在则按此顺序绘图
if (file.exists(sample_order_file)) {
  message("✅ 检测到外部样本顺序文件: ", sample_order_file)
  samples_order <- read.xlsx(sample_order_file,sheet = "Sheet1",colNames = F)[,1]
  samples_order <- rev(samples_order)
} else {
  message("⚠️ 未检测到外部样本顺序文件，使用 .fam 文件默认顺序")
  samples_order <- samples_fam
}

# -----------------------------
# ④ 手动指定 K 值范围
# -----------------------------
K_values <- 2:15  # 你想画的 K 值范围

plot_list <- list()

# -----------------------------
# ⑤ 读取并合并Q矩阵
# -----------------------------
for (K in K_values) {
  Q_file <- paste0(prefix, ".", K, ".Q")
  if (!file.exists(Q_file)) {
    warning(paste("⚠️ 文件不存在:", Q_file))
    next
  }
  
  Q_data <- read.table(Q_file, header = FALSE)
  colnames(Q_data) <- paste0("Cluster", 1:K)
  Q_data$Sample <- samples_fam  # 与fam样本对应
  
  
  # 转为长格式
  df_long <- reshape2::melt(Q_data, id.vars = "Sample",
                            variable.name = "Group", value.name = "Proportion")
  
  # 使用指定顺序
  df_long$Sample <- factor(df_long$Sample, levels = samples_order)
  
  plot_list[[paste0("K", K)]] <- StructurePlot(df_long, K)
}

# -----------------------------
# ⑥ 最下方样本标签
# -----------------------------
SampleLabel <- ggplot(
  data.frame(Sample = factor(samples_order, levels = samples_order)),
  aes(x = Sample, y = 0)) +
  geom_text(aes(label = Sample), angle = 90, hjust = 1, vjust = 0.5, size = 1) +
  scale_x_discrete(drop = FALSE, limits = samples_order, expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = unit(c(0,0,0,0), "mm"))

# -----------------------------
# ⑦ 拼接图层并输出
# -----------------------------
FinalPlot <- wrap_plots(plot_list, ncol = 1) / SampleLabel + plot_layout(heights = c(rep(1, length(plot_list)), 0.6))

ggsave("667XianSample.Admixture.K2-15.NoWild.Reroot.pdf", FinalPlot, width = 35, height = 1.5 * (length(plot_list) + 1))
