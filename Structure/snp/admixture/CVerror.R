rm(list = ls())

# 加载必要的包
library(ggplot2)
library(tidyr)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(minpack.lm)
library(broom)


mydata <- read.xlsx("667CVerror.SNP.xlsx",sheet = "Sheet1")

# 绘制折线图
ggplot(mydata, aes(x = K, y = Cverror)) +
  geom_line(linewidth = 0.4,color = "black") +
  geom_point(aes(x = K, y = Cverror), color = "black", size = 0.6) +
  labs(x = "Subgroup number", y = "CV error") +
  #geom_point(aes(x = 8, y = 0.40624), color = "red", size = 0.8) +
  #annotate("segment", x = 1, xend = 8, y = 0.40624, yend = 0.40624,linetype = "dashed", linewidth = 0.5, color = "red") +  # 横线
  #annotate("segment", x = 8, xend = 8, y = 0.2, yend = 0.40624,linetype = "dashed", linewidth = 0.5, color = "red") +  # 竖线
  #annotate("text", x = 5, y = 0.35, label = "K = 8\nCV = 0.40624", color = "red", size = 2.5, hjust = 0)+
  #ylim(0.2,0.55)+
  theme_bw() +  #使绘图背景为白色
  theme(panel.grid.major = element_blank(),  #移除主要网格线
        panel.grid.minor = element_blank(),  #移除次要网格线
        axis.text.x = element_text(family = "sans",size = 12,colour = "black",face = "plain",angle = 0),  #x轴标签字体大小及样式，可以用axis.title.x指定x轴标题字体、位置、颜色等
        axis.title.x = element_text(size = 12,margin = margin(t = 5)),  #轴标题大小及距离x轴标签距离,在 margin 函数中，t 参数表示上边距（距离 x 轴标题到 x 轴标签的距离），r、b 和 l 分别表示右边距、下边距和左边距
        axis.text.y = element_text(family = "sans",size = 12,colour = "black",face = "plain",angle = 0),  #y轴标签字体大小及样式
        axis.title.y = element_text(size = 12,margin = margin(r = 5)))

ggsave(filename = "667Sample.CVerror.SNP.pdf",units = "cm",width = 7,height = 6)
