rm(list = ls())

library(tidyverse)
library(openxlsx)

#读取数据
pca <- read.table("667XianSample_k31_sketch_Total0.1.PCA.eigenvec", header = F)
eigval <- read.table("667XianSample_k31_sketch_Total0.1.PCA.eigenval", header = F)

pcs <- paste0("PC", 1:nrow(eigval))
eigval[nrow(eigval),1] <- 0
percentage <- eigval$V1/sum(eigval$V1)*100
eigval_df <- as.data.frame(cbind(pcs, eigval[,1], percentage), stringsAsFactors = F) #最终获得每个主成分的方差解释率文件
names(eigval_df) <- c("PCs", "variance", "proportion")  #修改列名
eigval_df$variance <- as.numeric(eigval_df$variance)  #转化为数值数据
eigval_df$proportion <- as.numeric(eigval_df$proportion)  #转化为数值数据

#获取第一主成分的方差解释率，保留两位小数，round()函数是保留小数的
pc1_proportion <- paste0(round(eigval_df[1,3],2),"%")
pc2_proportion <- paste0(round(eigval_df[2,3],2),"%")
pc3_proportion <- paste0(round(eigval_df[3,3],2),"%")
pc4_proportion <- paste0(round(eigval_df[4,3],2),"%")
pc5_proportion <- paste0(round(eigval_df[5,3],2),"%")
pc6_proportion <- paste0(round(eigval_df[6,3],2),"%")
pc7_proportion <- paste0(round(eigval_df[7,3],2),"%")
pc8_proportion <- paste0(round(eigval_df[8,3],2),"%")
pc9_proportion <- paste0(round(eigval_df[9,3],2),"%")
pc10_proportion <- paste0(round(eigval_df[10,3],2),"%")

sample <- read.xlsx("667XianSampleGroup.xlsx", sheet = "Sheet1")  #读取样本的分群体文件
sample <- data.frame(V2=sample[,1],Population=sample[,2])
data <- left_join(sample,pca[,2:12],by="V2") 
#将pca[,1:4]数据和sample数据框合并，左连接，并以“V1”数据作为链接标准
#之所以取pca[,1:4],是因为我们只要前两个PC，如果你要前三个，则是pca[,1:5]
colnames(data) <- c("SampleID","Group","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
data$Group <- factor(data$Group,levels = c("XianLandrace1","XianLandrace2","XianCultivar1","XianCultivar2","Outgroup")) #将分组转化为因子变量
data <- data[!is.na(data$Group),,]


ggplot(data,aes(PC2,PC3))+
  geom_point(aes(color=Group), size=0.4)+
  stat_ellipse(aes(color=Group),level = 0.95, show.legend = FALSE, linewidth = 0.35)+
  scale_color_manual(values = c(XianLandrace1="#C5000B",XianLandrace2="#FF950E",XianCultivar2="#83CAFF",XianCultivar1="#579D1C"))+  #如果有4个分组的话，就要四个颜色
  labs(x=paste0("PC2(",pc2_proportion,")"),y=paste0("PC3 (",pc3_proportion,")"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  #移除主要网格线
        panel.grid.minor = element_blank(),  #移除次要网格线
        axis.text.x = element_text(family = "sans",size = 10,colour = "black",face = "plain",angle = 0),  #x轴标签字体大小及样式，可以用axis.title.x指定x轴标题字体、位置、颜色等
        axis.title.x = element_text(size = 10,margin = margin(t = 5)),  #轴标题大小及距离x轴标签距离,在 margin 函数中，t 参数表示上边距（距离 x 轴标题到 x 轴标签的距离），r、b 和 l 分别表示右边距、下边距和左边距
        axis.text.y = element_text(family = "sans",size = 10,colour = "black",face = "plain",angle = 0),  #y轴标签字体大小及样式
        axis.title.y = element_text(size = 10,margin = margin(t = 5)),
        legend.position = "right",  #可以是 "left", "right", "top", "bottom" 或 c(x, y)
        legend.background = element_rect(fill = NA, color = NA), # 图例背景  color控制边框颜色
        legend.title = element_blank(), # 图例标题
        legend.text = element_text(size = 6),
        legend.key.height = unit(10, "pt"))


ggsave(filename = "667XianSample.PCA.Kmer.PC2PC3.pdf",width = 9.5,height = 5.5,units = "cm")



