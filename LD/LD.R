rm(list = ls())

library(ggplot2)

r2_mean <- function(lddata,startcut=100,cutpoint=1000,bin1=10,bin2=100) {
  #startcut=100,意思是前100行不进行求平均运算
  #cutpoint是分割点，意味着101~1000（cutpoint=1000）行采用每100（bin1=100）行求平均，从2001行开始往后每1000（bin2=1000）行求平均
  lddata_start <- lddata[1:startcut,]
  colnames(lddata_start) <- c("Distance","r2")
  
  lddata_head <- lddata[startcut+1:cutpoint,]  #截取101~2000（cutpoint=2000）行
  lddata_head$V3 <- ceiling(lddata_head$V1 / bin1)# 创建分组因子，每100（bin1=100）行一组
  lddata_head_mean <- aggregate(V2 ~ V3,data = lddata_head,FUN = mean) # 计算每组的均值
  result_lddata_head <- data.frame(Distance = lddata_head_mean$V3 * bin1, r2 = lddata_head_mean$V2) # 创建新的数据框，Dist为每组的最大值，r2为均值
  
  #对于从（cutpoint=1000）往后的行也是如此
  lddata_other <- lddata[cutpoint+1:nrow(lddata),]
  lddata_other$V3 <- ceiling(lddata_other$V1 / bin2)
  lddata_other_mean <- aggregate(V2 ~ V3,data = lddata_other,FUN = mean)
  result_lddata_other <- data.frame(Distance = lddata_other_mean$V3 * bin2, r2 = lddata_other_mean$V2)
  
  Final_result <- rbind(lddata_start,result_lddata_head,result_lddata_other)
  return(Final_result)
}


#ld_group1 <- read.table("GengCultivar.LDdecay.final.stat",header = F,sep = "\t")[,1:2]
#ld_group2 <- read.table("GengLandrace1.LDdecay.final.stat",header = F,sep = "\t")[,1:2]
#ld_group3 <- read.table("GengLandrace2.LDdecay.final.stat",header = F,sep = "\t")[,1:2]
ld_group4 <- read.table("XianCultivar1.PopLDdecay.stat",header = F,sep = "\t")[,1:2]
ld_group5 <- read.table("XianCultivar2.PopLDdecay.stat",header = F,sep = "\t")[,1:2]
ld_group6 <- read.table("XianLandrace1.PopLDdecay.stat",header = F,sep = "\t")[,1:2]
ld_group7 <- read.table("XianLandrace2.PopLDdecay.stat",header = F,sep = "\t")[,1:2]


#ld_group1_mean <- r2_mean(lddata = ld_group1)
#ld_group2_mean <- r2_mean(lddata = ld_group2)
#ld_group3_mean <- r2_mean(lddata = ld_group3)
ld_group4_mean <- r2_mean(lddata = ld_group4)
ld_group5_mean <- r2_mean(lddata = ld_group5)
ld_group6_mean <- r2_mean(lddata = ld_group6)
ld_group7_mean <- r2_mean(lddata = ld_group7)


ld_data <- rbind(ld_group4_mean,ld_group5_mean,ld_group6_mean,ld_group7_mean)
ld_data$Subgroup <- rep(x = c("XianCultivar1","XianCultivar2","XianLandrace1","XianLandrace2"), each = nrow(ld_data)/4) 
#创造一列分亚群，分几个亚群就除以几


# 创建LD衰减图
ggplot(data = ld_data,aes(x = Distance/1000,y = r2,colour = Subgroup)) +
  geom_line(linewidth = 0.4) +  # 增加线条粗细
  scale_color_manual(values = c(XianLandrace2="#FF950E",XianCultivar2="#83CAFF",XianLandrace1="#C5000B",XianCultivar1="#579D1C")) +  # 自定义颜色
  labs(x = "Distance (kb)", y = expression(r^2),color = "") +
  theme_bw() +  #使绘图背景为白色
  theme(panel.grid.major = element_blank(),  #移除主要网格线
        panel.grid.minor = element_blank(),  #移除次要网格线
        #axis.ticks.length = unit(0.2,"cm"),      # 修改刻度线长度
        #axis.ticks = element_line(linewidth = 0.7),     # 修改刻度线粗细
        #panel.border = element_rect(fill = NA,colour = "black",linewidth = 2,linetype = 1), #添加图框线
        axis.text.x = element_text(family = "sans",size = 11,colour = "black",face = "plain",angle = 0),  #x轴标签字体大小及样式，可以用axis.title.x指定x轴标题字体、位置、颜色等
        axis.title.x = element_text(size = 11,margin = margin(t=5)),  #轴标题大小及距离x轴标签距离,在 margin 函数中，t 参数表示上边距（距离 x 轴标题到 x 轴标签的距离），r、b 和 l 分别表示右边距、下边距和左边距
        axis.text.y = element_text(family = "sans",size = 11,colour = "black",face = "plain",angle = 0),  #y轴标签字体大小及样式
        axis.title.y = element_text(size = 11,margin = margin(r = 5)))+
  theme(legend.position = "right",  #可以是 "left", "right", "top", "bottom" 或 c(x, y)
        legend.background = element_rect(fill = NA, color = NA), # 图例背景  color控制边框颜色
        legend.title = element_blank(), # 图例标题,居中对齐
        legend.key.size = unit(0.4, "cm"),  #图例每行之间的间距
        legend.text = element_text(size = 8))

ggsave(filename = "667XianSample.LDdecay.pdf",units = "cm",height = 6,width = 10.5)

