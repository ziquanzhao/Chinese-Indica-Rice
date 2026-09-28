#qpdstat函数只有在f4mode=FALSE 只有在第一个参数是 genotype prefix 时才真正能计算 D

library(admixtools)
library(dplyr)

extract_f2(pref = "../data/667XianSample.28Wild.PASS.missing0.2.maf0.01.Biallelic.OnlyGT.OnlySNP",outdir = "rice_f2",pops = c("XianLandrace1","XianLandrace2","XianCultivar1","XianCultivar2","OutGroup"),blgsize = 1000000,auto_only = FALSE,overwrite = TRUE)
f2_blocks <- f2_from_precomp("rice_f2") #这一步是将我们几个群体几百万SNP建库并保存到硬盘，后续可以直接读取进行所有find_graphs和qpgraph分析

qpdstat(data = f2_blocks,pop1 = "XianCultivar1",pop2 = "OutGroup",pop3 = "XianLandrace2",pop4 = "XianCultivar2")
qpdstat(data = f2_blocks,pop1 = "XianCultivar1",pop2 = "OutGroup",pop3 = "XianLandrace2",pop4 = "XianLandrace1")
qpdstat(data = f2_blocks,pop1 = "XianCultivar2",pop2 = "OutGroup",pop3 = "XianCultivar1",pop4 = "XianLandrace2")
qpdstat(data = f2_blocks,pop1 = "XianCultivar2",pop2 = "OutGroup",pop3 = "XianLandrace2",pop4 = "XianLandrace1")
qpdstat(data = f2_blocks,pop1 = "XianLandrace1",pop2 = "OutGroup",pop3 = "XianLandrace2",pop4 = "XianCultivar1")
qpdstat(data = f2_blocks,pop1 = "XianLandrace1",pop2 = "OutGroup",pop3 = "XianLandrace2",pop4 = "XianCultivar2")
