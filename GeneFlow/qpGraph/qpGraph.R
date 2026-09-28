#1.先从 0 个 admixture edge 开始
# 4个目标群体 + Oru，第一轮numadmixt=0，先看看单纯的 bifurcating population tree 能不能解释数据？
# ADMIXTOOLS 2 的 find_graphs() 可以自动搜索拓扑，并固定 outgroup。官方文档明确说明它会搜索与观测 f-statistics 相容的 admixture graphs。

library(admixtools)
#先准备 blocked f2，先建库
extract_f2(pref = "../data/667XianSample.28Wild.PASS.missing0.2.maf0.01.Biallelic.OnlyGT.OnlySNP",outdir = "rice_f2",pops = c("XianLandrace1","XianLandrace2","XianCultivar1","XianCultivar2","OutGroup"),blgsize = 1000000,auto_only = FALSE,overwrite = TRUE)
f2_blocks <- f2_from_precomp("rice_f2") #这一步是将我们几个群体几百万SNP建库并保存到硬盘，后续可以直接读取进行所有find_graphs和qpgraph分析




#然后自动搜索拓扑结构，自动寻找 A、B、C、D、Oru 应该如何连接,建议不要只跑一次。因为图搜索存在局部最优问题
#numadmix = 0 不允许 admixture edge，只搜索纯树状拓扑
#如果我们提出一个模型去试验qpgraph函数计算拟合值，实际上是我提出的这一张历史模型，能不能解释数据？
#而find_graphs则是在允许的模型空间里，哪些历史模型最能解释数据
#我们可以先使用find_graphs看看我们的数据更倾向于什么 topology？然后再使用qpgraph函数正式比较A|[B|(C,D)]，B|[A|(C,D)]，(A,B)|(C,D)，(A,C)|(B,D)  (A,D)|(B,C)
#先让 ADMIXTOOLS 完全不允许基因流，看看四个群体最稳定的纯树 topology 到底是什么。 find_graphs() 本身就是通过随机起始 graph 和持续修改 topology 自动寻找与 observed f-statistics 最相容的图，因此跑多次是官方推荐的处理局部最优的方法
library(dplyr)
g0_list <- vector("list", 20)
for (i in 1:20) {set.seed(i);g0_list[[i]] <- find_graphs(data = f2_blocks,numadmix = 0,outpop = "OutGroup", stop_gen = 20,numgraphs = 20)}  #numgraphs = 50 是每一代保留/处理的 graph 数量，不是“最终只搜索 50 张图”；stop_gen = 500 是搜索代数上限
names(g0_list) <- paste0("g", 1:20)
#每一轮提取最佳模型
best_results_g0 <- bind_rows(lapply(seq_along(g0_list), function(i) {g0_list[[i]] %>% filter(is.finite(score)) %>% slice_min(score, n = 1, with_ties = FALSE) %>% mutate(run = i)})) %>% arrange(score)
# 查看20次搜索的最佳 score
print(best_results_g0 %>% select(run, score))
# 所有20次搜索中的总体最佳模型
winner_g0 <- best_results_g0 %>% slice_min(score, n = 1, with_ties = FALSE)
cat("\n[RESULT] Best run:", winner_g0$run, "\n")
cat("[RESULT] Best score:", winner_g0$score, "\n")
# 查看最佳模型的边
print(winner_g0$edges[[1]])
write.csv(winner_g0$edges[[1]],"numadmix0.best.edges.csv")
# 绘制最佳模型
library(ggplot2)
p <- plot_graph(winner_g0$edges[[1]],textsize = 3.5)
p
ggsave("numadmix0.best.graph.pdf",plot = p,width = 7,height = 10)
#真正应该重点看三个东西
1. 20次最佳 score 是否收敛到相同或非常接近的值；
2. 最低 score 的模型 topology 是否反复出现；
3. 这些高拟合 topology 中，你关心的 A/B/C/D 分化顺序是否一致。
尤其第 3 点才是最终想回答“是否存在稳定时间序列”的核心。

#重新计算 winner 的 f4 residual
#因为 find_graphs() 默认主要优化 likelihood score，而 worst_residual 只有在计算 fitted f-statistics 后才能得到。ADMIXTOOLS 官方明确说明 qpgraph(return_fstats = TRUE) 会返回最大的 absolute \(f_4\) residual Z-score。
fit_g0 <- qpgraph(data = f2_blocks,graph = winner_g0$edges[[1]],return_fstats = TRUE,numstart = 100)
fit_g0$score
fit_g0$worst_residual
fit_g0$f4 %>% arrange(desc(abs(z)))
write.csv(fit_g0$f4,"numadmix0.best.f4_residual.csv")
fit_g0_summary <- tibble(numadmix = 0,score = fit_g0$score,worst_residual = fit_g0$worst_residual)
write.csv(fit_g0_summary,"numadmix0.best.summary.csv")
#如果score = 198.904，|worst_residual| = 1.82 < 3，
#那么这棵纯树拟合其实相当不错，这也就是说，我们这几个群体不存在显著的基因流，因为我们在numadmix = 0，强制排除基因流的情况下就能得到一个拟合的非常好的拓扑树
#但是如果score = 198.9038，|worst_residual| = 13.12946 >> 3，
#那么这就是一个很好的信号，暗示我们在排除基因流的情况下无法拟合一个很好的拓扑树，即使对于score = 198.904这个score值已经是最低的拓扑树来说，残差值仍然非常大，也就说明我们无法在不考虑基因流的情况下拟合现有数据




#2.前面从 0 个 admixture edge 开始，我们发现即使是最低score值的拓扑树，残差值仍然非常大，得到拓扑树仍然是不可信的
#我们这里就需要考虑1个或者2个admixture edge

#还是先固定在1个或者2个admixture edge情况下，然后使用find_graphs函数先自由搜索合适的拓扑结构
library(admixtools)
library(dplyr)
g1_list <- vector("list", 20)
for (i in 1:20) {set.seed(i);g1_list[[i]] <- find_graphs(data = f2_blocks,numadmix = 1,outpop = "OutGroup",stop_gen = 100,stop_gen2 = 20,numgraphs = 50)}
names(g1_list) <- paste0("run", 1:20)
#只看20个 winner，把所有搜索到的 graph 合并,因为真实拓扑树未必会是阈值最低的树
#可能有几棵树不存在统计值上的差异，属于数据波动引起的，但最低的那个不是最佳拓扑树
all_g1 <- bind_rows(lapply(seq_along(g1_list), function(i) {g1_list[[i]] %>% mutate(run = i)})) %>% filter(is.finite(score)) %>% arrange(score)
#先看一共有多少
nrow(all_g1)
#增加 graph hash
all_g1$hash <- vapply(all_g1$edges,function(x) {graph_hash(edges_to_igraph(x))},character(1))
#去掉同一 topology 的重复搜索结果,看看有多少真正不同的 graph
unique_g1 <- all_g1 %>% arrange(score) %>% distinct(hash, .keep_all = TRUE)
nrow(unique_g1)
write.csv(unique_g1,"numadmix1.unique_g1.topology.csv")
unique_g1 %>% select(run, score, hash) %>% print(n = 30)
#重点检查 near-optimal models,可以查看score<2的模型，这些模型其实都可以作为候选模型
near_g1 <- unique_g1 %>% filter(score < 2) %>% arrange(score)
near_g1 %>% select(run, score, hash) %>% print(n = Inf)  #可以看到有多少个模型符合score < 2
 



#然后是numadmix = 2的情况
library(admixtools)
library(dplyr)
g2_list <- vector("list", 20)
for (i in 1:20) {set.seed(i);g2_list[[i]] <- find_graphs(data = f2_blocks,numadmix = 2,outpop = "OutGroup",stop_gen = 200,stop_gen2 = 30,numgraphs = 100)}
names(g2_list) <- paste0("run", 1:20)

#只看20个 winner，把所有搜索到的 graph 合并,因为真实拓扑树未必会是阈值最低的树
#可能有几棵树不存在统计值上的差异，属于数据波动引起的，但最低的那个不是最佳拓扑树
best_results_g2 <- bind_rows(lapply(seq_along(g2_list), function(i) {g2_list[[i]] %>% filter(is.finite(score)) %>% slice_min(score, n = 1, with_ties = FALSE) %>% mutate(run = i)})) %>% arrange(score)
all_g2 <- bind_rows(lapply(seq_along(g2_list), function(i) {g2_list[[i]] %>% mutate(run = i)})) %>% filter(is.finite(score)) %>% arrange(score)

#先看一共有多少
nrow(all_g2)

#增加 graph hash
all_g2$hash <- vapply(all_g2$edges,function(x) {graph_hash(edges_to_igraph(x))},character(1))

#去掉同一 topology 的重复搜索结果,看看有多少真正不同的 graph
unique_g2 <- all_g2 %>% arrange(score) %>% distinct(hash, .keep_all = TRUE)
unique_g2_summary <- unique_g2 %>% mutate(topology_id = row_number()) %>% select(topology_id,generation,score,mutation,hash,lasthash,run)
write.csv(unique_g2_summary,"numadmix2.unique_g2.summary.csv",row.names = FALSE)

#重点检查 near-optimal models,可以查看score<2的模型，这些模型其实都可以作为候选模型
near_g2 <- unique_g2 %>% filter(score < 2) %>% arrange(score)
near_g2 %>% select(run, score, hash) %>% print(n = Inf)  #可以看到有多少个模型符合score < 2
write.csv(near_g2 %>% select(run, score, hash),"numadmix2.near_g2.score2.topology.csv")

#打印查看并绘图
for (i in seq_len(nrow(near_g2))) {print(near_g2$edges[[i]])}
for (i in seq_len(nrow(near_g2))) {p <- plot_graph(near_g2$edges[[i]],textsize = 3.5,title = paste0("Model ", i,"  score=", signif(near_g2$score[i], 4)));ggsave(paste0("numadmix2.nearoptimal.",i,".pdf"),plot = p,width = 7,height = 10)}
write.csv(near_g2$edges[[4]],"numadmix2.near_g2.model4.edges.csv")

#根据图像去找你认为的最佳模型
#此时发现model1的Smodel1≈0.00245，Smodel4≈1.4301
#两者都是 numadmix=2，所以模型复杂度相同。Model 1 在训练数据上确实拟合得更好。
#但是：1.43 本身也已经是非常小的 score。Model 4 相对于 Model 1 的拟合差异是否具有统计学意义？
#ADMIXTOOLS 官方推荐使用 block-bootstrap + out-of-sample score 来比较 graph，因为这样可以减少模型选择和过拟合造成的影响。
graph1 <- near_g2$edges[[1]]
graph4 <- near_g2$edges[[4]]
fits_14 <- qpgraph_resample_multi(f2_blocks,list(Model1 = graph1,Model4 = graph4),nboot = 200,numstart = 20)
compare_14 <- compare_fits(fits_14[[1]]$score_test,fits_14[[2]]$score_test)
print(compare_14) #compare_fits() 官方定义的 p_emp 是 bootstrap 中两个模型 score 差异的经验双侧 P 值；ci_low 和 ci_high 则是 score difference 的 95% bootstrap interval。
compare_14$p_emp
compare_14$ci_low
compare_14$ci_high
#结果表明p_emp=0.995，CI为[-5.598,11.9929]，Smodel4≈1.4301完全在CI范围内
#因此可以说model1和model4在统计学上没有差异，model4完全可以作为候选模型

#Model 4 自己也要重新精确拟合一次
graph4 <- near_g2$edges[[4]]
fit_g2_model4 <- qpgraph(data = f2_blocks,graph = graph4,return_fstats = TRUE,numstart = 100)
fit_g2_model4$score
fit_g2_model4$worst_residual
fit_g2_model4$f4 %>% arrange(desc(abs(z))) #注意这里的p值不能小于0.05，这里的p值的意思是观测到的 f4（est）和模型预测的 f4（fit）之间，是否存在显著差异。p大于0.05的意思是这一个 f4 statistic 的 observed value 与 Model 4 的 predicted value 并没有显著不一致，表明我们的拟合很好
write.csv(fit_g2_model4$f4 %>% arrange(desc(abs(z))),"numadmix2.best.model4.f4_residual.csv")
fit_g2_model4_summary <- tibble(numadmix = 2,score = fit_g2_model4$score,worst_residual = fit_g2_model4$worst_residual)
write.csv(fit_g2_model4_summary,"numadmix2.best.model4.summary.csv")