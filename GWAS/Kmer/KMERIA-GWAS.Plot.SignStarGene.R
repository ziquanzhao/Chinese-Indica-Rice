#!/usr/bin/env Rscript

rm(list = ls())

# ===================== 加载包 =====================
required_packages <- c("optparse", "ggplot2", "dplyr", "readxl", "ggrepel")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# ===================== 参数定义 =====================
option_list <- list(
  make_option(c("--input"), type = "character",
              help = "GWAS结果文件，必须包含 KmerID, Chr, Start, P-value 四列 (必需)"),

  make_option(c("--star-gene"), type = "character", default = NULL,
              help = "已知基因文件，支持 xlsx/xls/txt/tsv；至少包含 Name, Chr, Start, End 四列 [default: NULL]"),

  make_option(c("--ld-annovation-gene"), dest = "ld_annovation_gene",
              type = "integer", default = NULL,
              help = "基因注释距离(bp)；提供 --star-gene 时必需 [default: NULL]"),

  make_option(c("--ld-remove-alone-site"), dest = "ld_remove_alone_site",
              type = "integer", default = 100,
              help = "孤点判断和连锁组划分距离(bp) [default: %default]"),

  make_option(c("--ld-group-site-number"), dest = "ld_group_site_number",
              type = "integer", default = 10,
              help = "连锁组最少位点数；少于该阈值的连锁组将被删除 [default: %default]"),

  make_option(c("--threshold"), type = "double", default = 5,
              help = "显著性阈值(-log10P) [default: %default]"),
			  
  make_option(c("--y-start"), dest = "y_start", type = "double", default = 0,
              help = "y轴起点(-log10P)，须为大于或等于 0 的有限数值 [default: %default]"),
			  
  make_option(c("--phenotype"), type = "character",
              help = "表型名称 [default: 从 --input 文件名中移除 .Plot.list 后缀]"),

  make_option(c("--output-dir"), type = "character", default = ".",
            help = "输出目录 [default: 当前目录]"),

  make_option(c("--output_plot"), type = "character",
              help = "曼哈顿图输出文件名或前缀，将同时生成同名 PNG 和 PDF [default: <phenotype>.GWAS.manhattan.<threshold>]"),

  make_option(c("--output_SignificantSite"), type = "character",
              help = "显著位点输出文件名 [default: <phenotype>.GWAS.SignificantSite.<threshold>.list]"),

  make_option(c("--output_remove_alone_site"), type = "character",
              help = "过滤孤点和小连锁组后的数据文件名 [default: <phenotype>.GWAS.RemoveAloneSite.list]"),

  make_option(c("--gff-pep"), type = "character", default = NULL,
              help = "基因组总蛋白序列文件 (FASTA格式)；提供后将从中提取 GFFGene 命中的蛋白序列 [default: NULL]")
)

help_examples <- paste(
  "使用示例:",
  "",
  "  # 标注显著位点附近 100 kb 范围内的已知基因:",
  "     Rscript KMERIA-GWAS.Plot.SignStarGene.R --input HeadingDate.k31.gemma.blastn.Final.Plot.list --star-gene StarGene_HeadingDate_candidates.xlsx --ld-annovation-gene 100000",
  sep = "\n"
)

parser <- OptionParser(
  option_list = option_list,
  epilogue = help_examples
)

opt <- parse_args(parser)

# ===================== 参数检查 =====================
if (is.null(opt$input) || !nzchar(trimws(opt$input))) {
  stop("必须提供 --input 参数\n使用 --help 查看帮助")
}
if (!is.finite(opt$y_start) || opt$y_start < 0) {
  stop("--y-start 必须是大于或等于 0 的有限数值")
}

# 未提供 --phenotype 时，从输入文件名中移除 .plot.list 后缀作为表型名称
phenotype_inferred <- FALSE
if (is.null(opt$phenotype) || !nzchar(trimws(opt$phenotype))) {
  input_filename <- basename(opt$input)
  opt$phenotype <- sub("\\.plot\\.list$", "", input_filename, ignore.case = TRUE)

  if (!nzchar(opt$phenotype)) {
    stop("无法从 --input 文件名推断表型名称，请显式提供 --phenotype")
  }

  phenotype_inferred <- TRUE
}

# --ld-remove-alone-site 是孤点判断的必需参数
if (is.null(opt$ld_remove_alone_site)) {
  stop("必须提供 --ld-remove-alone-site 参数，用于判断孤点\n使用 --help 查看帮助")
}
if (is.na(opt$ld_remove_alone_site) || opt$ld_remove_alone_site < 0) {
  stop("--ld-remove-alone-site 必须是大于或等于 0 的整数")
}
if (is.na(opt$ld_group_site_number) || opt$ld_group_site_number < 1) {
  stop("--ld-group-site-number 必须是大于或等于 1 的整数")
}

annotate_known_gene <- !is.null(opt$`star-gene`)

if (!is.null(opt$ld_annovation_gene) &&
    (is.na(opt$ld_annovation_gene) || opt$ld_annovation_gene < 0)) {
  stop("--ld-annovation-gene 必须是大于或等于 0 的整数")
}
if (annotate_known_gene && is.null(opt$ld_annovation_gene)) {
  stop("提供 --star-gene 时必须同时提供 --ld-annovation-gene")
}

# ===================== 输出目录检查 =====================
if (!dir.exists(opt$`output-dir`)) {
  dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)
}

# 将运行日志同时输出到终端和 --output-dir 下的固定日志文件
log_file <- file.path(
  opt$`output-dir`,
  "KMERIA-GWAS.Plot.SignStarGene.log"
)
sink(log_file, append = FALSE, split = TRUE)

cat("日志文件:", log_file, "\n")
if (phenotype_inferred) {
  cat("未提供 --phenotype，自动使用表型名称:", opt$phenotype, "\n")
}

# ===================== 默认输出 =====================
if (is.null(opt$output_plot)) {
  opt$output_plot <- file.path(
    opt$`output-dir`,
    sprintf("%s.GWAS.manhattan.%g.png", opt$phenotype, opt$threshold)
  )
} else {
  opt$output_plot <- file.path(opt$`output-dir`, opt$output_plot)
}

if (is.null(opt$output_SignificantSite)) {
  opt$output_SignificantSite <- file.path(
    opt$`output-dir`,
    sprintf("%s.GWAS.SignificantSite.%g.list", opt$phenotype, opt$threshold)
  )
} else {
  opt$output_SignificantSite <- file.path(opt$`output-dir`, opt$output_SignificantSite)
}

if (is.null(opt$output_remove_alone_site)) {
  opt$output_remove_alone_site <- file.path(
    opt$`output-dir`,
    sprintf("%s.GWAS.RemoveAloneSite.list", opt$phenotype)
  )
} else {
  opt$output_remove_alone_site <- file.path(
    opt$`output-dir`,
    opt$output_remove_alone_site
  )
}

# 根据输出文件名或前缀生成同名 PNG 和 PDF 路径
build_plot_paths <- function(path) {
  image_extensions <- c("png", "pdf", "jpg", "jpeg", "tif", "tiff", "svg")
  extension <- tolower(tools::file_ext(path))
  prefix <- if (extension %in% image_extensions) {
    tools::file_path_sans_ext(path)
  } else {
    path
  }

  c(
    png = paste0(prefix, ".png"),
    pdf = paste0(prefix, ".pdf")
  )
}

manhattan_paths <- build_plot_paths(opt$output_plot)

# GFFGene ID 汇总文件和蛋白提取文件的默认路径（运行时按需生成）
opt$output_GFFGene_list <- file.path(
  opt$`output-dir`,
  sprintf("%s.GWAS.SignificantSite.GFFGene.%g.list", opt$phenotype, opt$threshold)
)
opt$output_GFFGene_pep <- file.path(
  opt$`output-dir`,
  sprintf("%s.GWAS.SignificantSite.GFFGene.%g.pep.fa", opt$phenotype, opt$threshold)
)

# ===================== 读取 GWAS 数据
# KMERIA 整理结果格式（制表符分隔，含表头）：
# KmerID  Chr  Start  P-value
# =====================
cat("\n读取文件:", opt$input, "\n")

gwas_raw <- read.table(opt$input,
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE,
                       check.names = FALSE,
                       quote = "",
                       comment.char = "")

required_gwas_cols <- c("KmerID", "Chr", "Start", "P-value")
miss_gwas_cols <- setdiff(required_gwas_cols, colnames(gwas_raw))
if (length(miss_gwas_cols) > 0) {
  stop(sprintf(
    "KMERIA GWAS结果文件缺少必要列: %s；应包含列: KmerID, Chr, Start, P-value",
    paste(miss_gwas_cols, collapse = ", ")
  ))
}

# ===================== 数据整理 =====================
my_data <- data.frame(
  Chr_label = gwas_raw$Chr,
  CHR       = suppressWarnings(as.integer(gsub("[^0-9]", "", gwas_raw$Chr))),
  BP        = suppressWarnings(as.numeric(gwas_raw$Start)),
  P         = suppressWarnings(as.numeric(gwas_raw[["P-value"]])),
  KmerID    = gwas_raw$KmerID,
  stringsAsFactors = FALSE
)

# 内部位点标识由标准化后的染色体和物理位置自动生成，不再依赖输入 ID 列
my_data$SNP <- paste0(
  sprintf("Chr%02d", my_data$CHR),
  "_",
  format(my_data$BP, scientific = FALSE, trim = TRUE)
)

# 去除 NA
my_data <- my_data[complete.cases(my_data), ]

# 去除非法 P 值
my_data <- my_data[!is.na(my_data$P) & my_data$P > 0 & my_data$P < 1, ]

if (nrow(my_data) == 0) {
  stop("过滤后没有可用于分析的有效位点，请检查输入文件中的P值")
}

# ===================== 按 Chr 和 Start 去重 =====================
# 同一物理位置仅保留 P-value 最小的记录；最小值并列时保留输入顺序最前的记录。
n_before_deduplication <- nrow(my_data)
duplicate_position_summary <- my_data %>%
  count(CHR, BP, name = "OccurrenceNumber")
n_duplicate_positions <- sum(duplicate_position_summary$OccurrenceNumber > 1)

my_data <- my_data %>%
  mutate(InputOrder = row_number()) %>%
  group_by(CHR, BP) %>%
  slice_min(order_by = P, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(InputOrder) %>%
  select(-InputOrder)

n_duplicate_rows_removed <- n_before_deduplication - nrow(my_data)

cat("\n按 Chr 和 Start 去重:\n")
cat("去重前有效行数:", n_before_deduplication, "\n")
cat("存在重复的物理位置数量:", n_duplicate_positions, "\n")
cat("删除的重复行数量:", n_duplicate_rows_removed, "\n")
cat("去重后有效行数:", nrow(my_data), "\n")

# 排序
my_data <- my_data[order(my_data$CHR, my_data$BP), ]

# 增加 -log10P，并先按 --threshold 区分显著点和背景点
my_data$logP <- -log10(my_data$P)
p_threshold <- 10^(-opt$threshold)
threshold_text <- sprintf("%g", opt$threshold)
plot_title <- sub(
  "(?<=\\.gemma\\.)[^.]+(?=\\.Final$)",
  threshold_text,
  opt$phenotype,
  ignore.case = TRUE,
  perl = TRUE
)

significant_data <- my_data[my_data$P <= p_threshold, ]
background_data <- my_data[my_data$P > p_threshold, ]

cat("\nGWAS显著性阈值 (-log10(P)):", opt$threshold, "\n")
cat("对应P-value阈值:", format(p_threshold, scientific = TRUE), "\n")
cat("阈值筛选前位点数量:", nrow(my_data), "\n")
cat("阈值线以上显著位点数量:", nrow(significant_data), "\n")
cat("阈值线以下背景位点数量:", nrow(background_data), "\n")

# ===================== 仅过滤显著位点中的孤点和小连锁组 =====================
# 同一染色体内，相邻位点间距小于或等于 --ld-remove-alone-site 时视为连续连接，
# 由此形成的每个连通区间定义为一个连锁组；背景点不参与距离和数量计算。
significant_data <- significant_data %>%
  group_by(CHR) %>%
  mutate(
    UpstreamDistance = BP - lag(BP),
    DownstreamDistance = lead(BP) - BP,
    NewLDGroup = is.na(UpstreamDistance) |
      UpstreamDistance > opt$ld_remove_alone_site,
    LDGroup = cumsum(NewLDGroup)
  ) %>%
  group_by(CHR, LDGroup) %>%
  mutate(
    LDGroupSiteNumber = n(),
    IsAlone = LDGroupSiteNumber == 1,
    RemoveLDGroup = !IsAlone &
      LDGroupSiteNumber < opt$ld_group_site_number
  ) %>%
  ungroup()

n_significant_before_linkage_filter <- nrow(significant_data)
n_alone_site <- sum(significant_data$IsAlone)
ld_group_summary <- significant_data %>%
  filter(!IsAlone) %>%
  distinct(CHR, LDGroup, LDGroupSiteNumber)
n_ld_group <- nrow(ld_group_summary)
n_small_ld_group <- sum(
  ld_group_summary$LDGroupSiteNumber < opt$ld_group_site_number
)
n_small_ld_group_site <- sum(significant_data$RemoveLDGroup)

significant_data_filtered <- significant_data[
  !significant_data$IsAlone & !significant_data$RemoveLDGroup,
]

# 曼哈顿图保留全部背景点，仅删除不满足连锁条件的显著点
plot_input_data <- bind_rows(background_data, significant_data_filtered) %>%
  arrange(CHR, BP)

remove_alone_out <- plot_input_data[, c("KmerID", "Chr_label", "BP", "P")]
colnames(remove_alone_out) <- c("KmerID", "Chr", "Start", "P-value")

write.table(remove_alone_out,
            file = opt$output_remove_alone_site,
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE,
            fileEncoding = "UTF-8")

cat("\n孤点判断和连锁组划分距离 (--ld-remove-alone-site):",
    opt$ld_remove_alone_site, "bp\n")
cat("参与连锁判定的显著位点数量:", n_significant_before_linkage_filter, "\n")
cat("移除的显著孤点数量:", n_alone_site, "\n")
cat("非孤点连锁组数量:", n_ld_group, "\n")
cat("连锁组最少位点数阈值:", opt$ld_group_site_number, "\n")
cat("删除的小连锁组数量:", n_small_ld_group, "\n")
cat("小连锁组中删除的显著位点数量:", n_small_ld_group_site, "\n")
cat("最终保留显著位点数量:", nrow(significant_data_filtered), "\n")
cat("曼哈顿图保留背景位点数量:", nrow(background_data), "\n")
cat("曼哈顿图最终绘制位点数量:", nrow(plot_input_data), "\n")
cat("过滤结果已保存:", opt$output_remove_alone_site, "\n")

if (nrow(plot_input_data) == 0) {
  stop("阈值及连锁过滤后没有剩余绘图位点，请检查 --threshold、--ld-remove-alone-site 和 --ld-group-site-number 设置")
}

my_data <- plot_input_data[, c("SNP", "CHR", "BP", "P", "logP")]

unique_chr <- sort(unique(my_data$CHR))
chr_labels <- sprintf("Chr%02d", unique_chr)

# ===================== 提取通过连锁过滤的显著位点 =====================
sig <- my_data[my_data$P <= p_threshold, c("SNP", "CHR", "BP", "P", "logP")]
sig$Chr <- sprintf("Chr%02d", sig$CHR)
sig$POS <- sig$BP

if (annotate_known_gene) {
  sig$LD_Left  <- sig$BP - opt$ld_annovation_gene
  sig$LD_Right <- sig$BP + opt$ld_annovation_gene
}

sig_out <- sig[, c("Chr", "POS", "P")]
colnames(sig_out) <- c("Chr", "POS", "Pvalue")

# StarGene 列将在基因匹配完成后填充，写文件操作也在基因匹配之后执行

# ===================== 如果提供 --star-gene，则读取并匹配
# 规则：显著位点 ±--ld-annovation-gene 与基因区间重叠
# =====================
known_gene_hits <- data.frame()
label_data <- data.frame()

if (annotate_known_gene) {

  gene_file <- opt$`star-gene`
  cat("\n读取已知基因文件:", gene_file, "\n")
  cat("基因注释距离 (--ld-annovation-gene):",
      opt$ld_annovation_gene, "bp\n")

  if (grepl("\\.xlsx$|\\.xls$", gene_file, ignore.case = TRUE)) {
    gene_data <- readxl::read_excel(gene_file, sheet = "StarGene")
    gene_data <- as.data.frame(gene_data, stringsAsFactors = FALSE)
  } else {
    gene_data <- read.table(gene_file,
                            header = TRUE,
                            sep = "\t",
                            stringsAsFactors = FALSE,
                            check.names = FALSE,
                            quote = "",
                            comment.char = "")
  }

  required_cols <- c("Name", "Chr", "Start", "End")
  miss_cols <- setdiff(required_cols, colnames(gene_data))
  if (length(miss_cols) > 0) {
    stop(sprintf("已知基因文件缺少必要列: %s",
                 paste(miss_cols, collapse = ", ")))
  }

  gene_data <- gene_data[, required_cols]
  gene_data$Start <- suppressWarnings(as.numeric(gene_data$Start))
  gene_data$End  <- suppressWarnings(as.numeric(gene_data$End))

  gene_data$Chr_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", gene_data$Chr)))
  gene_data <- gene_data[complete.cases(gene_data[, c("Name", "Chr_num", "Start", "End")]), ]

  overlap_list <- list()

  if (nrow(sig) > 0 && nrow(gene_data) > 0) {
    for (i in seq_len(nrow(sig))) {
      hit <- gene_data[
        gene_data$Chr_num == sig$CHR[i] &
          sig$LD_Left[i] <= gene_data$End &
          sig$LD_Right[i] >= gene_data$Start,
      ]

      if (nrow(hit) > 0) {
        tmp <- data.frame(
          SNP = sig$SNP[i],
          Sig_Chr = sprintf("Chr%02d", sig$CHR[i]),
          Sig_POS = sig$BP[i],
          Pvalue = sig$P[i],
          logP = sig$logP[i],
          GeneName = hit$Name,
          GeneChr = sprintf("Chr%02d", hit$Chr_num),
          GeneStart = hit$Start,
          GeneEnd = hit$End,
          stringsAsFactors = FALSE
        )
        overlap_list[[length(overlap_list) + 1]] <- tmp
      }
    }
  }

  if (length(overlap_list) > 0) {
    known_gene_hits <- dplyr::bind_rows(overlap_list)

    cat("关联到已知基因的显著位点数量:", length(unique(known_gene_hits$SNP)), "\n")
    cat("命中的已知基因数量:", length(unique(known_gene_hits$GeneName)), "\n")

    # 每个基因只保留最显著的那个位点用于图上标注
    label_data <- known_gene_hits %>%
      group_by(GeneName) %>%
      slice_min(order_by = Pvalue, n = 1, with_ties = FALSE) %>%
      ungroup()

  } else {
    cat("没有显著位点在 ±", opt$ld_annovation_gene,
        " bp 范围内命中已知基因\n", sep = "")
  }
}

# ===================== 填充 StarGene 列 =====================
if (annotate_known_gene && nrow(known_gene_hits) > 0) {
  # 按位点 SNP 聚合所有命中的基因名（;分隔）
  snp_to_gene <- known_gene_hits %>%
    group_by(SNP) %>%
    summarise(StarGene = paste(unique(GeneName), collapse = ";"), .groups = "drop")

  sig_snp <- paste0(sig$Chr, "_", sig$POS)
  merged <- merge(
    data.frame(SNP = sig_snp, Chr = sig$Chr, POS = sig$POS, Pvalue = sig$P,
               stringsAsFactors = FALSE),
    snp_to_gene,
    by = "SNP",
    all.x = TRUE
  )
  merged$StarGene[is.na(merged$StarGene)] <- ""
  sig_out <- merged[, c("Chr", "POS", "Pvalue", "StarGene")]
} else if (annotate_known_gene) {
  # 提供了 --star-gene 但无命中，StarGene 列全为空
  sig_out$StarGene <- ""
}

# ===================== 填充 GFFGene 列
# 规则：仅对 StarGene 列为空的显著位点，用与 StarGene 相同的
#       ±--ld-annovation-gene 窗口
#       与基因区间重叠进行注释
# GFFGene 表读自 --star-gene 文件的第二个工作表 "GFFGene"
# =====================
if (annotate_known_gene) {
  gene_file <- opt$`star-gene`

  gff_data <- tryCatch({
    if (grepl("\\.xlsx$|\\.xls$", gene_file, ignore.case = TRUE)) {
      sheets <- readxl::excel_sheets(gene_file)
      if ("GFFGene" %in% sheets) {
        df <- readxl::read_excel(gene_file, sheet = "GFFGene")
        as.data.frame(df, stringsAsFactors = FALSE)
      } else {
        cat("警告: --star-gene 文件中未找到工作表 'GFFGene'，跳过 GFFGene 注释\n")
        NULL
      }
    } else {
      cat("警告: --star-gene 不是 xlsx/xls 文件，跳过 GFFGene 注释\n")
      NULL
    }
  }, error = function(e) {
    cat("警告: 读取 GFFGene 工作表失败:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(gff_data)) {
    required_gff_cols <- c("Name", "Chr", "Start", "End")
    miss_gff <- setdiff(required_gff_cols, colnames(gff_data))
    if (length(miss_gff) > 0) {
      cat(sprintf("警告: GFFGene 工作表缺少列 %s，跳过 GFFGene 注释\n",
                  paste(miss_gff, collapse = ", ")))
      gff_data <- NULL
    }
  }

  if (!is.null(gff_data)) {
    gff_data$Start   <- suppressWarnings(as.numeric(gff_data$Start))
    gff_data$End     <- suppressWarnings(as.numeric(gff_data$End))
    gff_data$Chr_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", gff_data$Chr)))
    gff_data <- gff_data[complete.cases(gff_data[, c("Name", "Chr_num", "Start", "End")]), ]

    cat("读取 GFFGene 工作表，共", nrow(gff_data), "条基因记录\n")

    # 确保 sig_out 中有 StarGene 列（无论是否有命中）
    if (!"StarGene" %in% colnames(sig_out)) sig_out$StarGene <- ""

    # 对 StarGene 为空的位点做 GFFGene 注释
    sig_out$GFFGene <- ""

    # 重建 sig 的 CHR_num 和 BP 用于匹配
    sig_chr_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", sig_out$Chr)))
    sig_bp      <- suppressWarnings(as.numeric(sig_out$POS))

    for (i in seq_len(nrow(sig_out))) {
      # 只注释 StarGene 为空的行
      if (sig_out$StarGene[i] != "") next

      ld_left  <- sig_bp[i] - opt$ld_annovation_gene
      ld_right <- sig_bp[i] + opt$ld_annovation_gene
      hit_gff <- gff_data[
        gff_data$Chr_num == sig_chr_num[i] &
          ld_left  <= gff_data$End &
          ld_right >= gff_data$Start,
      ]

      if (nrow(hit_gff) > 0) {
        sig_out$GFFGene[i] <- paste(unique(hit_gff$Name), collapse = ";")
      }
    }

    n_gff_hit <- sum(sig_out$GFFGene != "")
    cat("GFFGene 注释命中位点数量:", n_gff_hit, "\n")

  } else {
    # 读取失败或无工作表，GFFGene 列留空
    if (!"StarGene" %in% colnames(sig_out)) sig_out$StarGene <- ""
    sig_out$GFFGene <- ""
  }
}

cat("显著位点数量:", nrow(sig_out), "\n")

write.table(sig_out,
            file = opt$output_SignificantSite,
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE,
            fileEncoding = "UTF-8")

# ===================== 汇总 GFFGene ID 并提取蛋白序列 =====================
if (annotate_known_gene && "GFFGene" %in% colnames(sig_out)) {

  # 拆分所有 GFFGene 字段，去重汇总
  gff_ids_raw <- sig_out$GFFGene[sig_out$GFFGene != ""]
  gff_ids_all <- unique(unlist(strsplit(gff_ids_raw, ";", fixed = TRUE)))
  gff_ids_all <- sort(gff_ids_all[nchar(trimws(gff_ids_all)) > 0])

  if (length(gff_ids_all) > 0) {
    cat("\nGFFGene 去重后基因ID数量:", length(gff_ids_all), "\n")

    # 写出 GFFGene ID 列表文件
    writeLines(gff_ids_all, con = opt$output_GFFGene_list)
    cat("GFFGene ID 列表已保存:", opt$output_GFFGene_list, "\n")

    # 如果提供了 --gff-pep，从总蛋白文件中提取对应序列
    if (!is.null(opt$`gff-pep`)) {
      pep_file <- opt$`gff-pep`
      if (!file.exists(pep_file)) {
        cat("警告: --gff-pep 指定的文件不存在:", pep_file, "\n")
      } else {
        cat("读取蛋白序列文件:", pep_file, "\n")

        # 逐行读取 FASTA，提取命中的序列
        lines <- readLines(pep_file, warn = FALSE)
        out_lines <- character(0)
        capture <- FALSE

        for (ln in lines) {
          if (startsWith(ln, ">")) {
            # 取 > 后第一个空白之前的 ID
            header_id <- sub("^>([^[:space:]]+).*", "\\1", ln)
            capture <- header_id %in% gff_ids_all
          }
          if (capture) out_lines <- c(out_lines, ln)
        }

        n_extracted <- sum(startsWith(out_lines, ">"))
        if (n_extracted > 0) {
          writeLines(out_lines, con = opt$output_GFFGene_pep)
          cat("提取到蛋白序列数量:", n_extracted, "\n")
          cat("蛋白序列已保存:", opt$output_GFFGene_pep, "\n")
          # 提示未找到的ID
          found_ids <- sub("^>([^[:space:]]+).*", "\\1",
                           out_lines[startsWith(out_lines, ">")])
          missing_ids <- setdiff(gff_ids_all, found_ids)
          if (length(missing_ids) > 0) {
            cat("警告: 以下", length(missing_ids), "个 ID 在蛋白文件中未找到:\n")
            cat(paste(" ", missing_ids, collapse = "\n"), "\n")
          }
        } else {
          cat("警告: 蛋白文件中未找到任何 GFFGene ID 对应的序列\n")
        }
      }
    }
  } else {
    cat("\n没有 GFFGene 注释命中，跳过 ID 汇总和蛋白提取\n")
  }
}

# ===================== 生成曼哈顿图坐标
# =====================
chr_info <- my_data %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP), .groups = "drop") %>%
  arrange(CHR)

chr_info$tot <- cumsum(chr_info$chr_len) - chr_info$chr_len

plot_data <- my_data %>%
  left_join(chr_info, by = "CHR") %>%
  mutate(BPcum = BP + tot)

axisdf <- plot_data %>%
  group_by(CHR) %>%
  summarise(center = (max(BPcum) + min(BPcum)) / 2, .groups = "drop")

# 染色体交替颜色
chr_colors <- rep(c("#579D1C", "#FF950E"), length.out = length(unique_chr))
names(chr_colors) <- unique_chr

plot_data$PointColor <- chr_colors[as.character(plot_data$CHR)]

cat("曼哈顿图绘制位点数量:", nrow(plot_data), "\n")

# 明星基因红点数据：保留全部关联明星基因的显著位点，
# 与后续“每个基因只选择一个文字标签”的 label_data 完全独立。
star_gene_point_data <- data.frame()
if (nrow(known_gene_hits) > 0) {
  hit_snps <- unique(known_gene_hits$SNP)
  star_gene_point_data <- plot_data %>%
    filter(SNP %in% hit_snps) %>%
    distinct(SNP, .keep_all = TRUE)

  cat("明星基因关联红色位点数量:", nrow(star_gene_point_data), "\n")
}

# 明星基因文字标签：每个基因仅标注最显著的一个关联位点
if (nrow(label_data) > 0) {
  label_data <- label_data %>%
    left_join(plot_data[, c("SNP", "BPcum")], by = "SNP") %>%
    distinct(GeneName, .keep_all = TRUE)

  cat("明星基因文字标签数量:", nrow(label_data), "\n")
}

# ===================== 绘制曼哈顿图
# =====================
cat("\n绘制曼哈顿图:\n",
    " PNG:", manhattan_paths[["png"]], "\n",
    " PDF:", manhattan_paths[["pdf"]], "\n")

# 起点高于所有位点时，仍保证 y 轴范围有效。
y_upper <- max(max(plot_data$logP, na.rm = TRUE), opt$threshold, opt$y_start) + 1
cat("y轴起点:", opt$y_start, "\n")

p1 <- ggplot(plot_data, aes(x = BPcum, y = logP)) +
  geom_point(aes(color = factor(CHR)), size = 0.25,shape = 16) +
  scale_color_manual(values = chr_colors) +
  geom_hline(yintercept = opt$threshold,color = "black", linetype = 2, linewidth = 0.2) +
  scale_x_continuous(label = sprintf("Chr%02d", axisdf$CHR),breaks = axisdf$center,expand = c(0.005, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(opt$y_start, y_upper)) +
  labs(x = "", y = expression(-log[10](p))) +
  theme_bw(base_family = "sans") +
  theme(
    text = element_text(family = "sans"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.17),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.17),
	axis.ticks.length = unit(0.06,"cm"),
    axis.text.x = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.title.y = element_text(size = 6, margin = margin(r = 2.5)),
    legend.position = "none"
  )

# 所有关联明星基因的显著位点均叠加绘制为红色点
if (nrow(star_gene_point_data) > 0) {
  p1 <- p1 +
    geom_point(data = star_gene_point_data,
      aes(x = BPcum, y = logP),color = "red",size = 0.4,shape = 16)
}

if (nrow(label_data) > 0) {
  p1 <- p1 +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(x = BPcum, y = logP, label = GeneName),
      color = "black",
      family = "sans",
      fontface = "italic",
      size = 1.65,
      box.padding = 0.2,
      point.padding = 0.1,
      force = 1.0,
      force_pull = 0.5,
      segment.color = "grey40",
      segment.size = 0.15,
      max.overlaps = 100,
      min.segment.length = 0,
      max.iter = 5000,         # 增加迭代次数，让排布更充分
      max.time = 5              # 最多花5秒优化布局
    )
}

ggsave(filename = manhattan_paths[["png"]],plot = p1,width = 4, height = 1.1, dpi = 600)

ggsave(filename = manhattan_paths[["pdf"]],plot = p1,width = 4, height = 1.1,device = cairo_pdf)


cat("\n完成:", opt$phenotype, "\n")
sink()
