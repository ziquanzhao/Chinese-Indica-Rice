#!/usr/bin/env Rscript

rm(list = ls())

# ===================== 依赖包 =====================
required_packages <- c("optparse", "ggplot2", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少 R 包: ", paste(missing_packages, collapse = ", "),
    "\n请先安装后重新运行。",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
})

# ===================== 参数定义 =====================
option_list <- list(
  make_option(
    c("--input-fst"), dest = "input_fst", type = "character",
    help = "Fst 结果文件（必需）"
  ),
  make_option(
    c("--star-gene"), dest = "star_gene", type = "character", default = NULL,
    help = paste0(
      "已知基因文件，仅支持 xlsx/xls；必须包含 StarGene sheet，",
      "且至少包含 Name、Chr、Start、End 四列 [default: NULL]"
    )
  ),
  make_option(
    c("--threshold"), dest = "threshold", type = "double", default = 0.05,
    help = "Fst 分数的 top 比例，取值范围(0,1) [default: %default]"
  ),
  make_option(
    c("--output-dir"), dest = "output_dir", type = "character", default = ".",
    help = "输出目录 [default: 当前目录]"
  ),
  make_option(
    c("--output_plot"), dest = "output_plot", type = "character", default = NULL,
    help = paste0(
      "输出Fst图像的文件前缀，自动追加.<threshold>并同时生成 PNG 和 PDF ",
      "[default: <输入文件去掉.windowed.weir.fst后缀>.<threshold>]"
    )
  ),
  make_option(
    c("--output_SignificantInterval"),
    dest = "output_significant_interval", type = "character", default = NULL,
    help = paste0(
      "显著区间输出文件，包含Chr、Start、End、Fst、GFFGene五列；",
      "Fst读取自WEIGHTED_FST列；如果 Excel 中存在 GFFGene sheet，",
      "GFFGene列将给出命中的GeneID，",
      "[default: <输入文件名去除.windowed.weir.fst>.<threshold>.list]，",
      "同时还会给出一个所有命中的GFF注释中的GeneID列表,",
      "[default: <输入文件名去除.windowed.weir.fst>.<threshold>.GFFGeneID.list]"
    )
  )
)

help_examples <- paste(
  "使用示例:",
  "",
  "  1.按默认 top 5% 阈值绘图，不标注明星基因:",
  "     Rscript Fst.Plot.R --input-fst result.windowed.weir.fst",
  "",
  "  2.按 top 1% 阈值绘图并标注明星基因:",
  paste0(
    "     Rscript Fst.Plot.R --input-fst result.windowed.weir.fst ",
    "--star-gene genes.xlsx --threshold 0.01"
  ),
  "",
  "  3.指定输出目录和文件名前缀:",
  paste0(
    "     Rscript Fst.Plot.R --input-fst result.windowed.weir.fst ",
    "--output-dir results --output_plot fst.top5"
  ),
  sep = "\n"
)

parser <- OptionParser(option_list = option_list, epilogue = help_examples)
opt <- parse_args(parser)

# ===================== 通用函数 =====================
is_blank <- function(x) {
  is.null(x) || length(x) == 0 || is.na(x) || !nzchar(trimws(x))
}

resolve_output_path <- function(path, output_dir) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path
  } else {
    file.path(output_dir, path)
  }
}

build_plot_paths <- function(path) {
  image_extensions <- c("png", "pdf", "jpg", "jpeg", "tif", "tiff", "svg")
  extension <- tolower(tools::file_ext(path))
  prefix <- if (extension %in% image_extensions) {
    tools::file_path_sans_ext(path)
  } else {
    path
  }

  c(png = paste0(prefix, ".png"), pdf = paste0(prefix, ".pdf"))
}

chr_to_number <- function(x) {
  suppressWarnings(as.integer(sub(".*?([0-9]+).*", "\\1", as.character(x))))
}

read_gene_sheet <- function(file, sheet, required = TRUE) {
  sheets <- readxl::excel_sheets(file)
  if (!sheet %in% sheets) {
    if (required) {
      stop("文件中缺少必须的 sheet: ", sheet, call. = FALSE)
    }
    return(NULL)
  }

  gene_data <- as.data.frame(
    readxl::read_excel(
      file,
      sheet = sheet,
      # 全部先按文本读取，可避免超长 sheet 中空列的类型猜测警告；
      # Start、End 会在下方显式转换为数值。
      col_types = "text",
      .name_repair = "minimal"
    ),
    stringsAsFactors = FALSE
  )
  required_cols <- c("Name", "Chr", "Start", "End")
  missing_cols <- setdiff(required_cols, colnames(gene_data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "%s sheet 缺少必要列: %s",
        sheet, paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  gene_data <- gene_data[, required_cols, drop = FALSE]
  gene_data$Name <- trimws(as.character(gene_data$Name))
  gene_data$Chr_num <- chr_to_number(gene_data$Chr)
  gene_data$Start <- suppressWarnings(as.numeric(gene_data$Start))
  gene_data$End <- suppressWarnings(as.numeric(gene_data$End))

  # 容忍 Start、End 顺序写反的记录，并统一为左、右边界。
  left <- pmin(gene_data$Start, gene_data$End)
  right <- pmax(gene_data$Start, gene_data$End)
  gene_data$Start <- left
  gene_data$End <- right

  valid <- complete.cases(gene_data[, c("Name", "Chr_num", "Start", "End")]) &
    nzchar(gene_data$Name)
  n_invalid <- sum(!valid)
  if (n_invalid > 0) {
    warning(sheet, " sheet 中有 ", n_invalid, " 行无效，已忽略", call. = FALSE)
  }

  gene_data[valid, , drop = FALSE]
}

# 为每个区间汇总与之重叠的基因名称。使用按染色体扫描的方式，避免构建巨大的笛卡尔积。
annotate_interval_genes <- function(intervals, genes) {
  annotations <- rep("", nrow(intervals))
  if (nrow(intervals) == 0 || is.null(genes) || nrow(genes) == 0) {
    return(annotations)
  }

  shared_chr <- intersect(unique(intervals$Chr_num), unique(genes$Chr_num))
  for (chr in shared_chr) {
    interval_index <- which(intervals$Chr_num == chr)
    interval_index <- interval_index[
      order(intervals$start[interval_index], intervals$stop[interval_index])
    ]

    chr_genes <- genes[genes$Chr_num == chr, , drop = FALSE]
    chr_genes <- chr_genes[order(chr_genes$Start, chr_genes$End), , drop = FALSE]
    next_gene <- 1L
    active <- integer(0)

    for (idx in interval_index) {
      interval_start <- intervals$start[idx]
      interval_end <- intervals$stop[idx]

      while (next_gene <= nrow(chr_genes) &&
             chr_genes$Start[next_gene] <= interval_end) {
        active <- c(active, next_gene)
        next_gene <- next_gene + 1L
      }

      if (length(active) > 0) {
        # 区间已按起点排序，终点早于当前起点的基因以后也不可能再命中。
        active <- active[chr_genes$End[active] >= interval_start]
      }

      if (length(active) > 0) {
        hit <- active[
          chr_genes$Start[active] <= interval_end &
            chr_genes$End[active] >= interval_start
        ]
        if (length(hit) > 0) {
          annotations[idx] <- paste(unique(chr_genes$Name[hit]), collapse = ";")
        }
      }
    }
  }

  annotations
}

# ===================== 参数检查与默认输出 =====================
if (is_blank(opt$input_fst)) {
  stop("必须提供 --input-fst 参数\n使用 --help 查看帮助", call. = FALSE)
}
if (!file.exists(opt$input_fst)) {
  stop("Fst 结果文件不存在: ", opt$input_fst, call. = FALSE)
}
if (!is.finite(opt$threshold) || opt$threshold <= 0 || opt$threshold > 1) {
  stop("--threshold 必须是大于 0 且小于等于 1 的数字", call. = FALSE)
}
if (is_blank(opt$output_dir)) {
  stop("--output-dir 不能为空", call. = FALSE)
}

annotate_star_gene <- !is_blank(opt$star_gene)
if (annotate_star_gene) {
  if (!file.exists(opt$star_gene)) {
    stop("已知基因文件不存在: ", opt$star_gene, call. = FALSE)
  }
  if (!grepl("\\.(xlsx|xls)$", opt$star_gene, ignore.case = TRUE)) {
    stop("--star-gene 仅支持 xlsx/xls 文件", call. = FALSE)
  }
  if (!requireNamespace("readxl", quietly = TRUE) ||
      !requireNamespace("ggrepel", quietly = TRUE)) {
    stop(
      "使用 --star-gene 时需要 R 包 readxl 和 ggrepel，请先安装后重新运行。",
      call. = FALSE
    )
  }
}

if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(opt$output_dir)) {
  stop("无法创建输出目录: ", opt$output_dir, call. = FALSE)
}

input_name <- basename(opt$input_fst)
input_prefix <- sub(
  "\\.windowed\\.weir\\.fst$", "", input_name,
  ignore.case = TRUE
)
# 对其他以 .fst 结尾的文件名保留兼容性。
if (identical(input_prefix, input_name)) {
  input_prefix <- sub("\\.fst$", "", input_name, ignore.case = TRUE)
}
threshold_text <- format(opt$threshold, scientific = FALSE, trim = TRUE)

plot_name <- if (is_blank(opt$output_plot)) input_prefix else opt$output_plot
plot_path <- resolve_output_path(plot_name, opt$output_dir)
# --output_plot 始终作为前缀使用；即使输入了图片扩展名，也先移除扩展名再追加阈值。
plot_extension <- tolower(tools::file_ext(plot_path))
image_extensions <- c("png", "pdf", "jpg", "jpeg", "tif", "tiff", "svg")
plot_prefix <- if (plot_extension %in% image_extensions) {
  tools::file_path_sans_ext(plot_path)
} else {
  plot_path
}
plot_paths <- build_plot_paths(paste0(plot_prefix, ".Fst.", threshold_text))

significant_name <- if (is_blank(opt$output_significant_interval)) {
  sprintf("%s.%s.list", input_prefix, threshold_text)
} else {
  opt$output_significant_interval
}
significant_path <- resolve_output_path(significant_name, opt$output_dir)
gff_gene_id_path <- file.path(
  opt$output_dir,
  sprintf("%s.%s.GFFGeneID.list", input_prefix, threshold_text)
)

output_parents <- unique(dirname(c(plot_paths, significant_path, gff_gene_id_path)))
for (directory in output_parents) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(directory)) {
    stop("无法创建输出目录: ", directory, call. = FALSE)
  }
}

# ===================== 读取与整理 Fst 数据 =====================
cat("读取 Fst 文件:", opt$input_fst, "\n")
fst_raw <- read.table(
  opt$input_fst,
  header = TRUE,
  sep = "",
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

required_fst_cols <- c("CHROM", "BIN_START", "BIN_END", "WEIGHTED_FST")
missing_fst_cols <- setdiff(required_fst_cols, colnames(fst_raw))
if (length(missing_fst_cols) > 0) {
  stop(
    "Fst 文件缺少必要列: ", paste(missing_fst_cols, collapse = ", "),
    call. = FALSE
  )
}

fst_all <- data.frame(
  chrom = as.character(fst_raw$CHROM),
  start = suppressWarnings(as.numeric(fst_raw$BIN_START)),
  stop = suppressWarnings(as.numeric(fst_raw$BIN_END)),
  # Fst 分数只允许来自 WEIGHTED_FST。
  fst = suppressWarnings(as.numeric(fst_raw$WEIGHTED_FST)),
  stringsAsFactors = FALSE
)
fst_all$Chr_num <- chr_to_number(fst_all$chrom)

# 只有区间本身不完整的行才从 top 比例的分母中剔除。
# fst 为 NA/NaN/Inf 的行仍属于有效区间，必须计入总区间数。
valid_interval <- complete.cases(
  fst_all[, c("chrom", "start", "stop", "Chr_num")]
) & nzchar(trimws(fst_all$chrom))
n_invalid_interval <- sum(!valid_interval)
if (n_invalid_interval > 0) {
  warning(
    "Fst 文件中有 ", n_invalid_interval,
    " 行缺少有效的区间信息，已从 top 比例分母中剔除",
    call. = FALSE
  )
}
fst_all <- fst_all[valid_interval, , drop = FALSE]

if (nrow(fst_all) == 0) {
  stop("过滤后没有包含有效区间信息的 Fst 记录", call. = FALSE)
}

# 容忍输入中的区间边界顺序写反。
fst_left <- pmin(fst_all$start, fst_all$stop)
fst_right <- pmax(fst_all$start, fst_all$stop)
fst_all$start <- fst_left
fst_all$stop <- fst_right
fst_all <- fst_all %>%
  arrange(Chr_num, start, stop) %>%
  mutate(RowIndex = row_number())

# 无计算结果的区间保留在分母中，但不能参与分数排序、绘图或显著区间输出。
valid_score <- is.finite(fst_all$fst)
n_no_result <- sum(!valid_score)
fst_data <- fst_all[valid_score, , drop = FALSE]

if (nrow(fst_data) == 0) {
  stop("Fst 文件中没有任何有限的 fst 分数可用于排序和绘图", call. = FALSE)
}

# ===================== top 比例阈值与显著区间 =====================
# top 区间目标数量由全部有效区间计算，包括没有 fst 结果的区间。
top_target_n <- ceiling(nrow(fst_all) * opt$threshold)
top_rank <- min(top_target_n, nrow(fst_data))
fst_threshold <- sort(fst_data$fst, decreasing = TRUE)[top_rank]
significant_data <- fst_data %>%
  filter(fst >= fst_threshold) %>%
  arrange(Chr_num, start, stop)

cat("计入 top 比例分母的区间数量:", nrow(fst_all), "\n")
cat("其中无 Fst 计算结果的区间数量:", n_no_result, "\n")
cat("具有有效 Fst 分数的区间数量:", nrow(fst_data), "\n")
if (top_target_n > nrow(fst_data)) {
  warning(
    "top 比例对应 ", top_target_n, " 个区间，但只有 ", nrow(fst_data),
    " 个区间具有有效 Fst 分数；将全部有效分数区间视为显著区间",
    call. = FALSE
  )
}
cat(
  sprintf(
    paste0(
      "top %g%% 目标区间数量: %d；阈值: %.10g；",
      "达到阈值的有效区间数量: %d\n"
    ),
    opt$threshold * 100, top_target_n,
    fst_threshold, nrow(significant_data)
  )
)

# ===================== 读取并匹配 StarGene / GFFGene =====================
star_gene_data <- NULL
gff_gene_data <- NULL
known_gene_hits <- data.frame()
label_data <- data.frame()

if (annotate_star_gene) {
  cat("读取已知基因文件:", opt$star_gene, "\n")
  star_gene_data <- read_gene_sheet(opt$star_gene, "StarGene", required = TRUE)
  cat("StarGene 有效记录数量:", nrow(star_gene_data), "\n")

  sheets <- readxl::excel_sheets(opt$star_gene)
  if ("GFFGene" %in% sheets) {
    gff_gene_data <- read_gene_sheet(opt$star_gene, "GFFGene", required = FALSE)
    cat("GFFGene 有效记录数量:", nrow(gff_gene_data), "\n")
  } else {
    cat("未找到 GFFGene sheet，显著区间结果中的 GFFGene 列将留空\n")
  }

  hit_list <- vector("list", nrow(star_gene_data))
  if (nrow(significant_data) > 0 && nrow(star_gene_data) > 0) {
    for (i in seq_len(nrow(star_gene_data))) {
      hit_idx <- which(
        significant_data$Chr_num == star_gene_data$Chr_num[i] &
          significant_data$start <= star_gene_data$End[i] &
          significant_data$stop >= star_gene_data$Start[i]
      )

      if (length(hit_idx) > 0) {
        hit_list[[i]] <- data.frame(
          RowIndex = significant_data$RowIndex[hit_idx],
          GeneName = star_gene_data$Name[i],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  hit_list <- Filter(Negate(is.null), hit_list)

  if (length(hit_list) > 0) {
    known_gene_hits <- bind_rows(hit_list) %>% distinct()
    label_data <- known_gene_hits %>%
      left_join(
        fst_data[, c("RowIndex", "fst"), drop = FALSE],
        by = "RowIndex"
      ) %>%
      group_by(GeneName) %>%
      slice_max(order_by = fst, n = 1, with_ties = FALSE) %>%
      ungroup()

    cat(
      "命中明星基因的显著区间数量:",
      length(unique(known_gene_hits$RowIndex)), "\n"
    )
    cat("命中的明星基因数量:", length(unique(known_gene_hits$GeneName)), "\n")
  } else {
    cat("没有显著区间与 StarGene 中的基因重叠\n")
  }
}

# 无论是否提供 GFFGene sheet，结果中始终保留 GFFGene 列。
significant_out <- data.frame(
  Chr = significant_data$chrom,
  # 禁用科学计数法，保证基因组坐标在文本结果中保持直观。
  Start = format(significant_data$start, scientific = FALSE, trim = TRUE),
  End = format(significant_data$stop, scientific = FALSE, trim = TRUE),
  Fst = significant_data$fst,
  GFFGene = annotate_interval_genes(significant_data, gff_gene_data),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.table(
  significant_out,
  file = significant_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  fileEncoding = "UTF-8"
)
cat("显著区间结果已保存:", significant_path, "\n")
if (!is.null(gff_gene_data)) {
  gff_gene_fields <- significant_out$GFFGene[
    !is.na(significant_out$GFFGene) & nzchar(significant_out$GFFGene)
  ]
  gff_gene_ids <- if (length(gff_gene_fields) > 0) {
    trimws(unlist(strsplit(gff_gene_fields, ";", fixed = TRUE)))
  } else {
    character(0)
  }
  gff_gene_ids <- sort(
    unique(gff_gene_ids[nzchar(gff_gene_ids)]),
    method = "radix"
  )

  writeLines(gff_gene_ids, con = gff_gene_id_path, useBytes = TRUE)

  cat("GFFGene 注释命中区间数量:", sum(nzchar(significant_out$GFFGene)), "\n")
  cat("GFFGene 注释命中唯一基因数量:", length(gff_gene_ids), "\n")
  cat("GFFGene ID 列表已保存:", gff_gene_id_path, "\n")
}

# ===================== 构建全基因组连续坐标 =====================
chr_info <- fst_all %>%
  group_by(Chr_num) %>%
  summarise(
    chr_min = min(start, na.rm = TRUE),
    chr_max = max(stop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Chr_num) %>%
  mutate(
    chr_len = chr_max,
    offset = lag(cumsum(chr_len), default = 0)
  )

plot_data <- fst_data %>%
  left_join(chr_info[, c("Chr_num", "offset")], by = "Chr_num") %>%
  mutate(Pos = (start + stop) / 2 + offset)

axis_data <- chr_info %>%
  transmute(
    Chr_num = Chr_num,
    center = offset + (chr_min + chr_max) / 2,
    label = sprintf("Chr%02d", Chr_num)
  )

unique_chr <- sort(unique(plot_data$Chr_num))
chr_colors <- rep(c("#579D1C", "#FF950E"), length.out = length(unique_chr))
names(chr_colors) <- as.character(unique_chr)
plot_data$KnownGeneHit <- plot_data$RowIndex %in% unique(known_gene_hits$RowIndex)

if (nrow(label_data) > 0) {
  label_data <- label_data %>%
    select(RowIndex, GeneName) %>%
    left_join(plot_data[, c("RowIndex", "Pos", "fst")], by = "RowIndex") %>%
    distinct(GeneName, .keep_all = TRUE)
}

# ===================== 绘图 =====================
cat("绘制 Fst 全基因组图:\n")
cat(" PNG:", plot_paths[["png"]], "\n")
cat(" PDF:", plot_paths[["pdf"]], "\n")

y_min <- min(c(0, plot_data$fst), na.rm = TRUE)
y_max <- max(c(plot_data$fst, fst_threshold), na.rm = TRUE)
y_range <- y_max - y_min
y_padding <- max(0.02, y_range * 0.08)
plot_font_family <- "sans"

p <- ggplot(plot_data, aes(x = Pos, y = fst)) +
  geom_point(aes(color = factor(Chr_num)), size = 0.35,shape=16) +
  scale_color_manual(values = chr_colors) +
  geom_hline(yintercept = fst_threshold,color = "black", linetype = 2, linewidth = 0.25) +
  scale_x_continuous(labels = axis_data$label,breaks = axis_data$center,expand = expansion(mult = c(0.005, 0.005))) +
  scale_y_continuous(
    limits = c(y_min, y_max + y_padding),
    expand = expansion(mult = c(0.005, 0))
  ) +
  labs(x = "", y = "Fst") +
  theme_bw(base_family = plot_font_family) +
  theme(
    text = element_text(family = plot_font_family),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.ticks = element_line(linewidth = 0.25),
	axis.ticks.length = ggplot2::unit(0.05,"cm"),
    axis.text.x = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(size = 7, colour = "black"),
    axis.title.y = element_text(size = 7, margin = margin(r = 2.5)),
    legend.position = "none"
  )

# 与明星基因重叠的显著窗口用红色覆盖；每个基因只标注 Fst 最高的窗口。
if (any(plot_data$KnownGeneHit)) {
  p <- p + geom_point(
    data = plot_data[plot_data$KnownGeneHit, , drop = FALSE],
    aes(x = Pos, y = fst),
    inherit.aes = FALSE,
    color = "red",
    size = 0.7,shape=16
  )
}

if (nrow(label_data) > 0) {
  p <- p + ggrepel::geom_text_repel(
    data = label_data,
    aes(x = Pos, y = fst, label = GeneName),
    inherit.aes = FALSE,
    color = "black",
    family = plot_font_family,
    fontface = "italic",
    size = 2.2,
    box.padding = 0.3,
    point.padding = 0.15,
    force = 1,
    force_pull = 0.5,
    segment.color = "grey40",
    segment.size = 0.2,
    max.overlaps = 100,
    min.segment.length = 0,
    max.iter = 5000,
    max.time = 5
  )
}

ggsave(
  filename = plot_paths[["png"]],
  plot = p,
  width = 4, height = 1.2, dpi = 600
)
ggsave(
  filename = plot_paths[["pdf"]],
  plot = p,
  width = 4, height = 1.2,
  device = grDevices::cairo_pdf
)

cat("完成。\n")
