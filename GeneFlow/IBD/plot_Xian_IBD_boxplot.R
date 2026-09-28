#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

#=========================================================
# 1. 输入文件与输出文件
#=========================================================
stats_files <- c(
  "XianLandrace1.vs.XianLandrace2.stats.csv",
  "XianCultivar1.vs.XianCultivar2.stats.csv",
  "XianLandrace1.vs.XianCultivar1.stats.csv",
  "XianLandrace2.vs.XianCultivar1.stats.csv",
  "XianLandrace1.vs.XianCultivar2.stats.csv",
  "XianLandrace2.vs.XianCultivar2.stats.csv"
)

plot_data_file <- "PlotData_Xian_IBD_total_length_per_sample_pair.tsv"
output_pdf <- "Boxplot_Xian_IBD_total_length_per_sample_pair.pdf"

#=========================================================
# 2. 群体缩写和比较顺序
#=========================================================
# 同时兼容 XianCultivar 和结果文件中的 XianCultivar2，两者均缩写为 XC2。
group_abbr <- c(
  XianLandrace1 = "XL1",
  XianLandrace2 = "XL2",
  XianCultivar1 = "XC1",
  XianCultivar2 = "XC2",
  XianCultivar = "XC2"
)

group_rank <- c(XL1 = 1L, XL2 = 2L, XC1 = 3L, XC2 = 4L)

pair_order <- c(
  "XL1_vs_XL2",
  "XL1_vs_XC1",
  "XL1_vs_XC2",
  "XL2_vs_XC1",
  "XL2_vs_XC2",
  "XC1_vs_XC2"
)

make_pair_key <- function(group1, group2) {
  rank1 <- unname(group_rank[group1])
  rank2 <- unname(group_rank[group2])

  ifelse(
    rank1 <= rank2,
    paste(group1, group2, sep = "_vs_"),
    paste(group2, group1, sep = "_vs_")
  )
}

#=========================================================
# 3. 读取并合并 6 个比较组的统计结果
#=========================================================
missing_files <- stats_files[!file.exists(stats_files)]
if (length(missing_files) > 0L) {
  stop("缺少输入文件：", paste(missing_files, collapse = ", "))
}

required_columns <- c(
  "Sample1", "Sample2", "Group1", "Group2", "Total_IBD_cM"
)

message("正在读取 6 个 IBD 统计文件……")

plot_data <- rbindlist(
  lapply(stats_files, function(file) {
    dat <- fread(file)

    missing_columns <- setdiff(required_columns, names(dat))
    if (length(missing_columns) > 0L) {
      stop(
        "文件 ", file, " 缺少必需列：",
        paste(missing_columns, collapse = ", ")
      )
    }

    group1_abbr <- unname(group_abbr[as.character(dat$Group1)])
    group2_abbr <- unname(group_abbr[as.character(dat$Group2)])

    unknown_groups <- unique(c(
      as.character(dat$Group1[is.na(group1_abbr)]),
      as.character(dat$Group2[is.na(group2_abbr)])
    ))
    unknown_groups <- unknown_groups[!is.na(unknown_groups)]

    if (length(unknown_groups) > 0L) {
      stop(
        "文件 ", file, " 中存在未定义缩写的群体：",
        paste(unknown_groups, collapse = ", ")
      )
    }

    if (!is.numeric(dat$Total_IBD_cM)) {
      converted_value <- suppressWarnings(as.numeric(dat$Total_IBD_cM))
      bad_value <- is.na(converted_value) & !is.na(dat$Total_IBD_cM)

      if (any(bad_value)) {
        stop("文件 ", file, " 的 Total_IBD_cM 列包含非数值内容。")
      }

      dat[, Total_IBD_cM := converted_value]
    }

    data.table(
      SourceFile = file,
      Sample1 = dat$Sample1,
      Sample2 = dat$Sample2,
      Group1 = dat$Group1,
      Group2 = dat$Group2,
      pair_group = make_pair_key(group1_abbr, group2_abbr),
      Total_IBD_cM = dat$Total_IBD_cM
    )
  }),
  use.names = TRUE
)

total_pair_count <- nrow(plot_data)
non_finite_count <- plot_data[, sum(!is.finite(Total_IBD_cM))]
zero_ibd_count <- plot_data[, sum(is.finite(Total_IBD_cM) & Total_IBD_cM == 0)]

# 保留 Total_IBD_cM = 0 的样本对，仅排除缺失值和无穷值。
plot_data <- plot_data[is.finite(Total_IBD_cM)]

message("原始样本对数量：", total_pair_count)
message("保留 Total_IBD_cM = 0 的样本对数量：", zero_ibd_count)
message("排除非有限值的样本对数量：", non_finite_count)

observed_pairs <- unique(plot_data$pair_group)
missing_pairs <- setdiff(pair_order, observed_pairs)
unexpected_pairs <- setdiff(observed_pairs, pair_order)

if (length(missing_pairs) > 0L) {
  stop("输入文件中缺少比较组：", paste(missing_pairs, collapse = ", "))
}

if (length(unexpected_pairs) > 0L) {
  stop("输入文件中出现非预期比较组：", paste(unexpected_pairs, collapse = ", "))
}

plot_data[, pair_group := factor(pair_group, levels = pair_order)]
setorder(plot_data, pair_group, Sample1, Sample2)

fwrite(plot_data, plot_data_file, sep = "\t", quote = FALSE, na = "NA")
message("用于统计和绘图的样本对数量：", nrow(plot_data))
message("已保存绘图数据：", plot_data_file)

#=========================================================
# 4. Tukey HSD 紧凑字母标注
#=========================================================
make_compact_letters <- function(tukey_table, group_means, alpha = 0.05) {
  groups <- names(sort(group_means, decreasing = TRUE))
  significant <- matrix(
    FALSE,
    nrow = length(groups),
    ncol = length(groups),
    dimnames = list(groups, groups)
  )

  for (comparison in rownames(tukey_table)) {
    parts <- strsplit(comparison, "-", fixed = TRUE)[[1]]
    if (length(parts) == 2L) {
      is_significant <- tukey_table[comparison, "p adj"] < alpha
      significant[parts[1], parts[2]] <- is_significant
      significant[parts[2], parts[1]] <- is_significant
    }
  }

  letter_pool <- c(letters, LETTERS, paste0("L", seq_len(100L)))
  letter_groups <- list()
  group_letters <- setNames(vector("list", length(groups)), groups)

  for (group in groups) {
    placed <- FALSE

    for (letter in names(letter_groups)) {
      holders <- letter_groups[[letter]]
      if (all(!significant[group, holders])) {
        letter_groups[[letter]] <- c(holders, group)
        group_letters[[group]] <- c(group_letters[[group]], letter)
        placed <- TRUE
      }
    }

    if (!placed) {
      new_letter <- letter_pool[length(letter_groups) + 1L]
      letter_groups[[new_letter]] <- group
      group_letters[[group]] <- c(group_letters[[group]], new_letter)
    }
  }

  for (i in seq_along(groups)) {
    for (j in seq_along(groups)) {
      if (j <= i) next

      group_i <- groups[i]
      group_j <- groups[j]
      shares_letter <- length(intersect(
        group_letters[[group_i]], group_letters[[group_j]]
      )) > 0L

      if (!significant[group_i, group_j] && !shares_letter) {
        new_letter <- letter_pool[length(letter_groups) + 1L]
        letter_groups[[new_letter]] <- c(group_i, group_j)
        group_letters[[group_i]] <- c(group_letters[[group_i]], new_letter)
        group_letters[[group_j]] <- c(group_letters[[group_j]], new_letter)
      }
    }
  }

  data.table(
    pair_group = names(group_letters),
    Letters = vapply(
      group_letters,
      paste0,
      collapse = "",
      FUN.VALUE = character(1)
    )
  )
}

get_letter_positions <- function(dat, alpha = 0.05) {
  dat_use <- dat[!is.na(Total_IBD_cM)]

  if (uniqueN(dat_use$pair_group) < 2L) {
    return(data.table())
  }

  if (uniqueN(dat_use$Total_IBD_cM) < 2L) {
    letters_dt <- data.table(
      pair_group = pair_order,
      Letters = "a"
    )
  } else {
    aov_result <- aov(Total_IBD_cM ~ pair_group, data = dat_use)
    tukey_result <- TukeyHSD(aov_result, "pair_group")$pair_group
    group_means <- tapply(
      dat_use$Total_IBD_cM,
      dat_use$pair_group,
      mean,
      na.rm = TRUE
    )
    letters_dt <- make_compact_letters(
      tukey_result,
      group_means,
      alpha = alpha
    )
  }

  group_max <- dat_use[, .(
    group_max = max(Total_IBD_cM)
  ), by = pair_group]

  value_range <- diff(range(dat_use$Total_IBD_cM))
  padding <- value_range * 0.06

  if (!is.finite(padding) || padding <= 0) {
    padding <- max(abs(dat_use$Total_IBD_cM)) * 0.06
  }
  if (!is.finite(padding) || padding <= 0) {
    padding <- 1
  }

  letters_dt <- merge(
    letters_dt,
    group_max,
    by = "pair_group",
    all.x = TRUE
  )
  letters_dt[, pair_group := factor(pair_group, levels = pair_order)]
  letters_dt[, label_y := group_max + padding]
  letters_dt[]
}

letter_positions <- tryCatch(
  get_letter_positions(plot_data),
  error = function(e) {
    warning("无法计算 Tukey HSD 字母标注：", conditionMessage(e))
    data.table()
  }
)

#=========================================================
# 5. 绘制箱线图
#=========================================================
p <- ggplot(
  plot_data,
  aes(x = pair_group, y = Total_IBD_cM)
) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.25,
    color = "black",
    linewidth = 0.6
  ) +
  geom_boxplot(
    notch = FALSE,
    width = 0.6,
    color = "black",
    fill = NA,
    outlier.size = 0.6,
    linewidth = 0.6
  ) +
  scale_x_discrete(
    limits = pair_order,
    labels = function(x) gsub("_vs_", ".vs.", x, fixed = TRUE)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    x = NULL,
    y = "Total IBD length per sample pair (cM)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      family = "sans",
      size = 16,
      colour = "black",
      angle = 35,
      hjust = 1
    ),
    axis.text.y = element_text(
      family = "sans",
      size = 18,
      colour = "black"
    ),
    axis.title.y = element_text(size = 18, margin = margin(r = 10)),
    legend.position = "none"
  )

if (nrow(letter_positions) > 0L) {
  p <- p + geom_text(
    data = letter_positions,
    aes(x = pair_group, y = label_y, label = Letters),
    inherit.aes = FALSE,
    vjust = -0.2,
    size = 8
  )
}

ggsave(
  filename = output_pdf,
  plot = p,
  width = 6,
  height = 5,
  units = "in"
)

message("已保存箱线图：", output_pdf)
message("完成。")
