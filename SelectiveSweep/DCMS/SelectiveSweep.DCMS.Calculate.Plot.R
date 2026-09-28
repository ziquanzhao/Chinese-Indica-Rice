#!/usr/bin/env Rscript

# 基于 FST、核苷酸多样性比值和 XPCLR 计算 MINOTAUR 风格的 DCMS，并绘制曼哈顿图。
# BED 是协方差背景的排除区间，不会删除正式分析区间。

options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "SelectiveSweep.DCMS.Calculate.Plot.R"
active_logger <- NULL

print_help <- function() {
  cat(paste0(
    "用法:\n",
    "  Rscript ", script_name, " --input-dir DIR --group-combination FILE [选项]\n\n",
    "功能:\n",
    "  按 Ref<Ref>.vs.Check<Check> 识别同一亚群对的 PI、FST 和 XPCLR 文件，\n",
    "  仅保留三种方法共有且统计值有限的窗口，将三个统计量转换为右尾经验p 值，并按 MINOTAUR::DCMS 的公式计算 DCMS。\n\n",
    "运行逻辑:\n",
    "  1. 从亚群对列表逐行构建 Ref<Ref>.vs.Check<Check> 文件前缀。\n",
    "  2. 在 --input-dir 顶层查找以该前缀开头、分别以下列后缀结尾的文件：\n",
    "       .log2_pi_ratio.pi      PI\n",
    "       .windowed.weir.fst     FST\n",
    "       .xpclr                 XPCLR\n",
    "  3. 按 CHROM/BIN_START/BIN_END 取三种方法的区间交集并删除非有限值。\n",
    "  4. 如果提供 --screen-bed，BED 中列出的是协方差背景的排除区域（黑名单）；与 BED 有任意交集的窗口仅在计算协方差矩阵和相关矩阵时被排除。这些窗口不会从正式分析删除，仍参与经验 p 值转换、获得 DCMS 并写入结果。\n",
    "  5. 每个亚群对输出 Ref<Ref>.vs.Check<Check>.DCMS.list，并绘制 DCMS 曼哈顿图。\n",
    "  6. 如果提供 --star-gene，在入选窗口中查找与已知基因重叠的区间并进行标注。\n",
    "  7. 带基因标注的 PNG、PDF 和 overlap.tsv 使用 --phenotype-id 作为表型标识。\n\n",
    "必需参数:\n",
    "  --input-dir DIR\n",
    "      输入文件所在目录，仅搜索该目录顶层。\n",
    "  --group-combination FILE\n",
    "      制表符分隔的亚群对列表；前两列必须依次命名为 Ref 和 Check。\n\n",
    "  --screen-bed FILE\n",
    "      协方差背景排除区域文件（黑名单），不是要从正式分析删除的区间列表。与 BED 区域存在任意交集的窗口，只从协方差矩阵和相关矩阵的计算背景中排除；它们仍参与经验 p 值和 DCMS 计算，也仍保留在最终输出中。未提供时，使用全部三种方法共有且统计值有限的窗口计算背景矩阵。\n",
    "  --output-dir DIR\n",
    "      输出目录，默认为当前工作目录。\n",
    "  --top-threshold NUMBER\n",
    "      标记 DCMS 最高比例，范围为 (0, 1]，默认 0.05。并列值全部保留\n",
    "  --star-gene FILE\n",
    "      可选的已知基因文件，支持 xlsx、xls、txt 和 tsv，至少包含 Name、Chr、Start 和 End 四列 [default: NULL]。Excel 文件优先读取名为 StarGene 的工作表；如果不存在，使用第一个工作表。未提供时绘制不带基因标注的图。\n",
    "  --phenotype-id ID\n",
    "      表型 ID，用于带基因标注的 PNG、PDF 和 overlap.tsv 文件名。仅能与 --star-gene 一起使用；未提供时，默认为 --star-gene 的文件名去掉最后一个扩展名。\n",
    "  -h, --help\n",
    "      显示本帮助信息。\n\n",
    "使用示例:\n",
    "  Rscript ", script_name, " --input-dir . --group-combination Combination.list --screen-bed OsNIP-T2T.centromere.telomere.bed --top-threshold 0.05 --star-gene StarGene_PlantHeight_candidates.xlsx --phenotype-id PlantHeight\n"
  ))
}

stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

parse_arguments <- function(args) {
  if (any(args %in% c("-h", "--help"))) {
    print_help()
    quit(save = "no", status = 0)
  }

  allowed <- c(
    "--input-dir",
    "--group-combination",
    "--screen-bed",
    "--output-dir",
    "--top-threshold",
    "--star-gene",
    "--phenotype-id"
  )
  values <- list(
    input_dir = NULL,
    group_combination = NULL,
    screen_bed = NULL,
    output_dir = getwd(),
    top_threshold = "0.05",
    star_gene = NULL,
    phenotype_id = NULL
  )
  option_to_name <- c(
    "--input-dir" = "input_dir",
    "--group-combination" = "group_combination",
    "--screen-bed" = "screen_bed",
    "--output-dir" = "output_dir",
    "--top-threshold" = "top_threshold",
    "--star-gene" = "star_gene",
    "--phenotype-id" = "phenotype_id"
  )

  seen <- character(0)
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    option <- token
    value <- NULL

    if (grepl("^--[^=]+=", token)) {
      option <- sub("=.*$", "", token)
      value <- sub("^[^=]*=", "", token)
    }

    if (!(option %in% allowed)) {
      stopf("无法识别的参数：%s。请使用 --help 查看帮助。", token)
    }
    if (option %in% seen) {
      stopf("参数重复提供：%s", option)
    }
    seen <- c(seen, option)

    if (is.null(value)) {
      if (i == length(args)) {
        stopf("参数 %s 缺少取值。", option)
      }
      value <- args[[i + 1L]]
      if (startsWith(value, "--")) {
        stopf("参数 %s 缺少取值。", option)
      }
      i <- i + 1L
    }
    if (!nzchar(value)) {
      stopf("参数 %s 的取值不能为空。", option)
    }

    values[[unname(option_to_name[[option]])]] <- value
    i <- i + 1L
  }

  if (is.null(values$input_dir)) {
    stopf("缺少必需参数 --input-dir。")
  }
  if (is.null(values$group_combination)) {
    stopf("缺少必需参数 --group-combination。")
  }

  threshold <- suppressWarnings(as.numeric(values$top_threshold))
  if (length(threshold) != 1L || is.na(threshold) || !is.finite(threshold) ||
      threshold <= 0 || threshold > 1) {
    stopf("--top-threshold 必须是大于 0 且不大于 1 的有限数值。")
  }
  values$top_threshold <- threshold
  if (!is.null(values$phenotype_id)) {
    values$phenotype_id <- trimws(values$phenotype_id)
    if (!nzchar(values$phenotype_id) ||
        grepl("[/\\\\]", values$phenotype_id) ||
        values$phenotype_id %in% c(".", "..")) {
      stopf("--phenotype-id 不能为空、.、..，也不能包含路径分隔符。")
    }
    if (is.null(values$star_gene)) {
      stopf("--phenotype-id 仅能与 --star-gene 一起使用。")
    }
  }
  values
}

normalize_existing_path <- function(path, label, directory = FALSE) {
  if (directory) {
    if (!dir.exists(path)) {
      stopf("%s不存在或不是目录：%s", label, path)
    }
  } else if (!file.exists(path) || dir.exists(path)) {
    stopf("%s不存在或不是普通文件：%s", label, path)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

resolve_phenotype_id <- function(star_gene, phenotype_id) {
  if (is.null(star_gene)) {
    return(NULL)
  }
  resolved <- if (is.null(phenotype_id)) {
    tools::file_path_sans_ext(basename(star_gene))
  } else {
    phenotype_id
  }
  resolved <- trimws(resolved)
  if (!nzchar(resolved) || grepl("[/\\\\]", resolved) || resolved %in% c(".", "..")) {
    stopf(
      "无法生成有效的表型 ID；请通过 --phenotype-id 提供不含路径分隔符的值。"
    )
  }
  resolved
}

make_logger <- function(log_file) {
  cat("", file = log_file)
  function(level, fmt, ...) {
    text <- sprintf(fmt, ...)
    line <- sprintf(
      "[%s] [%s] %s",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      level,
      text
    )
    message(line)
    cat(line, "\n", file = log_file, append = TRUE, sep = "")
  }
}

read_group_combinations <- function(file) {
  dat <- tryCatch(
    read.table(
      file,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(e) stopf("无法读取亚群对列表 %s：%s", file, conditionMessage(e))
  )

  if (ncol(dat) < 2L) {
    stopf("亚群对列表至少需要两列，前两列依次为 Ref 和 Check。")
  }
  names(dat)[1L] <- sub("^\\ufeff", "", names(dat)[1L])
  if (!identical(names(dat)[1:2], c("Ref", "Check"))) {
    stopf(
      "亚群对列表前两列表头必须依次为 Ref 和 Check；当前为：%s。",
      paste(names(dat)[1:2], collapse = ", ")
    )
  }
  if (nrow(dat) == 0L) {
    stopf("亚群对列表不包含任何数据行。")
  }

  groups <- data.frame(
    Ref = trimws(as.character(dat[[1L]])),
    Check = trimws(as.character(dat[[2L]])),
    stringsAsFactors = FALSE
  )
  invalid <- !nzchar(groups$Ref) | !nzchar(groups$Check) |
    grepl("[/\\\\]", groups$Ref) | grepl("[/\\\\]", groups$Check)
  if (any(invalid)) {
    stopf(
      "亚群对列表第 %s 个数据行包含空名称或路径分隔符。",
      paste(which(invalid), collapse = ", ")
    )
  }
  duplicated_pair <- duplicated(groups[c("Ref", "Check")])
  if (any(duplicated_pair)) {
    stopf(
      "亚群对列表包含重复组合：%s",
      paste(
        paste0("Ref", groups$Ref[duplicated_pair], ".vs.Check", groups$Check[duplicated_pair]),
        collapse = ", "
      )
    )
  }
  groups
}

read_screen_bed <- function(file, log_message) {
  bed <- tryCatch(
    read.table(
      file,
      header = FALSE,
      sep = "",
      quote = "",
      comment.char = "#",
      fill = TRUE,
      stringsAsFactors = FALSE
    ),
    error = function(e) stopf("无法读取 BED 文件 %s：%s", file, conditionMessage(e))
  )
  if (nrow(bed) == 0L || ncol(bed) < 3L) {
    stopf("BED 文件必须至少包含一行和三列：chrom、start、end。")
  }

  bed_start0 <- suppressWarnings(as.numeric(bed[[2L]]))
  bed_end0 <- suppressWarnings(as.numeric(bed[[3L]]))
  invalid <- !nzchar(as.character(bed[[1L]])) |
    !is.finite(bed_start0) | !is.finite(bed_end0) |
    bed_start0 < 0 | bed_end0 <= bed_start0 |
    bed_start0 != floor(bed_start0) | bed_end0 != floor(bed_end0)
  if (any(invalid)) {
    stopf(
      "BED 文件第 %s 行坐标无效；要求 0 <= start < end。",
      paste(which(invalid), collapse = ", ")
    )
  }

  result <- data.frame(
    CHROM = as.character(bed[[1L]]),
    START_1BASED = bed_start0 + 1,
    END_1BASED = bed_end0,
    stringsAsFactors = FALSE
  )
  log_message(
    "INFO",
    "读取 BED：%s；有效屏蔽区域 %d 个，涉及 %d 条染色体。",
    basename(file),
    nrow(result),
    length(unique(result$CHROM))
  )
  result
}

find_one_method_file <- function(files, pair_prefix, suffix, method) {
  basenames <- basename(files)
  matched <- files[
    startsWith(basenames, paste0(pair_prefix, ".")) & endsWith(basenames, suffix)
  ]
  matched <- sort(matched)
  if (length(matched) == 0L) {
    return(list(error = sprintf("%s 未找到以 %s 开头且以 %s 结尾的文件。", method, pair_prefix, suffix)))
  }
  if (length(matched) > 1L) {
    return(list(error = sprintf(
      "%s 找到多个候选文件，无法唯一确定：%s",
      method,
      paste(basename(matched), collapse = ", ")
    )))
  }
  list(file = matched[[1L]])
}

discover_file_groups <- function(input_dir, groups, log_message) {
  all_files <- list.files(input_dir, full.names = TRUE, recursive = FALSE, all.files = FALSE)
  all_files <- all_files[file.exists(all_files) & !dir.exists(all_files)]
  suffixes <- c(
    PI = ".log2_pi_ratio.pi",
    FST = ".windowed.weir.fst",
    XPCLR = ".xpclr"
  )

  discovered <- vector("list", nrow(groups))
  errors <- character(0)
  for (i in seq_len(nrow(groups))) {
    prefix <- paste0("Ref", groups$Ref[[i]], ".vs.Check", groups$Check[[i]])
    matches <- lapply(names(suffixes), function(method) {
      find_one_method_file(all_files, prefix, suffixes[[method]], method)
    })
    names(matches) <- names(suffixes)
    group_errors <- vapply(matches, function(x) x$error %||% "", character(1L))

    if (any(nzchar(group_errors))) {
      for (err in group_errors[nzchar(group_errors)]) {
        log_message("ERROR", "文件组 %s：%s", prefix, err)
      }
      errors <- c(errors, sprintf("%s 文件不完整或不唯一", prefix))
      next
    }

    discovered[[i]] <- list(
      Ref = groups$Ref[[i]],
      Check = groups$Check[[i]],
      prefix = prefix,
      pi = matches$PI$file,
      fst = matches$FST$file,
      xpclr = matches$XPCLR$file
    )
    log_message("INFO", "识别文件组：%s", prefix)
    log_message("INFO", "  PI    = %s", basename(matches$PI$file))
    log_message("INFO", "  FST   = %s", basename(matches$FST$file))
    log_message("INFO", "  XPCLR = %s", basename(matches$XPCLR$file))
  }

  if (length(errors) > 0L) {
    stopf("文件识别失败：%s。详见日志。", paste(errors, collapse = "；"))
  }
  discovered
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

read_required <- function(file, columns) {
  dat <- tryCatch(
    read.table(
      file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) stopf("无法读取 %s：%s", file, conditionMessage(e))
  )
  missing_cols <- setdiff(columns, names(dat))
  if (length(missing_cols) > 0L) {
    stopf(
      "%s 缺少必需列：%s",
      basename(file),
      paste(missing_cols, collapse = ", ")
    )
  }
  dat
}

standardize_method_data <- function(file, method) {
  if (method == "PI") {
    dat <- read_required(file, c("CHROM", "BIN_START", "BIN_END", "LOG2_PI_RATIO"))
    result <- data.frame(
      CHROM = as.character(dat$CHROM),
      BIN_START = suppressWarnings(as.numeric(dat$BIN_START)),
      BIN_END = suppressWarnings(as.numeric(dat$BIN_END)),
      LOG2_PI_RATIO = suppressWarnings(as.numeric(dat$LOG2_PI_RATIO)),
      stringsAsFactors = FALSE
    )
  } else if (method == "FST") {
    dat <- read_required(file, c("CHROM", "BIN_START", "BIN_END", "WEIGHTED_FST"))
    result <- data.frame(
      CHROM = as.character(dat$CHROM),
      BIN_START = suppressWarnings(as.numeric(dat$BIN_START)),
      BIN_END = suppressWarnings(as.numeric(dat$BIN_END)),
      WEIGHTED_FST = suppressWarnings(as.numeric(dat$WEIGHTED_FST)),
      stringsAsFactors = FALSE
    )
  } else if (method == "XPCLR") {
    dat <- read_required(file, c("chrom", "start", "stop", "xpclr"))
    result <- data.frame(
      CHROM = as.character(dat$chrom),
      BIN_START = suppressWarnings(as.numeric(dat$start)),
      BIN_END = suppressWarnings(as.numeric(dat$stop)),
      XPCLR = suppressWarnings(as.numeric(dat$xpclr)),
      stringsAsFactors = FALSE
    )
  } else {
    stopf("内部错误：未知方法 %s。", method)
  }

  invalid_interval <- !nzchar(result$CHROM) |
    !is.finite(result$BIN_START) | !is.finite(result$BIN_END) |
    result$BIN_START < 0 | result$BIN_END < result$BIN_START |
    result$BIN_START != floor(result$BIN_START) |
    result$BIN_END != floor(result$BIN_END)
  if (any(invalid_interval)) {
    stopf(
      "%s 中第 %s 个数据行的染色体或区间坐标无效。",
      basename(file),
      paste(which(invalid_interval), collapse = ", ")
    )
  }

  key_cols <- c("CHROM", "BIN_START", "BIN_END")
  duplicated_key <- duplicated(result[key_cols]) | duplicated(result[key_cols], fromLast = TRUE)
  if (any(duplicated_key)) {
    stopf(
      "%s 包含 %d 行重复区间；为避免合并时产生笛卡尔积，已停止处理。",
      basename(file),
      sum(duplicated_key)
    )
  }
  result
}

find_screen_overlaps <- function(chrom, start, end, bed) {
  hit <- rep(FALSE, length(chrom))
  if (is.null(bed) || nrow(bed) == 0L || length(chrom) == 0L) {
    return(hit)
  }

  shared_chromosomes <- intersect(unique(chrom), unique(bed$CHROM))
  for (chr in shared_chromosomes) {
    window_idx <- which(chrom == chr)
    bed_chr <- bed[bed$CHROM == chr, , drop = FALSE]
    ord <- order(bed_chr$START_1BASED, bed_chr$END_1BASED)
    bed_starts <- bed_chr$START_1BASED[ord]
    cumulative_end <- cummax(bed_chr$END_1BASED[ord])

    last_start_before_window_end <- findInterval(end[window_idx], bed_starts)
    possible <- last_start_before_window_end > 0L
    if (any(possible)) {
      local_idx <- which(possible)
      hit[window_idx[local_idx]] <-
        cumulative_end[last_start_before_window_end[local_idx]] >= start[window_idx[local_idx]]
    }
  }
  hit
}

stat_to_pvalue_minotaur <- function(dfv, right_tailed) {
  df_vars <- as.matrix(dfv)
  n <- nrow(df_vars)
  d <- ncol(df_vars)
  if (n < 2L) {
    stopf("至少需要两个共有且有效的窗口。")
  }
  if (length(right_tailed) != d) {
    stopf("right_tailed 的长度必须与统计量数量相同。")
  }

  df_p <- as.data.frame(matrix(0, nrow = n, ncol = d))
  names(df_p) <- paste0(names(dfv), "_P")
  for (i in seq_len(d)) {
    p <- (rank(df_vars[, i]) - 1) / (n - 1)
    if (right_tailed[[i]]) {
      p <- 1 - p
    }
    df_p[, i] <- (p * n + 1) / (n + 2)
  }
  df_p
}

calculate_background_matrices <- function(stats_df, background_idx) {
  if (length(background_idx) < 2L) {
    stopf("BED 筛选后少于两个背景窗口，无法计算协方差矩阵。")
  }
  background <- as.matrix(stats_df[background_idx, , drop = FALSE])
  covariance <- cov(background, use = "pairwise.complete.obs")
  if (any(is.na(covariance)) || any(!is.finite(covariance))) {
    stopf("背景协方差矩阵包含 NA 或非有限值。")
  }
  if (inherits(try(solve(covariance), silent = TRUE), "try-error")) {
    stopf("背景协方差矩阵为奇异矩阵，无法按 MINOTAUR 的规则计算 DCMS。")
  }
  correlation <- covariance / sqrt(outer(diag(covariance), diag(covariance)))
  if (any(is.na(correlation)) || any(!is.finite(correlation))) {
    stopf("背景相关矩阵包含 NA 或非有限值；请检查是否存在零方差统计量。")
  }
  list(covariance = covariance, correlation = correlation)
}

calculate_dcms <- function(pvals_df, correlation) {
  df_p <- as.matrix(pvals_df)
  if (ncol(df_p) != nrow(correlation)) {
    stopf("p 值列数与背景相关矩阵维度不一致。")
  }
  dcms <- numeric(nrow(df_p))
  for (i in seq_len(ncol(df_p))) {
    dcms <- dcms +
      (log(1 - df_p[, i]) - log(df_p[, i])) /
      sum(abs(correlation[i, ]))
  }
  as.numeric(dcms)
}

format_matrix_for_log <- function(mat) {
  lines <- capture.output(
    print(round(mat, digits = 6), quote = FALSE)
  )
  paste(lines, collapse = "\n")
}

write_table_atomic <- function(dat, file) {
  tmp <- tempfile(pattern = paste0(".", basename(file), "."), tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  write.table(dat, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (!file.rename(tmp, file)) {
    stopf("无法将临时输出文件移动到目标位置：%s", file)
  }
}

star_gene_extension <- function(file) {
  extension <- tolower(tools::file_ext(file))
  supported <- c("xlsx", "xls", "txt", "tsv")
  if (!(extension %in% supported)) {
    stopf(
      "--star-gene 仅支持 xlsx、xls、txt 和 tsv 文件：%s",
      file
    )
  }
  extension
}

check_plot_dependencies <- function(star_gene) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stopf("缺少绘图必需的 R 包 ggplot2。")
  }
  if (!is.null(star_gene)) {
    extension <- star_gene_extension(star_gene)
    if (extension %in% c("xlsx", "xls") &&
        !requireNamespace("readxl", quietly = TRUE)) {
      stopf("读取 xlsx/xls 格式的 --star-gene 时需要 R 包 readxl。")
    }
  }
}

chromosome_number <- function(x) {
  suppressWarnings(as.numeric(sub("^Chr", "", trimws(as.character(x)), ignore.case = TRUE)))
}

read_gene_annotations <- function(file, log_message) {
  extension <- star_gene_extension(file)
  source_description <- basename(file)
  if (extension %in% c("xlsx", "xls")) {
    sheets <- tryCatch(
      readxl::excel_sheets(file),
      error = function(e) stopf(
        "无法读取基因 Excel 文件 %s 的工作表信息：%s",
        file,
        conditionMessage(e)
      )
    )
    if (length(sheets) == 0L) {
      stopf("基因 Excel 文件不包含任何工作表：%s", file)
    }
    selected_sheet <- if ("StarGene" %in% sheets) "StarGene" else sheets[[1L]]
    if (identical(selected_sheet, "StarGene")) {
      log_message("INFO", "已在 %s 中找到 StarGene 工作表。", basename(file))
    } else {
      log_message(
        "WARN",
        "%s 中未找到 StarGene 工作表；自动使用第一个工作表 %s。",
        basename(file),
        selected_sheet
      )
    }
    gene_data <- tryCatch(
      readxl::read_excel(file, sheet = selected_sheet),
      error = function(e) stopf(
        "无法读取基因 Excel 文件 %s 的工作表 %s：%s",
        file,
        selected_sheet,
        conditionMessage(e)
      )
    )
    source_description <- sprintf("%s（sheet=%s）", basename(file), selected_sheet)
  } else {
    first_line <- tryCatch(
      readLines(file, n = 1L, warn = FALSE),
      error = function(e) stopf("无法读取基因文本文件 %s：%s", file, conditionMessage(e))
    )
    if (length(first_line) == 0L || !nzchar(first_line[[1L]])) {
      stopf("基因文本文件为空：%s", file)
    }
    separator <- if (extension == "tsv" || grepl("\t", first_line[[1L]], fixed = TRUE)) {
      "\t"
    } else if (grepl(",", first_line[[1L]], fixed = TRUE)) {
      ","
    } else if (grepl(";", first_line[[1L]], fixed = TRUE)) {
      ";"
    } else {
      ""
    }
    gene_data <- tryCatch(
      read.table(
        file,
        header = TRUE,
        sep = separator,
        quote = "\"",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      error = function(e) stopf(
        "无法读取基因文本文件 %s：%s",
        file,
        conditionMessage(e)
      )
    )
  }
  gene_data <- as.data.frame(gene_data, stringsAsFactors = FALSE)
  if (ncol(gene_data) > 0L) {
    names(gene_data)[1L] <- sub("^\\ufeff", "", names(gene_data)[1L])
  }
  required_cols <- c("Name", "Chr", "Start", "End")
  missing_cols <- setdiff(required_cols, names(gene_data))
  if (length(missing_cols) > 0L) {
    stopf(
      "基因表缺少必需列：%s",
      paste(missing_cols, collapse = ", ")
    )
  }

  result <- data.frame(
    Chr = as.character(gene_data$Chr),
    Start = suppressWarnings(as.numeric(gene_data$Start)),
    End = suppressWarnings(as.numeric(gene_data$End)),
    Name = trimws(as.character(gene_data$Name)),
    stringsAsFactors = FALSE
  )
  result$ChrNumber <- chromosome_number(result$Chr)
  valid <- nzchar(result$Name) & is.finite(result$ChrNumber) &
    is.finite(result$Start) & is.finite(result$End) &
    result$Start >= 0 & result$End >= result$Start
  invalid_count <- sum(!valid)
  result <- result[valid, , drop = FALSE]
  if (nrow(result) == 0L) {
    stopf("基因表中没有坐标和名称均有效的基因。")
  }
  log_message(
    "INFO",
    "读取已知基因：%s；有效基因=%d，排除无效记录=%d。",
    source_description,
    nrow(result),
    invalid_count
  )
  result
}

format_top_tag <- function(top_threshold) {
  percentage <- format(
    top_threshold * 100,
    scientific = FALSE,
    trim = TRUE,
    digits = 10
  )
  percentage <- sub("(\\.[0-9]*?)0+$", "\\1", percentage)
  percentage <- sub("\\.$", "", percentage)
  paste0("top", gsub("\\.", "p", percentage), "pct")
}

prepare_manhattan_data <- function(result, comparison) {
  plot_data <- result
  plot_data$ChrNumber <- chromosome_number(plot_data$CHROM)
  invalid <- !is.finite(plot_data$ChrNumber)
  if (any(invalid)) {
    stopf(
      "%s 的 CHROM 中有 %d 个值无法解析为 Chr<数字>，无法绘制曼哈顿图。",
      comparison,
      sum(invalid)
    )
  }
  # 绘图标签保留原始染色体 ID，数值编号仅用于排序和坐标计算。
  plot_data$chrom <- as.character(plot_data$CHROM)

  chromosome_lengths <- aggregate(
    plot_data$BIN_END,
    by = list(ChrNumber = plot_data$ChrNumber),
    FUN = max,
    na.rm = TRUE
  )
  names(chromosome_lengths)[2L] <- "max_pos"
  chromosome_lengths <- chromosome_lengths[order(chromosome_lengths$ChrNumber), , drop = FALSE]
  chromosome_lengths$offset <- c(0, head(cumsum(chromosome_lengths$max_pos), -1L))
  offset_map <- setNames(chromosome_lengths$offset, chromosome_lengths$ChrNumber)
  plot_data$Pos <- plot_data$BIN_START + unname(offset_map[as.character(plot_data$ChrNumber)])
  plot_data
}

find_selected_gene_overlaps <- function(plot_data, gene_data, comparison) {
  if (is.null(gene_data)) {
    return(data.frame())
  }
  selected <- plot_data[plot_data$DCMS_SELECTED, , drop = FALSE]
  overlap_list <- vector("list", nrow(gene_data))
  overlap_count <- 0L
  for (i in seq_len(nrow(gene_data))) {
    hit <- selected$ChrNumber == gene_data$ChrNumber[[i]] &
      gene_data$Start[[i]] <= selected$BIN_END &
      gene_data$End[[i]] >= selected$BIN_START
    if (any(hit)) {
      overlap_count <- overlap_count + 1L
      overlap <- selected[hit, , drop = FALSE]
      overlap$GeneName <- gene_data$Name[[i]]
      overlap$GeneStart <- gene_data$Start[[i]]
      overlap$GeneEnd <- gene_data$End[[i]]
      overlap$COMPARISON <- comparison
      overlap_list[[overlap_count]] <- overlap
    }
  }
  if (overlap_count == 0L) {
    return(data.frame())
  }
  do.call(rbind, overlap_list[seq_len(overlap_count)])
}

plot_dcms_result <- function(
    result,
    comparison,
    output_dir,
    top_threshold,
    gene_data,
    phenotype_id,
    log_message) {
  plot_data <- prepare_manhattan_data(result, comparison)
  selected_data <- plot_data[plot_data$DCMS_SELECTED, , drop = FALSE]
  if (nrow(selected_data) == 0L) {
    stopf("%s 没有 DCMS_SELECTED 窗口，无法确定绘图阈值。", comparison)
  }
  dcms_threshold <- min(selected_data$DCMS)
  selected_gene_windows <- find_selected_gene_overlaps(
    plot_data,
    gene_data,
    comparison
  )
  plot_data$KnownGeneSelected <- FALSE
  if (nrow(selected_gene_windows) > 0L) {
    key_all <- paste(plot_data$CHROM, plot_data$BIN_START, plot_data$BIN_END, sep = "_")
    key_selected <- paste(
      selected_gene_windows$CHROM,
      selected_gene_windows$BIN_START,
      selected_gene_windows$BIN_END,
      sep = "_"
    )
    plot_data$KnownGeneSelected <- key_all %in% key_selected
    label_order <- order(
      selected_gene_windows$GeneName,
      -selected_gene_windows$DCMS
    )
    label_candidates <- selected_gene_windows[label_order, , drop = FALSE]
    gene_labels_plot <- label_candidates[
      !duplicated(label_candidates$GeneName),
      ,
      drop = FALSE
    ]
  } else {
    gene_labels_plot <- data.frame()
  }

  chromosome_levels <- unique(plot_data$chrom[order(plot_data$ChrNumber)])
  chromosome_colours <- setNames(
    rep(c("#579D1C", "#FF950E"), length.out = length(chromosome_levels)),
    chromosome_levels
  )
  midpoint <- aggregate(
    plot_data$Pos,
    by = list(chrom = plot_data$chrom),
    FUN = median
  )
  names(midpoint)[2L] <- "Chr_midpoint"
  midpoint$ChrNumber <- chromosome_number(midpoint$chrom)
  midpoint <- midpoint[order(midpoint$ChrNumber), , drop = FALSE]

  y_range <- range(plot_data$DCMS, na.rm = TRUE)
  y_padding <- diff(y_range) * 0.08
  if (!is.finite(y_padding) || y_padding == 0) {
    y_padding <- 1
  }
  y_lower_limit <- y_range[[1L]] - y_padding * 0.2
  y_upper_limit <- y_range[[2L]] + y_padding
  regular_y_breaks <- pretty(y_range, n = 5)
  regular_y_breaks <- regular_y_breaks[
    regular_y_breaks >= y_lower_limit & regular_y_breaks <= y_upper_limit
  ]
  if (length(regular_y_breaks) > 0L) {
    nearest_break <- which.min(abs(regular_y_breaks - dcms_threshold))
    regular_y_breaks <- regular_y_breaks[-nearest_break]
  }
  y_breaks <- sort(unique(c(regular_y_breaks, dcms_threshold)))
  y_labels <- formatC(y_breaks, format = "fg", digits = 4)
  threshold_label_index <- which.min(abs(y_breaks - dcms_threshold))
  y_labels[[threshold_label_index]] <- formatC(dcms_threshold,format = "f",digits = 2)

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = plot_data,ggplot2::aes(x = Pos / 1e6, y = DCMS, color = chrom),size = 0.35,shape=16) +
    ggplot2::geom_point(data = plot_data[plot_data$KnownGeneSelected, , drop = FALSE],ggplot2::aes(x = Pos / 1e6, y = DCMS),color = "red",size = 0.7,shape=16) +
    ggplot2::geom_hline(yintercept = dcms_threshold,linetype = 2,color = "black",linewidth = 0.25) +
    ggplot2::scale_color_manual(values = chromosome_colours) +
    ggplot2::labs(x = "", y = "DCMS") +
    ggplot2::scale_x_continuous(labels = midpoint$chrom,breaks = midpoint$Chr_midpoint / 1e6,expand = c(0.005, 0)) +
    ggplot2::scale_y_continuous(breaks = y_breaks,labels = y_labels,expand = c(0.005, 0),limits = c(y_lower_limit, y_upper_limit)) +
    ggplot2::theme_bw(base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black",fill = NA,linewidth = 0.25),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(linewidth = 0.25),
	  axis.ticks.length = ggplot2::unit(0.05,"cm"),
      axis.text.x = ggplot2::element_text(family = "sans",size = 6,colour = "black"),
      axis.text.y = ggplot2::element_text(family = "sans",size = 7,colour = "black"),
      axis.title.y = ggplot2::element_text(family = "sans",size = 7,margin = ggplot2::margin(r = 2.5)),
      legend.position = "none"
    )

  if (nrow(gene_labels_plot) > 0L) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = gene_labels_plot,
        ggplot2::aes(x = Pos / 1e6, y = DCMS, label = GeneName),
        color = "black",
        family = "sans",
        fontface = "italic",
        size = 2.2,
        segment.color = "grey40",
        min.segment.length = 0,
        max.overlaps = 100,
        box.padding = 0.3,
        point.padding = 0.15,
		force = 1,
		force_pull = 0.5,
		segment.size = 0.2,
		max.iter = 5000,
		max.time = 5
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = gene_labels_plot,
        ggplot2::aes(x = Pos / 1e6, y = DCMS, label = GeneName),
        color = "black",
        family = "sans",
        fontface = "italic",
        size = 1.5,
        vjust = -0.6,
        check_overlap = TRUE
      )
    }
  }

  top_tag <- format_top_tag(top_threshold)
  annotation_tag <- if (is.null(gene_data)) "" else paste0(".", phenotype_id)
  output_stem <- file.path(
    output_dir,
    paste0(comparison, ".DCMS.", top_tag, annotation_tag)
  )
  output_png <- paste0(output_stem, ".png")
  output_pdf <- paste0(output_stem, ".pdf")
  ggplot2::ggsave(
    filename = output_png,
    plot = p,
    device = "png",
    width = 4,
    height = 1.2,
    dpi = 600
  )
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(
    filename = output_pdf,
    plot = p,
    device = pdf_device,
    width = 4,
    height = 1.2
  )

  overlap_file <- NA_character_
  if (nrow(selected_gene_windows) > 0L) {
    overlap_output <- selected_gene_windows[
      , c(
        "COMPARISON", "CHROM", "BIN_START", "BIN_END", "DCMS",
        "GeneName", "GeneStart", "GeneEnd"
      ),
      drop = FALSE
    ]
    overlap_file <- paste0(output_stem, ".overlap.tsv")
    write_table_atomic(overlap_output, overlap_file)
  }

  unique_gene_windows <- if (nrow(selected_gene_windows) == 0L) {
    0L
  } else {
    length(unique(paste(
      selected_gene_windows$CHROM,
      selected_gene_windows$BIN_START,
      selected_gene_windows$BIN_END,
      sep = "_"
    )))
  }
  unique_genes <- if (nrow(selected_gene_windows) == 0L) {
    0L
  } else {
    length(unique(selected_gene_windows$GeneName))
  }
  log_message(
    "INFO",
    paste0(
      "%s 绘图完成：入选阈值 DCMS >= %.6g，入选窗口=%d，",
      "重叠已知基因窗口=%d，已知基因=%d；PNG=%s；PDF=%s"
    ),
    comparison,
    dcms_threshold,
    nrow(selected_data),
    unique_gene_windows,
    unique_genes,
    output_png,
    output_pdf
  )

  summary <- data.frame(
    COMPARISON = comparison,
    WINDOWS = nrow(plot_data),
    DCMS_TOP_THRESHOLD = top_threshold,
    DCMS_SELECTION_CUTOFF = dcms_threshold,
    SELECTED_WINDOWS = nrow(selected_data),
    SELECTED_STARGENE_WINDOWS = unique_gene_windows,
    SELECTED_STARGENES = unique_genes,
    PNG = output_png,
    PDF = output_pdf,
    stringsAsFactors = FALSE
  )
  list(
    summary = summary,
    gene_overlaps = selected_gene_windows,
    overlap_file = overlap_file,
    png = output_png,
    pdf = output_pdf
  )
}

process_group <- function(
    group,
    screen_bed,
    output_dir,
    top_threshold,
    gene_data,
    phenotype_id,
    log_message) {
  log_message("INFO", "开始处理：%s", group$prefix)

  pi_dat <- standardize_method_data(group$pi, "PI")
  fst_dat <- standardize_method_data(group$fst, "FST")
  xpclr_dat <- standardize_method_data(group$xpclr, "XPCLR")
  log_message(
    "INFO",
    "%s 原始区间数：PI=%d，FST=%d，XPCLR=%d。",
    group$prefix,
    nrow(pi_dat),
    nrow(fst_dat),
    nrow(xpclr_dat)
  )

  key_cols <- c("CHROM", "BIN_START", "BIN_END")
  merged <- merge(pi_dat, fst_dat, by = key_cols, all = FALSE, sort = FALSE)
  merged <- merge(merged, xpclr_dat, by = key_cols, all = FALSE, sort = FALSE)
  common_before_finite <- nrow(merged)

  all_keys <- unique(rbind(
    pi_dat[key_cols],
    fst_dat[key_cols],
    xpclr_dat[key_cols]
  ))
  log_message(
    "INFO",
    paste0(
      "%s 三种方法共有区间=%d；各方法因不共有而排除：",
      "PI=%d，FST=%d，XPCLR=%d；全部方法区间并集=%d，其中非三者共有=%d。"
    ),
    group$prefix,
    common_before_finite,
    nrow(pi_dat) - common_before_finite,
    nrow(fst_dat) - common_before_finite,
    nrow(xpclr_dat) - common_before_finite,
    nrow(all_keys),
    nrow(all_keys) - common_before_finite
  )

  finite_rows <- is.finite(merged$LOG2_PI_RATIO) &
    is.finite(merged$WEIGHTED_FST) &
    is.finite(merged$XPCLR)
  dropped_nonfinite <- sum(!finite_rows)
  merged <- merged[finite_rows, , drop = FALSE]
  log_message(
    "INFO",
    "%s 共有区间中非有限统计值排除=%d；正式分析保留=%d。",
    group$prefix,
    dropped_nonfinite,
    nrow(merged)
  )
  if (nrow(merged) < 2L) {
    stopf("%s 仅剩 %d 个共有有效窗口，无法计算 DCMS。", group$prefix, nrow(merged))
  }

  screen_hit <- find_screen_overlaps(
    merged$CHROM,
    merged$BIN_START,
    merged$BIN_END,
    screen_bed
  )
  background_idx <- which(!screen_hit)
  if (is.null(screen_bed)) {
    log_message(
      "INFO",
      "%s 未提供 BED；协方差背景使用全部 %d 个共有有效窗口。",
      group$prefix,
      length(background_idx)
    )
  } else {
    log_message(
      "INFO",
      paste0(
        "%s BED 与 %d 个共有有效窗口存在交集；这些窗口仅从协方差背景排除。",
        "协方差背景保留=%d，正式 DCMS 分析仍包含=%d。"
      ),
      group$prefix,
      sum(screen_hit),
      length(background_idx),
      nrow(merged)
    )
  }

  stats_df <- data.frame(
    FST = merged$WEIGHTED_FST,
    LOG2_PI_RATIO = merged$LOG2_PI_RATIO,
    XPCLR = merged$XPCLR,
    stringsAsFactors = FALSE
  )
  matrices <- calculate_background_matrices(stats_df, background_idx)
  log_message(
    "INFO",
    "%s 背景协方差矩阵：\n%s",
    group$prefix,
    format_matrix_for_log(matrices$covariance)
  )
  log_message(
    "INFO",
    "%s 背景相关矩阵：\n%s",
    group$prefix,
    format_matrix_for_log(matrices$correlation)
  )

  pvals_df <- stat_to_pvalue_minotaur(
    stats_df,
    right_tailed = c(TRUE, TRUE, TRUE)
  )
  dcms <- calculate_dcms(pvals_df, matrices$correlation)

  rank_desc <- rank(-dcms, ties.method = "min")
  nominal_top_n <- ceiling(length(dcms) * top_threshold)
  selected <- rank_desc <= nominal_top_n

  result <- data.frame(
    Ref = group$Ref,
    Check = group$Check,
    COMPARISON = group$prefix,
    merged,
    SCREENED_FROM_COVARIANCE = screen_hit,
    pvals_df,
    DCMS = dcms,
    DCMS_RANK_DESC = rank_desc,
    DCMS_PERCENTILE = 1 - (rank_desc - 1) / length(dcms),
    DCMS_SELECTED = selected,
    DCMS_TOP_THRESHOLD = top_threshold,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  result <- result[order(result$CHROM, result$BIN_START, result$BIN_END), , drop = FALSE]

  output_file <- file.path(output_dir, paste0(group$prefix, ".DCMS.list"))
  write_table_atomic(result, output_file)
  log_message(
    "INFO",
    paste0(
      "%s DCMS 完成：分析窗口=%d，阈值=top %.6g，名义窗口数=%d，",
      "实际标记=%d；输出=%s"
    ),
    group$prefix,
    nrow(result),
    top_threshold,
    nominal_top_n,
    sum(selected),
    output_file
  )

  plot_result <- plot_dcms_result(
    result = result,
    comparison = group$prefix,
    output_dir = output_dir,
    top_threshold = top_threshold,
    gene_data = gene_data,
    phenotype_id = phenotype_id,
    log_message = log_message
  )
  invisible(list(dcms_file = output_file, plot = plot_result))
}

main <- function() {
  opts <- parse_arguments(commandArgs(trailingOnly = TRUE))
  opts$input_dir <- normalize_existing_path(opts$input_dir, "输入目录", directory = TRUE)
  opts$group_combination <- normalize_existing_path(
    opts$group_combination,
    "亚群对列表",
    directory = FALSE
  )
  if (!is.null(opts$screen_bed)) {
    opts$screen_bed <- normalize_existing_path(opts$screen_bed, "BED 文件", directory = FALSE)
  }
  if (!is.null(opts$star_gene)) {
    opts$star_gene <- normalize_existing_path(
      opts$star_gene,
      "已知基因文件",
      directory = FALSE
    )
  }
  opts$phenotype_id <- resolve_phenotype_id(opts$star_gene, opts$phenotype_id)

  if (!dir.exists(opts$output_dir) &&
      !dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stopf("无法创建输出目录：%s", opts$output_dir)
  }
  opts$output_dir <- normalizePath(opts$output_dir, winslash = "/", mustWork = TRUE)
  log_file <- file.path(opts$output_dir, "SelectiveSweep.DCMS.Calculate.Plot.log")
  log_message <- make_logger(log_file)
  active_logger <<- log_message

  log_message("INFO", "%s 启动。", script_name)
  log_message("INFO", "输入目录：%s", opts$input_dir)
  log_message("INFO", "亚群对列表：%s", opts$group_combination)
  log_message("INFO", "BED：%s", opts$screen_bed %||% "未提供")
  log_message("INFO", "输出目录：%s", opts$output_dir)
  log_message("INFO", "Top 阈值：%.6g", opts$top_threshold)
  log_message("INFO", "已知基因文件：%s", opts$star_gene %||% "未提供（仅绘制 DCMS 曼哈顿图）")
  log_message("INFO", "表型 ID：%s", opts$phenotype_id %||% "未使用")

  check_plot_dependencies(opts$star_gene)
  log_message(
    "INFO",
    "绘图字体：sans；ggrepel=%s。",
    if (requireNamespace("ggrepel", quietly = TRUE)) "可用" else "不可用，将使用 ggplot2::geom_text"
  )

  groups <- read_group_combinations(opts$group_combination)
  log_message("INFO", "亚群对列表包含 %d 个唯一比较。", nrow(groups))
  screen_bed <- if (is.null(opts$screen_bed)) {
    NULL
  } else {
    read_screen_bed(opts$screen_bed, log_message)
  }
  gene_data <- if (is.null(opts$star_gene)) {
    NULL
  } else {
    read_gene_annotations(opts$star_gene, log_message)
  }
  file_groups <- discover_file_groups(opts$input_dir, groups, log_message)

  outputs <- vector("list", length(file_groups))
  output_index <- 0L
  for (group in file_groups) {
    output_index <- output_index + 1L
    outputs[[output_index]] <- process_group(
      group,
      screen_bed,
      opts$output_dir,
      opts$top_threshold,
      gene_data,
      opts$phenotype_id,
      log_message
    )
  }

  plot_summary <- do.call(
    rbind,
    lapply(outputs, function(x) x$plot$summary)
  )
  summary_file <- file.path(
    opts$output_dir,
    paste0("DCMS.", format_top_tag(opts$top_threshold), ".plot.summary.tsv")
  )
  write_table_atomic(plot_summary, summary_file)
  log_message("INFO", "绘图汇总表：%s", summary_file)

  overlap_tables <- lapply(outputs, function(x) x$plot$gene_overlaps)
  overlap_tables <- overlap_tables[vapply(overlap_tables, nrow, integer(1L)) > 0L]
  if (length(overlap_tables) > 0L) {
    all_overlaps <- do.call(rbind, overlap_tables)
    all_overlaps <- all_overlaps[
      , c(
        "COMPARISON", "CHROM", "BIN_START", "BIN_END", "DCMS",
        "GeneName", "GeneStart", "GeneEnd"
      ),
      drop = FALSE
    ]
    combined_overlap_file <- file.path(
      opts$output_dir,
      paste0(
        "all_comparisons.DCMS.",
        format_top_tag(opts$top_threshold),
        ".",
        opts$phenotype_id,
        ".overlap.tsv"
      )
    )
    write_table_atomic(all_overlaps, combined_overlap_file)
    log_message("INFO", "所有比较的已知基因重叠汇总：%s", combined_overlap_file)
  }

  log_message(
    "INFO",
    "%s 全部完成：成功处理 %d 个亚群对，生成 %d 个 DCMS 文件、%d 张 PNG 和 %d 份 PDF。",
    script_name,
    length(file_groups),
    length(outputs),
    length(outputs),
    length(outputs)
  )
  invisible(outputs)
}

tryCatch(
  main(),
  error = function(e) {
    if (is.function(active_logger)) {
      active_logger("ERROR", "%s", conditionMessage(e))
    } else {
      message(sprintf("[ERROR] %s", conditionMessage(e)))
    }
    quit(save = "no", status = 1L)
  }
)
