#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

read_tsv <- function(path) {
  read.delim(path, header = TRUE, sep = "\t", check.names = FALSE,
             stringsAsFactors = FALSE, na.strings = "NA", quote = "", comment.char = "")
}

format_scientific <- function(x) {
  value <- suppressWarnings(as.numeric(x))
  output <- rep("NA", length(value))
  output[is.infinite(value) & value > 0] <- "Inf"
  output[is.infinite(value) & value < 0] <- "-Inf"
  finite <- is.finite(value)
  output[finite] <- sprintf("%.2e", value[finite])
  output[finite & value == 0] <- "0"
  output
}

format_decimal <- function(x) {
  value <- suppressWarnings(as.numeric(x))
  output <- rep("NA", length(value))
  output[is.infinite(value) & value > 0] <- "Inf"
  output[is.infinite(value) & value < 0] <- "-Inf"
  finite <- is.finite(value)
  output[finite] <- sprintf("%.2f", value[finite])
  output[finite] <- sub("[.]?0+$", "", output[finite])
  output[finite & value == 0] <- "0"
  output
}

write_result <- function(result, output_file) {
  decimal_columns <- grep("Mean_DSRperKb$|\\.FC$", names(result), value = TRUE)
  probability_columns <- grep("\\.Welch\\.t$|\\.Welch\\.t\\.FDR$", names(result), value = TRUE)
  for (column in decimal_columns) result[[column]] <- format_decimal(result[[column]])
  for (column in probability_columns) result[[column]] <- format_scientific(result[[column]])
  write.table(result, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "NA")
}

welch_p_value <- function(n1, mean1, variance1, n2, mean2, variance2) {
  if (any(!is.finite(c(n1, mean1, variance1, n2, mean2, variance2))) || n1 < 2 || n2 < 2) return(NA_real_)
  variance1 <- max(variance1, 0)
  variance2 <- max(variance2, 0)
  term1 <- variance1 / n1
  term2 <- variance2 / n2
  standard_error_squared <- term1 + term2
  if (standard_error_squared <= .Machine$double.eps) {
    return(if (isTRUE(all.equal(mean1, mean2, tolerance = 1e-12))) 1 else 0)
  }
  denominator <- term1^2 / (n1 - 1) + term2^2 / (n2 - 1)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  degrees_of_freedom <- standard_error_squared^2 / denominator
  statistic <- (mean1 - mean2) / sqrt(standard_error_squared)
  2 * stats::pt(-abs(statistic), df = degrees_of_freedom)
}

fold_change <- function(numerator_mean, denominator_mean) {
  numerator_mean <- suppressWarnings(as.numeric(numerator_mean))
  denominator_mean <- suppressWarnings(as.numeric(denominator_mean))
  output <- rep(NA_real_, length(numerator_mean))
  valid <- is.finite(numerator_mean) & is.finite(denominator_mean) & numerator_mean >= 0 & denominator_mean >= 0
  ordinary <- valid & denominator_mean > 0
  output[ordinary] <- numerator_mean[ordinary] / denominator_mean[ordinary]
  output[valid & denominator_mean == 0 & numerator_mean > 0] <- Inf
  output
}

if (length(args) == 5L && identical(args[[1L]], "--raw")) {
  summary_file <- args[[2L]]
  group_order_file <- args[[3L]]
  comparison_group_file <- args[[4L]]
  output_file <- args[[5L]]
  groups <- readLines(group_order_file, warn = FALSE)
  groups <- groups[nzchar(groups)]
  if (length(groups) < 2L || anyDuplicated(groups)) stop("Group顺序列表至少需要两个不重复Group")
  comparison_groups <- read.delim(comparison_group_file, header = FALSE, sep = "\t", check.names = FALSE,
                                  stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (ncol(comparison_groups) != 2L || nrow(comparison_groups) < 1L) stop("标准化比较组表必须是至少一行的Ref/Check两列表")
  names(comparison_groups) <- c("Ref", "Check")
  if (any(!comparison_groups$Ref %in% groups) || any(!comparison_groups$Check %in% groups)) {
    stop("标准化比较组表包含不在Group顺序列表中的Group")
  }
  if (any(comparison_groups$Ref == comparison_groups$Check)) stop("比较组Ref与Check不能相同")
  dat <- read_tsv(summary_file)
  gene_columns <- c("GeneID", "Chr", "Start", "End", "RegionLength")
  if (!all(gene_columns %in% names(dat))) stop("组间统计中间表缺少基因信息列")
  result <- dat[gene_columns]
  for (group in groups) {
    mean_column <- paste0(group, ".Mean_DSRperKb")
    if (!mean_column %in% names(dat)) stop("组间统计中间表缺少列: ", mean_column)
    result[[mean_column]] <- dat[[mean_column]]
  }
  for (row in seq_len(nrow(comparison_groups))) {
    ref <- comparison_groups$Ref[[row]]
    check <- comparison_groups$Check[[row]]
    ref_required <- c(
      paste0(ref, ".ValidPairCount"), paste0(ref, ".Mean_DSRperKb"), paste0(ref, ".Variance_DSRperKb")
    )
    check_required <- c(
      paste0(check, ".ValidPairCount"), paste0(check, ".Mean_DSRperKb"), paste0(check, ".Variance_DSRperKb")
    )
    required <- c(ref_required, check_required)
    if (!all(required %in% names(dat))) stop("组间统计中间表缺少Welch t检验所需列: ", paste(setdiff(required, names(dat)), collapse = ","))

    result[[paste0(ref, ".vs.", check, ".FC")]] <- fold_change(dat[[ref_required[[2L]]]], dat[[check_required[[2L]]]])
    result[[paste0(ref, ".vs.", check, ".Welch.t")]] <- mapply(
      welch_p_value,
      dat[[ref_required[[1L]]]], dat[[ref_required[[2L]]]], dat[[ref_required[[3L]]]],
      dat[[check_required[[1L]]]], dat[[check_required[[2L]]]], dat[[check_required[[3L]]]],
      USE.NAMES = FALSE
    )
  }
  write_result(result, output_file)
  cat(sprintf("完成逐基因Group内均值及指定Group间FC/双侧Welch t检验：genes=%d, groups=%d, specified_comparisons=%d, tests=%d\n",
              nrow(result), length(groups), nrow(comparison_groups), nrow(comparison_groups)))
  quit(save = "no", status = 0L)
}

if (length(args) == 3L && identical(args[[1L]], "--fdr")) {
  raw_file <- args[[2L]]
  output_file <- args[[3L]]
  raw_result <- read_tsv(raw_file)
  p_columns <- grep("\\.Welch\\.t$", names(raw_result), value = TRUE)
  if (length(p_columns) < 1L) stop("原始Group比较表中没有.Welch.t列")
  gene_columns <- c("GeneID", "Chr", "Start", "End", "RegionLength")
  mean_columns <- grep("Mean_DSRperKb$", names(raw_result), value = TRUE)
  result <- raw_result[c(gene_columns, mean_columns)]
  for (p_column in p_columns) {
    prefix <- sub("\\.Welch\\.t$", "", p_column)
    fc_column <- paste0(prefix, ".FC")
    if (!fc_column %in% names(raw_result)) stop("原始Group比较表缺少与p-value对应的FC列: ", fc_column)
    p_values <- suppressWarnings(as.numeric(raw_result[[p_column]]))
    invalid <- !is.na(p_values) & (!is.finite(p_values) | p_values < 0 | p_values > 1)
    if (any(invalid)) stop("原始Welch t检验列包含[0,1]之外的p-value: ", p_column)
    result[[fc_column]] <- raw_result[[fc_column]]
    result[[paste0(p_column, ".FDR")]] <- stats::p.adjust(p_values, method = "BH")
  }
  write_result(result, output_file)
  cat(sprintf("从原始p-value表按每个Group比较独立完成BH-FDR校正：genes=%d, comparisons=%d\n",
              nrow(result), length(p_columns)))
  quit(save = "no", status = 0L)
}

stop("用法: ggComp.plus.group.compare.R --raw <summary.tsv> <group-order.list> <comparison-group.tsv> <raw-output.tsv>；或--fdr <raw-output.tsv> <fdr-output.tsv>")
