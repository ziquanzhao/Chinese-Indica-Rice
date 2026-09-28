#!/usr/bin/env Rscript

# 绘制 f4 值散点误差图；默认读取当前目录的 f4.xlsx。
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args)) args[[1]] else "f4.xlsx"
out_prefix <- if (length(args) >= 2) args[[2]] else sub("\\.[^.]+$", "", basename(input_file))
if (!file.exists(input_file)) stop("找不到输入文件：", input_file)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("需要 R 包 ggplot2")
# 让无用户字体缓存目录的 WSL2 环境也能安静地导出图形。

# 读取 xlsx 的第一个工作表。返回原始单元格内容，稍后自动定位真正的表头。
read_first_xlsx_sheet <- function(file) {
  if (requireNamespace("readxl", quietly = TRUE)) {
    return(as.data.frame(readxl::read_excel(
      file, sheet = 1, col_names = FALSE, col_types = "text"
    ), stringsAsFactors = FALSE))
  }
  if (!requireNamespace("xml2", quietly = TRUE)) stop("读取 xlsx 需要 R 包 readxl 或 xml2")
  entries <- utils::unzip(file, list = TRUE)$Name
  sheet_file <- "xl/worksheets/sheet1.xml"
  if (!sheet_file %in% entries) stop("xlsx 中找不到第一个工作表：", sheet_file)

  shared_strings <- character()
  if ("xl/sharedStrings.xml" %in% entries) {
    shared_xml <- xml2::read_xml(unz(file, "xl/sharedStrings.xml"))
    shared_strings <- xml2::xml_text(xml2::xml_find_all(shared_xml, ".//*[local-name()='si']"))
  }
  sheet_xml <- xml2::read_xml(unz(file, sheet_file))
  rows <- xml2::xml_find_all(sheet_xml, ".//*[local-name()='sheetData']/*[local-name()='row']")
  if (length(rows) < 2) stop("第一个工作表没有数据行")

  column_number <- function(reference) {
    letters <- strsplit(gsub("[0-9]", "", reference), "")[[1]]
    Reduce(function(value, letter) value * 26L + match(letter, LETTERS), letters, init = 0L)
  }
  cells <- lapply(rows, xml2::xml_find_all, xpath = "./*[local-name()='c']")
  max_column <- max(vapply(unlist(cells, recursive = FALSE), function(cell) {
    column_number(xml2::xml_attr(cell, "r"))
  }, integer(1)))
  values <- matrix(NA_character_, nrow = length(rows), ncol = max_column)
  for (row_index in seq_along(cells)) {
    for (cell in cells[[row_index]]) {
      col_index <- column_number(xml2::xml_attr(cell, "r"))
      cell_type <- xml2::xml_attr(cell, "t")
      value <- xml2::xml_text(xml2::xml_find_first(cell, "./*[local-name()='v']"))
      if (identical(cell_type, "s")) {
        value <- shared_strings[as.integer(value) + 1L]
      } else if (identical(cell_type, "inlineStr")) {
        value <- xml2::xml_text(cell)
      }
      values[row_index, col_index] <- value
    }
  }
  as.data.frame(values, stringsAsFactors = FALSE)
}

required_cols <- c("pop1", "pop2", "pop3", "pop4", "est", "se", "z", "p")

# Excel 文件可能在正式表格前包含空行、类型说明或旧表片段。
# 寻找同时包含全部必要字段的最后一个表头，以其后的非空行为数据。
extract_f4_table <- function(raw_df) {
  normalized <- as.data.frame(lapply(raw_df, function(column) {
    value <- trimws(as.character(column))
    value[value == ""] <- NA_character_
    value
  }), stringsAsFactors = FALSE)

  header_rows <- which(vapply(seq_len(nrow(normalized)), function(row_index) {
    row_values <- unname(unlist(normalized[row_index, , drop = TRUE], use.names = FALSE))
    all(required_cols %in% row_values)
  }, logical(1)))
  if (!length(header_rows)) {
    stop("找不到包含以下字段的完整表头：", paste(required_cols, collapse = ", "))
  }

  header_row <- tail(header_rows, 1)
  header_values <- unname(unlist(normalized[header_row, , drop = TRUE], use.names = FALSE))
  column_index <- match(required_cols, header_values)
  if (header_row == nrow(normalized)) stop("完整表头后没有数据行")

  result <- normalized[(header_row + 1L):nrow(normalized), column_index, drop = FALSE]
  names(result) <- required_cols
  result <- result[rowSums(!is.na(result)) > 0L, , drop = FALSE]
  rownames(result) <- NULL
  result
}

extension <- tolower(tools::file_ext(input_file))
df <- switch(extension,
  xlsx = extract_f4_table(read_first_xlsx_sheet(input_file)),
  tsv = utils::read.delim(input_file, check.names = FALSE),
  stop("仅支持 .xlsx 或 .tsv 输入文件")
)
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols)) stop("缺少必要列：", paste(missing_cols, collapse = ", "))
if (nrow(df) == 0L) stop("输入文件没有数据行")

pop_map <- c(
  XianLandrace1 = "L1", XianLandrace2 = "L2",
  XianCultivar1 = "C1", XianCultivar2 = "C2", OutGroup = "O"
)
pop_cols <- c("pop1", "pop2", "pop3", "pop4")
unknown_pops <- setdiff(unique(unlist(df[pop_cols], use.names = FALSE)), names(pop_map))
if (length(unknown_pops)) stop("未知群体名称：", paste(unknown_pops, collapse = ", "))
for (column in pop_cols) df[[column]] <- unname(pop_map[df[[column]]])
for (column in c("est", "se", "z", "p")) {
  df[[column]] <- as.numeric(df[[column]])
  if (any(!is.finite(df[[column]]))) stop("列 ", column, " 含非数值或缺失值")
}
if (any(df$se < 0)) stop("se 不能为负数")
if (any(df$p < 0 | df$p > 1)) stop("p 必须在 0 到 1 之间")

# 保持输入文件的原始行序，并按实际数据行数绘制。
df$row_id <- seq_len(nrow(df))
df$lower <- df$est - df$se
df$upper <- df$est + df$se
df$group_label <- sprintf("f4(%s,%s;%s,%s)", df$pop1, df$pop2, df$pop3, df$pop4)
df$significance <- ifelse(df$p < 0.001, "***", ifelse(df$p > 0.05, "ns", ""))
value_span <- diff(range(c(0, df$lower, df$upper)))
df$label_y <- if (value_span > 0) -value_span * 0.10 else -0.01
df$annotation <- sprintf("|Z|=%.2f\n%s", abs(df$z), df$significance)

p <- ggplot2::ggplot(df, ggplot2::aes(x = row_id, y = est)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color="black",linewidth = 0.5) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), color = "black", width = 0.25, linewidth = 0.5) +
  ggplot2::geom_point(color = "black", size = 2.5) +
  ggplot2::geom_text(
    ggplot2::aes(y = label_y, label = annotation),
    color = "black", size = 3.2, angle = 0,
    lineheight = 0.85, hjust = 0.5, vjust = 0.5,
    family = "DejaVu Sans"
  ) +
  ggplot2::scale_x_continuous(
    breaks = df$row_id, labels = df$group_label,
    expand = ggplot2::expansion(add = c(0.6, 0.6))
  ) +
  ggplot2::scale_y_continuous(
    breaks = function(limits) {
      candidate_breaks <- pretty(c(0, limits[2]))
      candidate_breaks[candidate_breaks >= 0]
    },
    expand = ggplot2::expansion(mult = c(0.12, 0.10))
  ) +
  ggplot2::labs(x = NULL, y = "f4 statistic") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 15, angle = 45,hjust = 1, vjust = 1,color="black",family = "DejaVu Sans"),
    axis.text.y = ggplot2::element_text(size = 18,color="black",family = "DejaVu Sans"),
    axis.title = ggplot2::element_text(size = 18,color="black",family = "DejaVu Sans"))

output_pdf <- paste0(out_prefix, "_scatter_errorbar.pdf")
output_png <- paste0(out_prefix, "_scatter_errorbar.png")
ggplot2::ggsave(output_pdf, p, width = 6, height = 5, units = "in", device = grDevices::cairo_pdf)
ggplot2::ggsave(output_png, p, width = 6, height = 5, units = "in",
                dpi = 600, device = grDevices::png, type = "cairo")
message("已生成：", output_pdf, "；", output_png)
