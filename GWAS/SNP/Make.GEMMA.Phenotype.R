#!/usr/bin/env Rscript

rm(list = ls())

# ===================== 加载包 =====================
required_packages <- c("readxl", "optparse")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# ===================== 参数定义 =====================
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              metavar = "FILE",
              help = "输入 Excel 文件路径（必需）"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              metavar = "DIR",
              help = "输出目录 [default: 当前目录]")
)

parser <- OptionParser(
  usage = "Rscript %prog -i FILE [-o DIR]",
  description = paste0(
    "\033[31m",
    paste(
      "1. 将 Excel 工作簿中的表型数据转换为 GEMMA 可用的单性状文件。",
      "2. 脚本会依次遍历工作簿内所有 Sheet；每个 Sheet 的第一列视为样本 ID，从第二列开始，每个性状分别输出为一个不含表头的单列 .GEMMA.list 文件。",
      "3. 输出文件名格式为：Sheet名称-性状名.GEMMA.list。",
      "4. 特别需要注意：Excel 工作簿中的样本顺序务必和基因型数据中的样本顺序要不多不少完全一致。",
      sep = "\n"
    ),
    "\033[0m"
  ),
  epilogue = paste(
    "使用示例：",
    "  # 输出到当前目录",
    "  Rscript Make.GEMMA.Phenotype.R -i phenotype.xlsx",
    "",
    "  # 输出到指定目录",
    "  Rscript Make.GEMMA.Phenotype.R -i phenotype.xlsx -o GEMMA_Phenotype",
    sep = "\n"
  ),
  option_list = option_list
)

opt <- parse_args(parser)

if (is.null(opt$input)) {
  print_help(parser)
  stop("必须提供输入 Excel 文件路径。", call. = FALSE)
}

if (!file.exists(opt$input)) {
  stop("输入 Excel 文件不存在：", opt$input, call. = FALSE)
}

# ===================== 创建输出目录 =====================
if (!dir.exists(opt$outdir)) {
  if (!dir.create(opt$outdir, recursive = TRUE)) {
    stop("无法创建输出目录：", opt$outdir, call. = FALSE)
  }
}

# ===================== 获取所有sheet名称 =====================
sheet_names <- excel_sheets(opt$input)

cat("输出目录：", normalizePath(opt$outdir), "\n")
cat("检测到以下 Sheet：\n")
print(sheet_names)

# ===================== 逐个sheet处理 =====================
for (sheet_index in seq_along(sheet_names)) {
  sheet <- sheet_names[[sheet_index]]
  cat(sprintf(
    "[%d/%d] 正在处理 Sheet：%s\n",
    sheet_index, length(sheet_names), sheet
  ))
  
  # 读取当前sheet
  df <- read_excel(opt$input, sheet = sheet, na = c("NA", "", "NaN"))
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  
  # 检查列数
  if (ncol(df) < 2) {
    warning(paste("数据表", sheet, "少于2列，跳过。"))
    next
  }
  
  # 从第二列开始逐个表型输出
  for (j in 2:ncol(df)) {
    phenotype_id <- colnames(df)[j]
    phenotype_value <- df[[j]]
    
    out_df <- data.frame(
#      ID1 = sample_id,
#      ID2 = sample_id,
      Trait = phenotype_value,
      stringsAsFactors = FALSE
    )
    
    # 输出文件名：Sheet名称-性状名.GEMMA.list
    out_file <- file.path(
      opt$outdir,
      paste0(sheet, "-", phenotype_id, ".GEMMA.list")
    )
    
    write.table(
      out_df,
      file = out_file,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      na = "NA"
    )
    
    cat("  已输出：", out_file, "\n")
  }
}

cat("全部处理完成！\n")
