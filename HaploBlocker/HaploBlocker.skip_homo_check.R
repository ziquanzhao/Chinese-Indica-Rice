# 仅在当前 R 会话中跳过 VCF 导入时的高度纯合提示检查。
# 不修改安装目录，不改变 inbred 设置，也不填补或删除缺失基因型。
local({
  suppressPackageStartupMessages({
    library(RandomFieldsUtils)
    library(HaploBlocker)
  })

  import_fun <- get("data_import", envir = asNamespace("HaploBlocker"))
  if (isTRUE(attr(import_fun, "homo_check_skipped"))) {
    message("当前会话已经跳过 VCF 高度纯合检查。")
  } else {
    source_lines <- deparse(body(import_fun), width.cutoff = 500L)
    old_check <- "if (mean(haplo1 == haplo2) > 0.95)"
    if (sum(grepl(old_check, source_lines, fixed = TRUE)) != 1L) {
      stop("当前版本的导入函数与预期不符，未修改任何函数。")
    }

    source_lines <- sub(old_check, "if (FALSE)", source_lines, fixed = TRUE)
    body(import_fun) <- parse(text = source_lines)[[1L]]
    attr(import_fun, "homo_check_skipped") <- TRUE
    assignInNamespace("data_import", import_fun, ns = "HaploBlocker")
    message("已跳过 VCF 高度纯合检查；修改仅对当前 R 会话有效。")
  }
})
