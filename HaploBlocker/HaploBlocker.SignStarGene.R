#!/usr/bin/env Rscript
# 将 HaploBlocker 区块列表转换为区块 × 样本的 0/0.5/1 表。
# 默认仅依赖基础 R；Excel 注释需要 readxl。
# 适用于 inbred = FALSE、每个样本保留两条单倍型的结果。

show_help <- function() {
    cat("HaploBlocker 区块样本携带矩阵导出工具\n\n",
        "用法：\n",
        "  Rscript HaploBlocker.SignStarGene.R --input-rds FILE --input-vcf FILE [--output-file FILE] [--star-gene FILE]\n\n",
        "参数：\n",
        "  --input-rds FILE    必填，HaploBlocker 输出的区块列表 RDS 文件。\n",
        "  --input-vcf FILE    必填，生成该结果时使用的原始 VCF，支持 .vcf 和 .vcf.gz。\n",
        "  --output-file FILE 可选，输出文件名，也可以指定完整路径。\n",
        "                     默认取 RDS 文件的基本文件名，去掉末尾 .rds，\n",
        "                     加上 .haplotype.list，保存到当前工作目录。\n",
        "  --star-gene FILE   可选，Excel 基因注释文件；使用此选项需要 readxl 包。\n",
        "                     StarGene 表必须包含 Chr、Start、End、AGIS、Name。\n",
        "                     GFFGene 表必须包含 Chr、Start、End、AGIS。\n",
        "                     表名不区分大小写，兼容 GFFgene；列名区分大小写。\n",
        "  -h, --help         显示本帮助并退出。\n\n",
        "功能与运行逻辑：\n",
        "  1. 从 VCF 表头提取样本，保持原始顺序。\n",
        "  2. 分批扫描全部变异记录的 CHROM 字段，自动确定染色体 ID。\n",
        "     不依据 ##contig 元信息判断；若变异记录包含多条染色体，报错退出。\n",
        "  3. 每个样本仅输出一列，列名为原始 SampleID，保持 VCF 样本顺序。\n",
        "  4. [[6]] 中的编号 2k-1、2k 对应 VCF 第 k 个样本的 Hap1、Hap2。\n",
        "     每个样本的输出值 = 两条单倍型中属于该区块的条数 / 2。\n",
        "     0：两条均不携带；0.5：仅一条携带；1：两条均携带。\n",
        "     这些数值表示区块成员剂量，不是序列相似度；0 不表示 DNA 缺失。\n",
        "  5. Start、End 取自区块 [[2]]$bp、[[3]]$bp，采用 1-based 闭区间。\n",
        "     bp 为 0 或异常时退出；本脚本不会自动修复坐标。\n",
        "  6. Block_ID 后输出 AnyHaplotypeContainRate 和 BothHaplotypeContainRate。\n",
        "     Any：剂量 >= 0.5 的样本数 / 全部样本数；Both：剂量 = 1 的样本数 / 全部样本数。\n",
        "     两列均以百分数表示，保留两位小数，例如 66.67%。\n",
        "  7. 按原始区块顺序写出制表符分隔文件，包含表头，不包含行名。\n\n",
        "可选基因注释逻辑：\n",
        "  仅匹配相同 Chr，且 Block.Start <= Gene.Start、Gene.End <= Block.End。\n",
        "  坐标包含两端；部分重叠不标注，不扩展基因区间，不转换参考组装。\n",
        "  先将完整覆盖的明星基因 Name 写入 StarGene，AGIS 写入 GFFGene。\n",
        "  再补充 GFFGene 表中完整覆盖基因的 AGIS，跳过所有出现在 StarGene 表中的 AGIS。\n",
        "  多项以分号分隔并去重，明星基因优先，各表内保持原行顺序；无命中填 .。\n",
        "  基因表必须与 VCF 使用同一参考组装和染色体命名；必需字段不能为空。\n\n",
        "输出列：\n",
        "  CHROM  Start  End  Block_ID  AnyHaplotypeContainRate  BothHaplotypeContainRate  Sample1  Sample2 ...\n\n",
        "  指定 --star-gene 时，在 Block_ID 后插入 StarGene、GFFGene 两列，\n",
        "  然后是两个比例列，样本矩阵从第九列开始；未指定时仍从第七列开始。\n\n",
        "使用条件：\n",
        "  RDS 须为 inbred = FALSE 生成的区块列表本身。\n",
        "  必须使用相应的原始 VCF，不能更换或重排样本。\n",
        "  成员编号本身不能验证样本身份，脚本无法识别所有输入错配情况。\n",
        "  相对路径基于运行命令时的当前工作目录；指定的输出目录须已存在。\n",
        "  同名输出文件会被覆盖。\n\n",
        "示例：\n",
        "  Rscript HaploBlocker.SignStarGene.R --input-rds Chr01.rds --input-vcf Chr01.vcf.gz --star-gene StarGene.xlsx --output-file Chr01.annotated.haplotype.list\n",
        sep = "")
}

parse_args <- function(args) {
    if (any(args %in% c("-h", "--help"))) {
        show_help()
        return(NULL)
    }
    allowed <- c("--input-rds", "--input-vcf", "--output-file", "--star-gene")
    values <- list()
    i <- 1L
    while (i <= length(args)) {
        token <- args[i]
        key <- sub("=.*$", "", token)
        if (!(key %in% allowed)) {
            stop(sprintf("未知参数：%s；使用 -h 查看帮助。", token))
        }
        if (!is.null(values[[key]])) {
            stop(sprintf("参数重复：%s", key))
        }
        if (grepl("=", token, fixed = TRUE)) {
            value <- substring(token, nchar(key) + 2L)
        } else {
            i <- i + 1L
            if (i > length(args) || startsWith(args[i], "--")) {
                stop(sprintf("参数 %s 缺少文件路径。", key))
            }
            value <- args[i]
        }
        if (!nzchar(value)) {
            stop(sprintf("参数 %s 的文件路径不能为空。", key))
        }
        values[[key]] <- value
        i <- i + 1L
    }
    for (key in c("--input-rds", "--input-vcf")) {
        if (is.null(values[[key]])) {
            stop(sprintf("缺少必填参数 %s；使用 -h 查看帮助。", key))
        }
    }
    if (is.null(values[["--output-file"]])) {
        values[["--output-file"]] <- paste0(
            sub("\\.rds$", "", basename(values[["--input-rds"]]),
                ignore.case = TRUE),
            ".haplotype.list"
        )
    }
    values
}

# 分批扫描变异记录，避免将整个 VCF 基因型矩阵加载到内存。
read_vcf_metadata <- function(vcf_file) {
    con <- gzfile(vcf_file, open = "rt")
    on.exit(close(con))
    repeat {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (!length(line)) stop("没有找到 VCF 的 #CHROM 表头。")
        if (startsWith(line, "#CHROM\t")) {
            fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
            if (length(fields) < 10L) stop("VCF 中没有样本列。")
            samples <- fields[10:length(fields)]
            break
        }
        if (!startsWith(line, "##")) stop("VCF 表头格式异常。")
    }
    if (anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) {
        stop("VCF 样本名称存在缺失、空值或重复。")
    }
    chromosome <- NULL
    n_variants <- 0
    repeat {
        lines <- readLines(con, n = 1000L, warn = FALSE)
        if (!length(lines)) break
        tabs <- regexpr("\t", lines, fixed = TRUE)
        if (any(tabs < 2L) || any(startsWith(lines, "#"))) {
            stop("VCF 变异记录格式异常，CHROM 字段缺失或存在意外表头。")
        }
        chromosomes <- unique(substr(lines, 1L, tabs - 1L))
        if (any(chromosomes == ".")) stop("VCF 变异记录的染色体 ID 缺失。")
        observed <- unique(c(chromosome, chromosomes))
        if (length(observed) > 1L) {
            stop(sprintf(
                "警报：VCF 变异记录包含多条染色体（%s），本脚本仅接受单染色体输入，已退出。",
                paste(observed, collapse = ", ")
            ))
        }
        chromosome <- observed
        n_variants <- n_variants + length(lines)
    }
    if (n_variants == 0) stop("VCF 中没有变异记录，无法确定染色体 ID。")
    list(samples = samples, chromosome = chromosome, n_variants = n_variants)
}

# 只在用户请求注释时加载 Excel 依赖；所有必需字段先校验再计算。
read_gene_annotations <- function(path) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("--star-gene 需要 readxl 包，请在运行脚本的 R 环境中安装：install.packages(\"readxl\")")
    }
    sheets <- readxl::excel_sheets(path)
    read_sheet <- function(sheet_name, required) {
        index <- which(tolower(trimws(sheets)) == tolower(sheet_name))
        if (length(index) != 1L) {
            stop(sprintf("Excel 必须有且仅有一个 %s 工作表（表名不区分大小写）。", sheet_name))
        }
        tab <- as.data.frame(readxl::read_excel(
            path, sheet = sheets[index], col_types = "text", .name_repair = "minimal"
        ), check.names = FALSE, stringsAsFactors = FALSE)
        missing <- setdiff(required, names(tab))
        if (length(missing)) {
            stop(sprintf("工作表 %s 缺少必需列：%s", sheets[index], paste(missing, collapse = ", ")))
        }
        if (any(vapply(required, function(key) sum(names(tab) == key) != 1L, logical(1)))) {
            stop(sprintf("工作表 %s 的必需列名存在重复。", sheets[index]))
        }
        tab <- tab[, required, drop = FALSE]
        for (key in required) {
            value <- trimws(tab[[key]])
            bad <- which(is.na(value) | !nzchar(value) | value == ".")
            if (length(bad)) {
                stop(sprintf("工作表 %s 第 %d 条数据的 %s 缺失。", sheets[index], bad[1], key))
            }
            if (any(grepl("[\t\r\n]", value))) {
                stop(sprintf("工作表 %s 的 %s 包含制表符或换行符，无法安全写入 TSV。", sheets[index], key))
            }
            tab[[key]] <- value
        }
        for (key in c("Start", "End")) {
            value <- suppressWarnings(as.numeric(tab[[key]]))
            bad <- which(!is.finite(value) | value < 1 | value != floor(value))
            if (length(bad)) {
                stop(sprintf("工作表 %s 第 %d 条数据的 %s 必须是正整数坐标。", sheets[index], bad[1], key))
            }
            tab[[key]] <- value
        }
        bad <- which(tab$Start > tab$End)
        if (length(bad)) {
            stop(sprintf("工作表 %s 第 %d 条数据的 Start 大于 End。", sheets[index], bad[1]))
        }
        tab
    }
    star <- read_sheet("StarGene", c("Chr", "Start", "End", "AGIS", "Name"))
    gff <- read_sheet("GFFGene", c("Chr", "Start", "End", "AGIS"))
    # 明星基因以 StarGene 表为准，不从 GFFGene 表重复补入。
    gff <- gff[!gff$AGIS %in% star$AGIS, , drop = FALSE]
    list(star = star, gff = gff)
}

annotate_blocks <- function(chr_id, start_bp, end_bp, genes) {
    star <- genes$star[genes$star$Chr == chr_id, , drop = FALSE]
    gff <- genes$gff[genes$gff$Chr == chr_id, , drop = FALSE]
    collapse_hits <- function(values) {
        if (!length(values)) return(".")
        paste(unique(values), collapse = ";")
    }
    star_names <- gff_ids <- rep(".", length(start_bp))
    for (i in seq_along(start_bp)) {
        star_hit <- star$Start >= start_bp[i] & star$End <= end_bp[i]
        gff_hit <- gff$Start >= start_bp[i] & gff$End <= end_bp[i]
        star_names[i] <- collapse_hits(star$Name[star_hit])
        gff_ids[i] <- collapse_hits(c(star$AGIS[star_hit], gff$AGIS[gff_hit]))
    }
    data.frame(StarGene = star_names, GFFGene = gff_ids, stringsAsFactors = FALSE)
}

main <- function() {
    args <- parse_args(commandArgs(trailingOnly = TRUE))
    if (is.null(args)) return(invisible(NULL))
    rds_file <- args[["--input-rds"]]
    vcf_file <- args[["--input-vcf"]]
    output_file <- args[["--output-file"]]
    star_gene_file <- args[["--star-gene"]]

    for (path in c(rds_file, vcf_file, star_gene_file)) {
        if (!file.exists(path) || dir.exists(path)) {
            stop(sprintf("输入文件不存在或不是普通文件：%s", path))
        }
    }
    if (!dir.exists(dirname(output_file))) stop("指定的输出目录不存在。")
    if (dir.exists(output_file)) stop("输出路径指向目录，请指定文件名。")
    input_paths <- normalizePath(c(rds_file, vcf_file, star_gene_file), mustWork = TRUE)
    output_path <- file.path(normalizePath(dirname(output_file), mustWork = TRUE),
                             basename(output_file))
    if (file.exists(output_file)) {
        output_path <- normalizePath(output_file, mustWork = TRUE)
    }
    if (output_path %in% input_paths) stop("输出文件不能与输入文件相同。")

    genes <- NULL
    if (!is.null(star_gene_file)) {
        genes <- read_gene_annotations(star_gene_file)
    }

    cat("正在扫描 VCF 样本顺序及全部变异记录的染色体 ID……\n")
    metadata <- read_vcf_metadata(vcf_file)
    samples <- metadata$samples
    chr_id <- metadata$chromosome
    cat("染色体：", chr_id, "；变异记录数：", metadata$n_variants, "\n")

    blocklist <- readRDS(rds_file)
    if (!is.list(blocklist) || length(blocklist) == 0L) {
        stop("blocklist 必须是非空的区块列表。")
    }
    # ==================== 4. 构建区块 × 样本矩阵 ====================

    n_blocks <- length(blocklist)
    n_samples <- length(samples)
    n_haps <- 2L * n_samples

    # 使用数值矩阵以保存 0.5；样本列严格遵循 VCF 表头，不排序。
    sample_dosage <- matrix(
        0,
        nrow = n_blocks,
        ncol = n_samples,
        dimnames = list(NULL, samples)
    )

    start_bp <- numeric(n_blocks)
    end_bp <- numeric(n_blocks)

    for (i in seq_len(n_blocks)) {
        block <- blocklist[[i]]

        if (!is.list(block) || length(block) < 6L ||
            !is.list(block[[2]]) || !is.list(block[[3]])) {
            stop(sprintf("第 %d 个区块的结构异常。", i))
        }

        # 成员编号对应输入 VCF 中交替排列的两条单倍型
        members <- unname(block[[6]])

        if (!is.numeric(members) ||
            anyNA(members) ||
            any(!is.finite(members)) ||
            any(members < 1 | members > n_haps |
                members != floor(members)) ||
            anyDuplicated(members)) {
            stop(sprintf("第 %d 个区块的单倍型编号异常。", i))
        }

        if (!is.numeric(block[[5]]) || length(block[[5]]) != 1L ||
            is.na(block[[5]]) ||
            block[[5]] != length(members)) {
            stop(sprintf("第 %d 个区块的成员数量与 [[5]] 不一致。", i))
        }

        # 已核对 HaploBlocker 1.7.2 的 data_import()：
        # 先 cbind(haplo1, haplo2)，再按以下表达式重排：
        # c(0, ncol(haplo1)) + rep(1:ncol(haplo1), each = 2)
        # 因此编号 2k-1 和 2k 均属于 VCF 第 k 个样本。
        # blockinfo_calculation()/block_merging() 使用全局列位置记录成员。
        # 上面已验证 members 无重复；每个样本最多计入两条单倍型。
        sample_index <- as.integer(ceiling(members / 2))
        sample_dosage[i, ] <- tabulate(sample_index, nbins = n_samples) / 2
        if (sum(sample_dosage[i, ]) * 2 != length(members)) {
            stop(sprintf("第 %d 个区块的样本剂量总和与成员数量不一致。", i))
        }

        start <- block[[2]]$bp
        end <- block[[3]]$bp

        if (!is.numeric(start) || !is.numeric(end) ||
            length(start) != 1L || length(end) != 1L) {
            stop(sprintf("第 %d 个区块的 bp 坐标格式异常。", i))
        }

        if (!is.finite(start) || !is.finite(end) ||
            start <= 0 || end < start ||
            start != floor(start) || end != floor(end)) {
            stop(sprintf(
                paste0(
                    "第 %d 个区块的 bp 坐标缺失或异常。",
                    "请先根据起止 snp 序号回查原始 VCF，补齐坐标后再导出。"
                ),
                i
            ))
        }

        start_bp[i] <- start
        end_bp[i] <- end
    }

    # ==================== 5. 添加区块信息 ====================

    # 起止坐标采用 VCF 的 1-based 坐标，包含两端
    result <- data.frame(
        CHROM = rep(chr_id, n_blocks),
        Start = start_bp,
        End = end_bp,
        Block_ID = paste0("Block_", seq_len(n_blocks)),
        # 分母为全部样本数；使用尚未格式化的样本剂量计算。
        AnyHaplotypeContainRate = sprintf(
            "%.2f%%", 100 * rowSums(sample_dosage >= 0.5) / n_samples
        ),
        BothHaplotypeContainRate = sprintf(
            "%.2f%%", 100 * rowSums(sample_dosage == 1) / n_samples
        ),
        stringsAsFactors = FALSE
    )

    if (!is.null(genes)) {
        cat("正在按基因完整包含条件进行注释……\n")
        annotation <- annotate_blocks(chr_id, start_bp, end_bp, genes)
        result <- cbind(result[, 1:4, drop = FALSE], annotation,
                        result[, 5:6, drop = FALSE])
    }
    metadata_names <- names(result)
    result <- cbind(
        result,
        as.data.frame(sample_dosage, optional = TRUE)
    )

    # 保留样本原始名称中的连字符等字符
    colnames(result) <- c(metadata_names, samples)

    # ==================== 6. 导出 TSV 文件 ====================

    write_result <- function(result, output_file) {
        # 避免坐标以科学计数法输出；结束后恢复原设置
        old_options <- options(scipen = 999)
        on.exit(options(old_options))

        write.table(
            result,
            file = output_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE
        )
    }

    write_result(result, output_file)

    # ==================== 7. 查看部分结果 ====================

    # 查看前 6 个区块、前 3 个样本
    print(result[
        seq_len(min(6L, n_blocks)),
        seq_len(min(length(metadata_names) + 3L, ncol(result))),
        drop = FALSE
    ])

    cat("区块数：", n_blocks, "\n")
    cat("样本数：", length(samples), "\n")
    cat("样本列数：", n_samples, "\n")
    cat("输出编码：0 = 两条均不携带；0.5 = 仅一条携带；1 = 两条均携带\n")
    cat("总列数：", ncol(result), "\n")
    cat("输出文件：", output_file, "\n")
}

tryCatch(
    main(),
    error = function(e) {
        cat("错误：", conditionMessage(e), "\n", file = stderr(), sep = "")
        quit(save = "no", status = 1L)
    }
)
