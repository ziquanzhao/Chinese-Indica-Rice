#!/usr/bin/env bash

set -Eeuo pipefail

PROGRAM_NAME=${0##*/}
VERSION="2.2.0"

# 仅在交互式终端中为流程步骤添加青色；重定向日志时不写入 ANSI 控制符。
if [[ -t 1 && -z ${NO_COLOR:-} ]]; then
    LOG_CYAN=$'\033[36m'
    LOG_RESET=$'\033[0m'
else
    LOG_CYAN=""
    LOG_RESET=""
fi

show_help() {
    cat <<'EOF'

GWAS-GEMMA.sh：使用 PLINK2、BCFtools、GCTA 和 GEMMA 完成 GWAS 分析

功能:
  1. 将已完成 QC/MISS/MAF 过滤的总 VCF（SNP + INDEL）转换为 BED。
  2. 从总 VCF 提取 SNP，进行 LD 剪枝，并基于剪枝位点计算 GCTA PCA 以协变量 cov 和 GEMMA 亲缘关系矩阵。
  3. 读取多 Sheet Excel 表型文件，审核表型值和样本 ID，并按 VCF 样本顺序生成 *.GEMMA.list。
  4. 遍历本次生成的所有 *.GEMMA.list 文件，生成 GEMMA_result/gemma.sh。
  5. 默认立即执行 gemma.sh；指定 --dry-run 时只生成、不执行最终 GWAS。

用法:
  bash GWAS-GEMMA.sh --input-vcf FILE --input-phenotype FILE [选项]

必需参数:
  --input-vcf FILE
      输入总 VCF 文件。支持 .vcf 和 .vcf.gz（后缀不区分大小写）。输入文件应已完成 QC、缺失率和 MAF 过滤。

  --input-phenotype FILE
      输入包含表型数据的 Excel 文件。支持 .xlsx 和 .xls（后缀不区分大小写）。每个 Sheet 的第一列必须是样本 ID，第二列起每列代表一个表型，第一行是性状名称。

可选参数:
  --output-dir DIR
      输出目录。默认：当前目录。

  --autosome-num INT
      GCTA 的常染色体数目。默认自动统计输入 VCF 数据记录中 CHROM 列。去重后的唯一值个数。手动提供该参数时使用指定值覆盖自动推断值。

  --indep-pairwise "WINDOW STEP R2"
      PLINK2 LD 剪枝参数。默认："50 10 0.2"。

  --pc-number INT
      写入协变量文件的 PC 数目。默认：3。

  --gk {1|2}
      GEMMA 亲缘关系矩阵算法：1=centered，2=standardized。默认：2。

  --lmm {1|2|3|4}
      GEMMA LMM 检验：1=Wald，2=LRT，3=Score，4=全部。默认：1。

  --miss FLOAT
      GEMMA 位点缺失率阈值，范围 [0,1]。默认：0.1。

  --maf FLOAT
      GEMMA 最小等位基因频率阈值，范围 [0,0.5]。默认：0.01。

  --r2 FLOAT
      GEMMA 相邻位点 LD 阈值，范围 [0,1]。默认：1。

  --parallel INT
      gemma.sh 同时并行运行的性状数。默认：1。该参数不改变单个 GEMMA 任务自身的线程数。

  --dry-run
      完成基因型预处理、PCA、亲缘关系矩阵计算并生成 gemma.sh，但不执行最终逐性状 GWAS。之后可手动运行：bash OUTPUT_DIR/GEMMA_result/gemma.sh

  -h, --help
      显示本帮助信息。

  --version
      显示版本号。

输出目录:
  OUTPUT_DIR/GEMMA_genotype/:	总VCF的BED格式、 OnlySNP的VCF、 LD剪枝的BED/VCF、 GCTA的GRM/PCA、 PCA协变量文件、VCF样本顺序与CHROM唯一值清单，以及GEMMA亲缘关系矩阵。
  OUTPUT_DIR/GEMMA_phenotype/:	从 Excel 各 Sheet 生成并按 VCF 样本顺序排列的 *.GEMMA.list。文件名格式为：Sheet名称-表型列名.GEMMA.list。
  OUTPUT_DIR/GEMMA_result/:	所有性状 GWAS 结果和 gemma.sh。结果前缀格式为：Sheet名称-表型列名.PCA<pc-number>.gk<gk>。

Excel 表型要求与审核:
  * 每个 Sheet 第一列为样本 ID；ID 不得为空或重复。
  * 第二列起为表型；列名不得为空或重复，非空值必须是有限数值。
  * 空单元格以及 NA/NaN 文本会统一输出为 NA。
  * 每个 Sheet 的样本 ID 必须与 VCF 样本 ID 一一对应。
  * 若 ID 集合一致但顺序不同，脚本会自动按 VCF 样本顺序重排。
  * 若 Excel 存在多余 ID 或缺少 VCF ID，脚本会列出不匹配 ID 并退出。

依赖:
  gemma、plink2、bcftools、gcta（也兼容命令名 gcta64）、Rscript、R 包 readxl

使用示例:
  # 使用默认参数并立即执行
  conda activate GWAS_GEMMA
  bash GWAS-GEMMA.sh --input-vcf Xian.missing0.9.maf0.01.Biallelic.OnlyGT.vcf.gz --input-phenotype Phenotype.xlsx

  # 手动指定常染色体数、3个PC、标准亲缘关系矩阵、4个并行性状、仅生成运行脚本但不执行
  bash GWAS-GEMMA.sh --input-vcf Xian.missing0.9.maf0.01.Biallelic.OnlyGT.vcf.gz --input-phenotype Phenotype.xlsx --autosome-num 12 --pc-number 3 --gk 2 --parallel 4 --dry-run
EOF
}

log() {
    printf '%s[%s] %s%s\n' \
        "$LOG_CYAN" "$(date '+%F %T')" "$*" "$LOG_RESET"
}

warn() {
    printf '警告: %s\n' "$*" >&2
}

die() {
    printf '错误: %s\n' "$*" >&2
    exit 1
}

run_command() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
    "$@"
}

require_option_value() {
    local option=$1
    local value=${2-}
    [[ -n "$value" && "$value" != --* ]] ||
        die "参数 $option 缺少取值。使用 --help 查看帮助。"
}

is_positive_integer() {
    [[ $1 =~ ^[1-9][0-9]*$ ]]
}

is_decimal() {
    [[ $1 =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]
}

decimal_in_range() {
    local value=$1
    local lower=$2
    local upper=$3
    is_decimal "$value" &&
        awk -v value="$value" -v lower="$lower" -v upper="$upper" \
            'BEGIN { exit !(value >= lower && value <= upper) }'
}

require_nonempty_file() {
    [[ -s $1 ]] || die "预期输出文件不存在或为空：$1"
}

check_second_column_id_order() {
    local reference_file=$1
    local table_file=$2
    local description=$3

    if ! paste "$reference_file" <(awk '{ print $2 }' "$table_file") |
        awk -v description="$description" '
            NF != 2 || $1 != $2 {
                printf "%s 第 %d 个样本不一致：VCF=%s，当前文件=%s\n",
                       description, NR, $1, $2 > "/dev/stderr"
                exit 1
            }
        '; then
        die "$description 的样本数或样本顺序与 VCF 不一致。"
    fi
}

generate_phenotype_files() {
    local excel_file=$1
    local sample_file=$2
    local phenotype_dir=$3
    local manifest_file=$4

    Rscript --vanilla - \
        "$excel_file" \
        "$sample_file" \
        "$phenotype_dir" \
        "$manifest_file" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
excel_file <- args[[1]]
sample_file <- args[[2]]
phenotype_dir <- args[[3]]
manifest_file <- args[[4]]

suppressPackageStartupMessages(library(readxl))

abort <- function(...) {
    message("错误: ", paste0(..., collapse = ""))
    quit(save = "no", status = 1)
}

show_ids <- function(ids, limit = 50L) {
    shown <- head(ids, limit)
    text <- paste(shown, collapse = ", ")
    if (length(ids) > limit) {
        text <- paste0(text, sprintf(" ...（另有 %d 个）", length(ids) - limit))
    }
    text
}

invalid_filename <- function(value) {
    grepl("[/\\\\\r\n]", value) || !nzchar(trimws(value))
}

vcf_ids <- readLines(sample_file, warn = FALSE, encoding = "UTF-8")
if (length(vcf_ids) == 0L) {
    abort("VCF 样本清单为空：", sample_file)
}
if (anyNA(vcf_ids) || any(!nzchar(vcf_ids))) {
    abort("VCF 样本清单中存在空样本 ID。")
}
if (anyDuplicated(vcf_ids)) {
    duplicated_ids <- unique(vcf_ids[duplicated(vcf_ids)])
    abort(
        "VCF 中存在重复样本 ID（", length(duplicated_ids), " 个）：",
        show_ids(duplicated_ids)
    )
}

sheet_names <- tryCatch(
    excel_sheets(excel_file),
    error = function(error) abort("无法读取 Excel 工作簿：", conditionMessage(error))
)
if (length(sheet_names) == 0L) {
    abort("Excel 工作簿不包含任何 Sheet：", excel_file)
}

phenotype_outputs <- list()
output_basenames <- character()
reordered_sheets <- character()
total_missing <- 0L

for (sheet_index in seq_along(sheet_names)) {
    sheet_name <- sheet_names[[sheet_index]]
    cat(sprintf(
        "[表型审核 %d/%d] Sheet：%s\n",
        sheet_index, length(sheet_names), sheet_name
    ))

    data <- tryCatch(
        read_excel(
            excel_file,
            sheet = sheet_name,
            col_types = "text",
            trim_ws = TRUE,
            na = c("", "NA", "NaN", "nan"),
            .name_repair = "minimal"
        ),
        error = function(error) {
            abort(
                "读取 Sheet [", sheet_name, "] 失败：",
                conditionMessage(error)
            )
        }
    )
    data <- as.data.frame(
        data,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    if (ncol(data) < 2L) {
        abort(
            "Sheet [", sheet_name,
            "] 至少需要两列：第一列样本 ID，第二列起为表型。"
        )
    }

    column_names <- names(data)
    if (anyNA(column_names) || any(!nzchar(trimws(column_names)))) {
        bad_columns <- which(is.na(column_names) | !nzchar(trimws(column_names)))
        abort(
            "Sheet [", sheet_name, "] 存在空列名，列号：",
            paste(bad_columns, collapse = ", ")
        )
    }

    trait_names <- trimws(column_names[-1L])
    if (anyDuplicated(trait_names)) {
        duplicated_traits <- unique(trait_names[duplicated(trait_names)])
        abort(
            "Sheet [", sheet_name, "] 存在重复表型列名：",
            show_ids(duplicated_traits)
        )
    }
    if (invalid_filename(sheet_name)) {
        abort(
            "Sheet 名称不能包含 /、\\ 或换行符，当前名称：[",
            sheet_name, "]"
        )
    }
    invalid_traits <- trait_names[vapply(
        trait_names,
        invalid_filename,
        logical(1)
    )]
    if (length(invalid_traits) > 0L) {
        abort(
            "Sheet [", sheet_name,
            "] 的以下表型列名为空或包含 /、\\、换行符：",
            show_ids(invalid_traits)
        )
    }

    sample_ids <- trimws(as.character(data[[1L]]))
    blank_sample_rows <- which(is.na(sample_ids) | !nzchar(sample_ids))
    if (length(blank_sample_rows) > 0L) {
        abort(
            "Sheet [", sheet_name, "] 存在空样本 ID，Excel 行号：",
            paste(blank_sample_rows + 1L, collapse = ", ")
        )
    }
    if (anyDuplicated(sample_ids)) {
        duplicated_ids <- unique(sample_ids[duplicated(sample_ids)])
        abort(
            "Sheet [", sheet_name, "] 存在重复样本 ID（",
            length(duplicated_ids), " 个）：", show_ids(duplicated_ids)
        )
    }

    extra_ids <- setdiff(sample_ids, vcf_ids)
    missing_ids <- setdiff(vcf_ids, sample_ids)
    if (length(extra_ids) > 0L || length(missing_ids) > 0L) {
        mismatch_message <- paste0(
            "Sheet [", sheet_name,
            "] 的样本 ID 无法与 VCF 样本一一对应。"
        )
        if (length(extra_ids) > 0L) {
            mismatch_message <- paste0(
                mismatch_message,
                "\n  Excel 中存在、VCF 中不存在（", length(extra_ids),
                " 个）：", show_ids(extra_ids)
            )
        }
        if (length(missing_ids) > 0L) {
            mismatch_message <- paste0(
                mismatch_message,
                "\n  VCF 中存在、Excel 中不存在（", length(missing_ids),
                " 个）：", show_ids(missing_ids)
            )
        }
        abort(mismatch_message)
    }

    order_index <- match(vcf_ids, sample_ids)
    if (!identical(sample_ids, vcf_ids)) {
        reordered_sheets <- c(reordered_sheets, sheet_name)
        cat("  样本集合一致但顺序不同：将自动按照 VCF 样本顺序重排。\n")
    } else {
        cat("  样本 ID 和顺序与 VCF 完全一致。\n")
    }

    for (trait_index in seq_along(trait_names)) {
        source_column <- trait_index + 1L
        trait_name <- trait_names[[trait_index]]
        value_text <- trimws(as.character(data[[source_column]]))
        missing_value <- is.na(value_text) |
            !nzchar(value_text) |
            tolower(value_text) %in% c("na", "nan")
        numeric_value <- suppressWarnings(as.numeric(value_text))
        invalid_value <- !missing_value &
            (is.na(numeric_value) | !is.finite(numeric_value))

        if (any(invalid_value)) {
            invalid_rows <- which(invalid_value)
            detail_rows <- head(invalid_rows, 20L)
            details <- paste(
                sprintf(
                    "Excel第%d行=%s",
                    detail_rows + 1L,
                    shQuote(value_text[detail_rows])
                ),
                collapse = "；"
            )
            if (length(invalid_rows) > length(detail_rows)) {
                details <- paste0(
                    details,
                    sprintf(
                        "；另有%d个非数字值",
                        length(invalid_rows) - length(detail_rows)
                    )
                )
            }
            abort(
                "Sheet [", sheet_name, "] 表型 [", trait_name,
                "] 存在非数字或非有限值：", details
            )
        }

        numeric_value[missing_value] <- NA_real_
        numeric_value <- numeric_value[order_index]
        missing_count <- sum(is.na(numeric_value))
        total_missing <- total_missing + missing_count

        output_basename <- paste0(
            sheet_name, "-", trait_name, ".GEMMA.list"
        )
        if (output_basename %in% output_basenames) {
            abort(
                "不同 Sheet/表型生成了重复输出文件名：",
                output_basename
            )
        }
        output_basenames <- c(output_basenames, output_basename)
        output_path <- file.path(phenotype_dir, output_basename)
        phenotype_outputs[[length(phenotype_outputs) + 1L]] <- list(
            path = output_path,
            values = numeric_value
        )

        if (missing_count > 0L) {
            cat(sprintf(
                "  表型 [%s]：%d 个空值将写为 NA。\n",
                trait_name, missing_count
            ))
        }
    }
}

if (length(phenotype_outputs) == 0L) {
    abort("Excel 工作簿中没有可输出的表型列。")
}

dir.create(phenotype_dir, recursive = TRUE, showWarnings = FALSE)
output_paths <- character(length(phenotype_outputs))
for (output_index in seq_along(phenotype_outputs)) {
    output <- phenotype_outputs[[output_index]]
    write.table(
        data.frame(Trait = output$values),
        file = output$path,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        na = "NA"
    )
    output_paths[[output_index]] <- output$path
}
writeLines(output_paths, manifest_file, useBytes = TRUE)

cat(sprintf(
    "表型审核完成：%d 个 Sheet，生成 %d 个表型文件，NA 共 %d 个。\n",
    length(sheet_names), length(output_paths), total_missing
))
if (length(reordered_sheets) > 0L) {
    cat(
        "已自动重排样本顺序的 Sheet：",
        paste(reordered_sheets, collapse = ", "),
        "\n",
        sep = ""
    )
}
RSCRIPT
}

INPUT_VCF=""
INPUT_PHENOTYPE=""
OUTPUT_DIR="."
INDEP_PAIRWISE="50 10 0.2"
AUTOSOME_NUM=""
PC_NUMBER=3
GK=2
LMM=1
MISS=0.1
MAF=0.01
GEMMA_R2=1
PARALLEL=1
DRY_RUN=false

while (($# > 0)); do
    case $1 in
        --input-vcf)
            require_option_value "$1" "${2-}"
            INPUT_VCF=$2
            shift 2
            ;;
        --output-dir)
            require_option_value "$1" "${2-}"
            OUTPUT_DIR=$2
            shift 2
            ;;
        --input-phenotype)
            require_option_value "$1" "${2-}"
            INPUT_PHENOTYPE=$2
            shift 2
            ;;
        --indep-pairwise)
            require_option_value "$1" "${2-}"
            INDEP_PAIRWISE=$2
            shift 2
            ;;
        --autosome-num)
            require_option_value "$1" "${2-}"
            AUTOSOME_NUM=$2
            shift 2
            ;;
        --pc-number)
            require_option_value "$1" "${2-}"
            PC_NUMBER=$2
            shift 2
            ;;
        --gk)
            require_option_value "$1" "${2-}"
            GK=$2
            shift 2
            ;;
        --lmm)
            require_option_value "$1" "${2-}"
            LMM=$2
            shift 2
            ;;
        --miss)
            require_option_value "$1" "${2-}"
            MISS=$2
            shift 2
            ;;
        --maf)
            require_option_value "$1" "${2-}"
            MAF=$2
            shift 2
            ;;
        --r2)
            require_option_value "$1" "${2-}"
            GEMMA_R2=$2
            shift 2
            ;;
        --parallel)
            require_option_value "$1" "${2-}"
            PARALLEL=$2
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        --version)
            printf '%s %s\n' "$PROGRAM_NAME" "$VERSION"
            exit 0
            ;;
        *)
            die "未知参数：$1。使用 --help 查看帮助。"
            ;;
    esac
done

[[ -n $INPUT_VCF ]] || die "必须提供 --input-vcf。"
[[ -n $INPUT_PHENOTYPE ]] || die "必须提供 --input-phenotype。"
[[ -f $INPUT_VCF ]] || die "输入 VCF 文件不存在：$INPUT_VCF"
[[ -f $INPUT_PHENOTYPE ]] ||
    die "输入 Excel 表型文件不存在：$INPUT_PHENOTYPE"

if [[ -n $AUTOSOME_NUM ]]; then
    is_positive_integer "$AUTOSOME_NUM" ||
        die "--autosome-num 必须是正整数，当前值：$AUTOSOME_NUM"
fi
is_positive_integer "$PC_NUMBER" ||
    die "--pc-number 必须是正整数，当前值：$PC_NUMBER"
is_positive_integer "$PARALLEL" ||
    die "--parallel 必须是正整数，当前值：$PARALLEL"
[[ $GK == 1 || $GK == 2 ]] ||
    die "--gk 只能是 1 或 2，当前值：$GK"
[[ $LMM =~ ^[1-4]$ ]] ||
    die "--lmm 只能是 1、2、3 或 4，当前值：$LMM"
decimal_in_range "$MISS" 0 1 ||
    die "--miss 必须是 [0,1] 范围内的数值，当前值：$MISS"
decimal_in_range "$MAF" 0 0.5 ||
    die "--maf 必须是 [0,0.5] 范围内的数值，当前值：$MAF"
decimal_in_range "$GEMMA_R2" 0 1 ||
    die "--r2 必须是 [0,1] 范围内的数值，当前值：$GEMMA_R2"

read -r -a LD_PARAMS <<<"$INDEP_PAIRWISE"
((${#LD_PARAMS[@]} == 3)) ||
    die "--indep-pairwise 必须包含 3 个值：WINDOW STEP R2。"
is_positive_integer "${LD_PARAMS[0]}" ||
    die "LD WINDOW 必须是正整数，当前值：${LD_PARAMS[0]}"
is_positive_integer "${LD_PARAMS[1]}" ||
    die "LD STEP 必须是正整数，当前值：${LD_PARAMS[1]}"
decimal_in_range "${LD_PARAMS[2]}" 0 1 ||
    die "LD R2 必须是 [0,1] 范围内的数值，当前值：${LD_PARAMS[2]}"

for tool in gemma plink2 bcftools Rscript; do
    command -v "$tool" >/dev/null 2>&1 ||
        die "当前环境中找不到软件：$tool。请安装并加入 PATH。"
done

Rscript --vanilla -e \
    'quit(save="no", status=if (requireNamespace("readxl", quietly=TRUE)) 0 else 1)' ||
    die "当前 R 环境中缺少 readxl 包。请先安装，例如：conda install -c conda-forge r-readxl"

if command -v gcta >/dev/null 2>&1; then
    GCTA_COMMAND=$(command -v gcta)
elif command -v gcta64 >/dev/null 2>&1; then
    GCTA_COMMAND=$(command -v gcta64)
    warn "未找到 gcta，改用兼容命令：gcta64"
else
    die "当前环境中找不到软件：gcta 或 gcta64。请安装并加入 PATH。"
fi

mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
INPUT_VCF=$(realpath "$INPUT_VCF")
INPUT_PHENOTYPE=$(realpath "$INPUT_PHENOTYPE")

PHENOTYPE_FILENAME=${INPUT_PHENOTYPE##*/}
shopt -s nocasematch
if [[ ! $PHENOTYPE_FILENAME =~ \.(xlsx|xls)$ ]]; then
    die "--input-phenotype 必须以 .xlsx 或 .xls 结尾：$PHENOTYPE_FILENAME"
fi
shopt -u nocasematch

GENOTYPE_DIR="$OUTPUT_DIR/GEMMA_genotype"
PHENOTYPE_DIR="$OUTPUT_DIR/GEMMA_phenotype"
RESULT_DIR="$OUTPUT_DIR/GEMMA_result"
mkdir -p "$GENOTYPE_DIR" "$PHENOTYPE_DIR" "$RESULT_DIR"

VCF_FILENAME=${INPUT_VCF##*/}
shopt -s nocasematch
if [[ $VCF_FILENAME =~ \.vcf\.gz$ ]]; then
    VCF_PREFIX=${VCF_FILENAME:0:${#VCF_FILENAME}-7}
elif [[ $VCF_FILENAME =~ \.vcf$ ]]; then
    VCF_PREFIX=${VCF_FILENAME:0:${#VCF_FILENAME}-4}
else
    die "--input-vcf 必须以 .vcf 或 .vcf.gz 结尾：$VCF_FILENAME"
fi
shopt -u nocasematch
[[ -n $VCF_PREFIX ]] || die "无法从输入文件名获得输出前缀。"

FULL_BFILE="$GENOTYPE_DIR/$VCF_PREFIX"
ONLY_SNP_VCF="$GENOTYPE_DIR/$VCF_PREFIX.OnlySNP.vcf.gz"
LD_PREFIX="$GENOTYPE_DIR/$VCF_PREFIX.OnlySNP.LDfilter"
GCTA_PREFIX="$LD_PREFIX.GCTA"
PCA_PREFIX="$LD_PREFIX.PCA"
COV_FILE="$LD_PREFIX.PCA${PC_NUMBER}.cov"
SAMPLE_LIST="$GENOTYPE_DIR/$VCF_PREFIX.samples.txt"
CHROM_LIST="$GENOTYPE_DIR/$VCF_PREFIX.chromosomes.txt"
KINSHIP_OUTPUT_NAME="$VCF_PREFIX.OnlySNP.LDfilter.GEMMA.gk${GK}"
GEMMA_SCRIPT="$RESULT_DIR/gemma.sh"
PHENOTYPE_MANIFEST="$PHENOTYPE_DIR/.GWAS-GEMMA.current-files.txt"

log "检测软件依赖"
printf '  gemma   : %s\n' "$(command -v gemma)"
printf '  plink2  : %s\n' "$(command -v plink2)"
printf '  bcftools: %s\n' "$(command -v bcftools)"
printf '  gcta    : %s\n' "$GCTA_COMMAND"
printf '  Rscript : %s\n' "$(command -v Rscript)"
printf '  R readxl: available\n'

log "统计 VCF CHROM 唯一值并推断常染色体数"
if ! bcftools query -f '%CHROM\n' "$INPUT_VCF" |
    awk 'NF > 0 && !seen[$0]++ { print $0 }' >"$CHROM_LIST"; then
    die "无法读取 VCF 的 CHROM 列：$INPUT_VCF"
fi
[[ -s $CHROM_LIST ]] ||
    die "VCF 中没有任何变异记录，无法推断 --autosome-num：$INPUT_VCF"
INFERRED_CHROM_COUNT=$(awk 'END { print NR }' "$CHROM_LIST")
CHROM_PREVIEW=$(awk '
    NR <= 20 {
        values = values (NR == 1 ? "" : ",") $0
    }
    END {
        if (NR > 20) {
            values = values sprintf(",...（另有%d个）", NR - 20)
        }
        print values
    }
' "$CHROM_LIST")

if [[ -z $AUTOSOME_NUM ]]; then
    AUTOSOME_NUM=$INFERRED_CHROM_COUNT
    AUTOSOME_NUM_SOURCE="VCF CHROM 自动推断"
else
    AUTOSOME_NUM_SOURCE="用户指定"
    if [[ $AUTOSOME_NUM != "$INFERRED_CHROM_COUNT" ]]; then
        warn "用户指定的 --autosome-num ($AUTOSOME_NUM) 与 VCF CHROM 唯一值数量 ($INFERRED_CHROM_COUNT) 不一致，将使用用户指定值。"
    fi
fi
printf '  CHROM 唯一值数 : %s\n' "$INFERRED_CHROM_COUNT"
printf '  CHROM 唯一值   : %s\n' "$CHROM_PREVIEW"
printf '  CHROM 清单     : %s\n' "$CHROM_LIST"

log "读取 VCF 样本顺序"
bcftools query -l "$INPUT_VCF" >"$SAMPLE_LIST"
require_nonempty_file "$SAMPLE_LIST"
SAMPLE_COUNT=$(awk 'END { print NR }' "$SAMPLE_LIST")
((PC_NUMBER < SAMPLE_COUNT)) ||
    die "--pc-number ($PC_NUMBER) 必须小于样本数 ($SAMPLE_COUNT)。"
PCA_COMPONENTS=10
((PC_NUMBER > PCA_COMPONENTS)) && PCA_COMPONENTS=$PC_NUMBER
((PCA_COMPONENTS < SAMPLE_COUNT)) || PCA_COMPONENTS=$((SAMPLE_COUNT - 1))

log "审核 Excel 表型并生成 GEMMA 表型文件"
if ! generate_phenotype_files \
    "$INPUT_PHENOTYPE" \
    "$SAMPLE_LIST" \
    "$PHENOTYPE_DIR" \
    "$PHENOTYPE_MANIFEST"; then
    die "Excel 表型审核或转换失败，未开始基因型预处理。"
fi
require_nonempty_file "$PHENOTYPE_MANIFEST"
mapfile -t PHENOTYPE_FILES <"$PHENOTYPE_MANIFEST"
((${#PHENOTYPE_FILES[@]} > 0)) ||
    die "没有从 Excel 生成任何 *.GEMMA.list 文件。"

log "复核 ${#PHENOTYPE_FILES[@]} 个 GEMMA 表型文件"
for phenotype_file in "${PHENOTYPE_FILES[@]}"; do
    [[ -f $phenotype_file ]] ||
        die "表型文件清单中的文件不存在：$phenotype_file"
    [[ ${phenotype_file##*/} == *.GEMMA.list ]] ||
        die "生成的表型文件后缀异常：$phenotype_file"

    if ! awk '
        function valid(value) {
            return value == "NA" ||
                   value ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
        }
        NF != 1 {
            printf "第 %d 行不是单列数据\n", NR > "/dev/stderr"
            exit 1
        }
        !valid($1) {
            printf "第 %d 行不是有效数值或 NA：%s\n", NR, $1 > "/dev/stderr"
            exit 1
        }
    ' "$phenotype_file"; then
        die "表型文件格式错误：$phenotype_file"
    fi

    phenotype_rows=$(awk 'END { print NR }' "$phenotype_file")
    [[ $phenotype_rows == "$SAMPLE_COUNT" ]] ||
        die "表型行数与 VCF 样本数不一致：$phenotype_file（表型 $phenotype_rows 行，VCF $SAMPLE_COUNT 个样本）"
done

printf '\n'
log "运行参数"
printf '  输入 VCF       : %s\n' "$INPUT_VCF"
printf '  输入 Excel     : %s\n' "$INPUT_PHENOTYPE"
printf '  输出目录       : %s\n' "$OUTPUT_DIR"
printf '  样本数         : %s\n' "$SAMPLE_COUNT"
printf '  表型文件数     : %s\n' "${#PHENOTYPE_FILES[@]}"
printf '  LD 参数        : %s %s %s\n' "${LD_PARAMS[@]}"
printf '  VCF CHROM 数   : %s\n' "$INFERRED_CHROM_COUNT"
printf '  常染色体数     : %s（%s）\n' "$AUTOSOME_NUM" "$AUTOSOME_NUM_SOURCE"
printf '  GCTA PCA 数    : %s\n' "$PCA_COMPONENTS"
printf '  协变量 PC 数   : %s\n' "$PC_NUMBER"
printf '  gk / lmm       : %s / %s\n' "$GK" "$LMM"
printf '  miss/maf/r2    : %s / %s / %s\n' "$MISS" "$MAF" "$GEMMA_R2"
printf '  性状并行数     : %s\n' "$PARALLEL"
printf '  dry-run        : %s\n\n' "$DRY_RUN"

log "[1/8] 将总 VCF 转换为 BED"
run_command plink2 \
    --vcf "$INPUT_VCF" \
    --allow-extra-chr \
    --make-bed \
    --out "$FULL_BFILE"
for extension in bed bim fam; do
    require_nonempty_file "$FULL_BFILE.$extension"
done
check_second_column_id_order "$SAMPLE_LIST" "$FULL_BFILE.fam" "总 VCF BED/FAM"

log "[2/8] 从总 VCF 提取 SNP"
run_command bcftools view \
    -v snps \
    --threads 3 \
    "$INPUT_VCF" \
    -Oz \
    -o "$ONLY_SNP_VCF"
require_nonempty_file "$ONLY_SNP_VCF"
run_command bcftools index -f -t "$ONLY_SNP_VCF"

log "[3/8] 对 SNP 执行 LD 剪枝"
TMP_PREFIX="$LD_PREFIX.tmp"
run_command plink2 \
    --vcf "$ONLY_SNP_VCF" \
    --make-bed \
    --allow-extra-chr \
    --out "$TMP_PREFIX"
run_command plink2 \
    --bfile "$TMP_PREFIX" \
    --indep-pairwise "${LD_PARAMS[@]}" \
    --allow-extra-chr \
    --out "$LD_PREFIX"
require_nonempty_file "$LD_PREFIX.prune.in"
run_command plink2 \
    --bfile "$TMP_PREFIX" \
    --extract "$LD_PREFIX.prune.in" \
    --make-bed \
    --allow-extra-chr \
    --out "$LD_PREFIX"
run_command plink2 \
    --bfile "$LD_PREFIX" \
    --recode vcf bgz \
    --allow-extra-chr \
    --out "$LD_PREFIX"
run_command bcftools index -f -t "$LD_PREFIX.vcf.gz"
for extension in bed bim fam vcf.gz vcf.gz.tbi; do
    require_nonempty_file "$LD_PREFIX.$extension"
done
check_second_column_id_order "$SAMPLE_LIST" "$LD_PREFIX.fam" "LD 剪枝 BED/FAM"
rm -f \
    "$TMP_PREFIX.bed" \
    "$TMP_PREFIX.bim" \
    "$TMP_PREFIX.fam" \
    "$TMP_PREFIX.log" \
    "$LD_PREFIX.prune.in" \
    "$LD_PREFIX.prune.out"

log "[4/8] 使用 GCTA 计算 GRM"
run_command "$GCTA_COMMAND" \
    --bfile "$LD_PREFIX" \
    --make-grm \
    --autosome-num "$AUTOSOME_NUM" \
    --out "$GCTA_PREFIX"
for extension in grm.bin grm.N.bin grm.id; do
    require_nonempty_file "$GCTA_PREFIX.$extension"
done

log "[5/8] 使用 GCTA 计算 PCA 并生成协变量文件"
run_command "$GCTA_COMMAND" \
    --grm "$GCTA_PREFIX" \
    --pca "$PCA_COMPONENTS" \
    --out "$PCA_PREFIX"
require_nonempty_file "$PCA_PREFIX.eigenvec"
check_second_column_id_order "$SAMPLE_LIST" "$PCA_PREFIX.eigenvec" "GCTA PCA"

awk -v pc_number="$PC_NUMBER" '
    BEGIN {
        FS = "[[:space:]]+"
        OFS = "\t"
    }
    NF < pc_number + 2 {
        printf "PCA 文件第 %d 行只有 %d 列，无法提取 %d 个 PC\n",
               NR, NF, pc_number > "/dev/stderr"
        exit 1
    }
    {
        printf "1"
        for (column = 3; column <= pc_number + 2; column++) {
            printf "%s%s", OFS, $column
        }
        printf "\n"
    }
' "$PCA_PREFIX.eigenvec" >"$COV_FILE"
require_nonempty_file "$COV_FILE"
cov_rows=$(awk 'END { print NR }' "$COV_FILE")
[[ $cov_rows == "$SAMPLE_COUNT" ]] ||
    die "协变量文件行数异常：预期 $SAMPLE_COUNT，实际 $cov_rows。"

log "[6/8] 将 LD 剪枝 FAM 文件第 6 列统一设置为 1"
FAM_TMP=$(mktemp "$LD_PREFIX.fam.tmp.XXXXXX")
awk '
    BEGIN { OFS = "\t" }
    NF < 6 {
        printf "FAM 文件第 %d 行少于 6 列\n", NR > "/dev/stderr"
        exit 1
    }
    {
        $6 = 1
        print
    }
' "$LD_PREFIX.fam" >"$FAM_TMP"
mv "$FAM_TMP" "$LD_PREFIX.fam"

log "[7/8] 使用 GEMMA 计算亲缘关系矩阵"
run_command gemma \
    -bfile "$LD_PREFIX" \
    -gk "$GK" \
    -o "$KINSHIP_OUTPUT_NAME" \
    -outdir "$GENOTYPE_DIR"

if [[ $GK == 1 ]]; then
    KINSHIP_FILE="$GENOTYPE_DIR/$KINSHIP_OUTPUT_NAME.cXX.txt"
else
    KINSHIP_FILE="$GENOTYPE_DIR/$KINSHIP_OUTPUT_NAME.sXX.txt"
fi
require_nonempty_file "$KINSHIP_FILE"

log "[8/8] 生成逐性状 GWAS 脚本"
{
    cat <<EOF
#!/usr/bin/env bash

set -uo pipefail

# 此脚本由 $PROGRAM_NAME 自动生成。
# 手动运行前请激活包含 gemma 的软件环境，例如：
#   conda activate GWAS_GEMMA

PARALLEL=$PARALLEL
failed_jobs=0
pids=()

command -v gemma >/dev/null 2>&1 || {
    printf '错误: 当前环境中找不到 gemma。请先安装或激活相应环境。\n' >&2
    exit 1
}

run_job() {
    local trait=\$1
    local status
    shift
    printf '[开始] %s\n' "\$trait"
    if "\$@"; then
        printf '[完成] %s\n' "\$trait"
    else
        status=\$?
        printf '[失败] %s（退出码：%s）\n' "\$trait" "\$status" >&2
        return "\$status"
    fi
}

wait_batch() {
    local pid
    for pid in "\${pids[@]}"; do
        if ! wait "\$pid"; then
            failed_jobs=\$((failed_jobs + 1))
        fi
    done
    pids=()
}

printf '共 %s 个性状，最大并行数：%s\n' '${#PHENOTYPE_FILES[@]}' "\$PARALLEL"
EOF

    for phenotype_file in "${PHENOTYPE_FILES[@]}"; do
        phenotype_name=${phenotype_file##*/}
        trait_name=${phenotype_name%.GEMMA.list}
        result_prefix="${trait_name}.PCA${PC_NUMBER}.gk${GK}"
        printf '\nrun_job %q gemma' "$trait_name"
        printf ' %q' \
            -bfile "$FULL_BFILE" \
            -k "$KINSHIP_FILE" \
            -lmm "$LMM" \
            -c "$COV_FILE" \
            -p "$phenotype_file" \
            -o "$result_prefix" \
            -outdir "$RESULT_DIR" \
            -miss "$MISS" \
            -maf "$MAF" \
            -r2 "$GEMMA_R2"
        printf ' &\npids+=("$!")\n'
        printf 'if ((${#pids[@]} >= PARALLEL)); then wait_batch; fi\n'
    done

    cat <<'EOF'

((${#pids[@]} == 0)) || wait_batch

if ((failed_jobs > 0)); then
    printf 'GWAS 完成，但有 %d 个性状失败。\n' "$failed_jobs" >&2
    exit 1
fi

printf '全部 GWAS 性状分析完成。\n'
EOF
} >"$GEMMA_SCRIPT"
chmod +x "$GEMMA_SCRIPT"
require_nonempty_file "$GEMMA_SCRIPT"

printf '\n'
log "预处理和命令生成完成"
printf '  GEMMA 亲缘关系矩阵：%s\n' "$KINSHIP_FILE"
printf '  PCA 协变量文件     ：%s\n' "$COV_FILE"
printf '  VCF 样本顺序       ：%s\n' "$SAMPLE_LIST"
printf '  VCF CHROM 唯一值   ：%s\n' "$CHROM_LIST"
printf '  GWAS 运行脚本      ：%s\n' "$GEMMA_SCRIPT"

if [[ $DRY_RUN == true ]]; then
    printf '\n'
    log "--dry-run 已启用，未执行最终逐性状 GWAS。"
    printf '后续手动运行：bash %q\n' "$GEMMA_SCRIPT"
else
    printf '\n'
    log "开始执行逐性状 GWAS"
    bash "$GEMMA_SCRIPT"
    log "全部流程完成"
fi
