#!/usr/bin/env bash

set -euo pipefail

PROGRAM_NAME="$(basename "$0")"

usage() {
    cat <<EOF
用法：
  bash ${PROGRAM_NAME} --input-vcf FILE --sample-group FILE --combination FILE [选项]

功能：
  根据样本分群和群体比较组合，生成并运行Fst、Pi(log2(Pi ratio))和/或XP-CLR
  选择清除分析。使用--dry-run时仅生成分析脚本，不执行分析命令。
  运行前会输出完整参数、亚群体个体数及各分析的比较组合，便于核对。
  终端输出和错误信息会同步追加到输出目录中的SelectiveSweep.log。

必需参数：
  --input-vcf FILE
      输入VCF文件。支持未压缩的.vcf和gzip/bgzip压缩的.vcf.gz/.vcf.bgz文件。

  --sample-group FILE
      样本分群文件。第一行为表头，第一列为样本ID，第二列为亚群体ID。列之间可使用Tab或空白字符分隔。

  --combination FILE
      两两比较组合文件。第一行为表头，第一列为参考组Ref，第二列为检验组Check。XP-CLR中Check对应--samplesA，Ref对应--samplesB；log2(Pi ratio)中Ref为分子，Check为分母。

可选参数：
  --output-dir DIR
      输出目录，默认：运行脚本时的当前目录。脚本将在其中创建 subgroup/，并根据分析类型创建 fst/、pi/ 和/或 xpclr/。

  --analysis-project {fst|pi|xpclr|all}
      运行的分析项目，默认：all。

  --window-size INT
      滑动窗口大小，默认：50000。

  --window-step INT
      滑动窗口步长，默认：10000。

  --parallel-xpclr INT
      XP-CLR 子命令的全局并行数。默认：输入 VCF 中的染色体数。

  --dry-run
      仅创建目录、亚群体样本列表和分析脚本，不执行 fst.sh、pi.sh 或 xpclr.sh。生成的脚本可分别单独运行。

  -h, --help
      显示本帮助文档并退出。

输出文件：
  SelectiveSweep.log
  subgroup/<亚群体ID>.list
  fst/fst.sh
  pi/pi.sh pi/calculate_pi_log2_ratio.R
  xpclr/xpclr.sh

命名规则：
  Ref<参考组>.vs.Check<检验组>.win<窗口>.step<步长>...
  其中 <窗口> 和 <步长> 分别取自 --window-size 和 --window-step

示例：
  # 生成并运行全部分析
  bash ${PROGRAM_NAME} --input-vcf rice.vcf.gz --output-dir ./ --sample-group SampleGroup.list --combination Combination.tsv --analysis-project all --window-size 50000 --window-step 10000 --parallel-xpclr 12

  # 仅生成 Fst 命令，随后手动运行
  bash ${PROGRAM_NAME} --input-vcf rice.vcf.gz --sample-group SampleGroup.list --combination Combination.tsv --analysis-project fst --dry-run

EOF
}

die() {
    printf '错误：%s\n' "$*" >&2
    exit 1
}

info() {
    printf '[%s] %s\n' "$PROGRAM_NAME" "$*"
}

require_value() {
    local option="$1"
    local value="${2-}"
    [[ -n "$value" && "$value" != --* ]] ||
        die "参数 ${option} 需要提供一个值。"
}

absolute_existing_file() {
    local path="$1"
    [[ -f "$path" ]] || die "文件不存在或不是普通文件：${path}"
    realpath "$path"
}

shell_quote() {
    printf '%q' "$1"
}

print_subgroup_summary() {
    local group_id sample_count group_count total_sample_count

    group_count="$(awk 'END { print NR + 0 }' "$GROUP_IDS_FILE")"
    total_sample_count="$(awk 'END { print NR + 0 }' "$NORMALIZED_SAMPLE_GROUP")"
    printf '  亚群体及个体数（共 %s 个亚群体、%s 个个体）：\n' \
        "$group_count" "$total_sample_count"
    while IFS= read -r group_id; do
        sample_count="$(awk 'END { print NR }' "${SUBGROUP_DIR}/${group_id}.list")"
        printf '    - %s：%s 个个体\n' "$group_id" "$sample_count"
    done < "$GROUP_IDS_FILE"
}

print_combination_summary() {
    local heading="${1:-比较组合（第一列为 Ref，第二列为 Check）}"
    local ref_group check_group
    local combination_index=0 combination_count

    combination_count="$(awk 'END { print NR + 0 }' "$NORMALIZED_COMBINATIONS")"
    printf '  %s，共 %s 个组合：\n' "$heading" "$combination_count"
    while IFS=$'\t' read -r ref_group check_group; do
        combination_index=$((combination_index + 1))
        printf '    %d. Ref=%s，Check=%s\n' \
            "$combination_index" "$ref_group" "$check_group"
    done < "$NORMALIZED_COMBINATIONS"
}

print_analysis_summary() {
    local analysis_name="$1"

    case "$analysis_name" in
        fst)
            printf '\n[Fst 执行计划]\n'
            print_subgroup_summary
            print_combination_summary
            ;;
        pi)
            printf '\n[Pi 执行计划]\n'
            print_subgroup_summary
            printf '  将为上述每个亚群体计算 Pi。\n'
            print_combination_summary '以下组合用于计算 log2(Ref Pi / Check Pi)'
            ;;
        xpclr)
            printf '\n[XP-CLR 执行计划]\n'
            print_subgroup_summary
            print_combination_summary
            printf '  XP-CLR方向：samplesA=Check（目标/检验群体），samplesB=Ref（参考群体）\n'
            if [[ "$PARALLEL_XPCLR" == "auto" ]]; then
                printf '  XP-CLR 并行数：自动（运行时取 VCF 染色体数）\n'
            else
                printf '  XP-CLR 并行数：%s\n' "$PARALLEL_XPCLR"
            fi
            ;;
    esac
}

print_run_summary() {
    local dry_run_text="否"

    if [[ "$DRY_RUN" == true ]]; then
        dry_run_text="是"
    fi

    printf '\n========== SelectiveSweep 参数确认 ==========\n'
    printf '  input VCF：%s\n' "$INPUT_VCF"
    printf '  output dir：%s\n' "$OUTPUT_DIR"
    printf '  analysis project：%s\n' "$ANALYSIS_PROJECT"
    printf '  sample-group：%s\n' "$SAMPLE_GROUP_FILE"
    printf '  combination：%s\n' "$COMBINATION_FILE"
    printf '  window size：%s\n' "$WINDOW_SIZE"
    printf '  window step：%s\n' "$WINDOW_STEP"
    if [[ "$PARALLEL_XPCLR" == "auto" ]]; then
        printf '  parallel XP-CLR：自动（VCF 染色体数）\n'
    else
        printf '  parallel XP-CLR：%s\n' "$PARALLEL_XPCLR"
    fi
    printf '  dry-run：%s\n' "$dry_run_text"

    case "$ANALYSIS_PROJECT" in
        fst|pi|xpclr)
            print_analysis_summary "$ANALYSIS_PROJECT"
            ;;
        all)
            print_analysis_summary fst
            print_analysis_summary pi
            print_analysis_summary xpclr
            ;;
    esac
    printf '\n=============================================\n\n'
}

write_script_header() {
    local script_path="$1"
    local analysis_name="$2"
    local group_id sample_count ref_group check_group
    local combination_index=0 group_count total_sample_count combination_count

    group_count="$(awk 'END { print NR + 0 }' "$GROUP_IDS_FILE")"
    total_sample_count="$(awk 'END { print NR + 0 }' "$NORMALIZED_SAMPLE_GROUP")"
    combination_count="$(awk 'END { print NR + 0 }' "$NORMALIZED_COMBINATIONS")"

    {
        printf '%s\n' '#!/usr/bin/env bash'
        printf '\n'
        printf '%s\n' 'set -euo pipefail'
        printf '\n'
        printf '# 自动生成的 %s 分析脚本\n' "$analysis_name"
        printf '# input VCF: %s\n' "$INPUT_VCF"
        printf '# sample-group file: %s\n' "$SAMPLE_GROUP_FILE"
        printf '# combination file: %s\n' "$COMBINATION_FILE"
        printf 'readonly INPUT_VCF=%s\n' "$(shell_quote "$INPUT_VCF")"
        printf 'readonly SAMPLE_GROUP_FILE=%s\n' "$(shell_quote "$SAMPLE_GROUP_FILE")"
        printf 'readonly COMBINATION_FILE=%s\n' "$(shell_quote "$COMBINATION_FILE")"
        printf 'readonly SUBGROUP_DIR=%s\n' "$(shell_quote "$SUBGROUP_DIR")"
        printf 'readonly ANALYSIS_DIR=%s\n' "$(shell_quote "$(dirname "$script_path")")"
        printf 'readonly LOG_FILE=%s\n' "$(shell_quote "$LOG_FILE")"
        printf 'readonly WINDOW_SIZE=%q\n' "$WINDOW_SIZE"
        printf 'readonly WINDOW_STEP=%q\n' "$WINDOW_STEP"
        printf '\n'
        printf '%s\n' 'if [[ "${SELECTIVE_SWEEP_LOG_FILE_ACTIVE:-}" != "$LOG_FILE" ]]; then'
        printf '%s\n' '    command -v tee >/dev/null 2>&1 || { echo "错误：找不到 tee，无法写入日志文件。" >&2; exit 127; }'
        printf '%s\n' '    exec > >(tee -a "$LOG_FILE") 2>&1'
        printf '%s\n' '    export SELECTIVE_SWEEP_LOG_FILE_ACTIVE="$LOG_FILE"'
        printf '%s\n' 'fi'
        printf 'echo %s\n' "$(shell_quote "[$analysis_name] 日志文件：${LOG_FILE}")"
        printf '\n'
        printf 'echo %s\n' "$(shell_quote "========== ${analysis_name} 参数确认 ==========")"
        printf 'echo %s\n' "$(shell_quote "  input VCF：${INPUT_VCF}")"
        printf 'echo %s\n' "$(shell_quote "  analysis dir：$(dirname "$script_path")")"
        printf 'echo %s\n' "$(shell_quote "  sample-group：${SAMPLE_GROUP_FILE}")"
        printf 'echo %s\n' "$(shell_quote "  combination：${COMBINATION_FILE}")"
        printf 'echo %s\n' "$(shell_quote "  window size：${WINDOW_SIZE}")"
        printf 'echo %s\n' "$(shell_quote "  window step：${WINDOW_STEP}")"
        printf 'echo %s\n' \
            "$(shell_quote "  亚群体及个体数（共 ${group_count} 个亚群体、${total_sample_count} 个个体）：")"
        while IFS= read -r group_id; do
            sample_count="$(awk 'END { print NR }' "${SUBGROUP_DIR}/${group_id}.list")"
            printf 'echo %s\n' "$(shell_quote "    - ${group_id}：${sample_count} 个个体")"
        done < "$GROUP_IDS_FILE"
        if [[ "$analysis_name" == "Pi" ]]; then
            printf 'echo %s\n' "$(shell_quote "  上述所有亚群体均会分别计算 Pi。")"
            printf 'echo %s\n' \
                "$(shell_quote "  以下组合计算 log2(Ref Pi / Check Pi)，共 ${combination_count} 个组合：")"
        else
            printf 'echo %s\n' \
                "$(shell_quote "  比较组合（第一列为 Ref，第二列为 Check），共 ${combination_count} 个组合：")"
        fi
        while IFS=$'\t' read -r ref_group check_group; do
            combination_index=$((combination_index + 1))
            printf 'echo %s\n' \
                "$(shell_quote "    ${combination_index}. Ref=${ref_group}，Check=${check_group}")"
        done < "$NORMALIZED_COMBINATIONS"
        if [[ "$analysis_name" == "XP-CLR" ]]; then
            printf 'echo %s\n' \
                "$(shell_quote "  XP-CLR方向：samplesA=Check（目标/检验群体），samplesB=Ref（参考群体）")"
            if [[ "$PARALLEL_XPCLR" == "auto" ]]; then
                printf 'echo %s\n' "$(shell_quote "  XP-CLR 并行数：自动（VCF 染色体数）")"
            else
                printf 'echo %s\n' "$(shell_quote "  XP-CLR 并行数：${PARALLEL_XPCLR}")"
            fi
        fi
        printf 'echo %s\n' "$(shell_quote "=============================================")"
        printf 'echo\n'
        printf '\n'
    } > "$script_path"
}

write_fst_script() {
    local script_path="${FST_DIR}/fst.sh"
    local ref_group check_group output_prefix

    write_script_header "$script_path" "Fst"
    {
        printf '%s\n' 'command -v vcftools >/dev/null 2>&1 || { echo "错误：找不到 vcftools。" >&2; exit 127; }'
        printf '\n'
        printf '%s\n' 'echo "[fst.sh] 开始运行 Fst 分析。"'
    } >> "$script_path"

    while IFS=$'\t' read -r ref_group check_group; do
        output_prefix="${FST_DIR}/Ref${ref_group}.vs.Check${check_group}.win${WINDOW_SIZE}.step${WINDOW_STEP}"
        {
            printf '\n'
            printf 'echo %s\n' "$(shell_quote "[fst.sh] ${ref_group} vs ${check_group}")"
            printf 'vcftools %s %s --weir-fst-pop %s --weir-fst-pop %s --fst-window-size %q --fst-window-step %q --out %s\n' \
                "$VCFTOOLS_INPUT_OPTION" \
                "$(shell_quote "$INPUT_VCF")" \
                "$(shell_quote "${SUBGROUP_DIR}/${ref_group}.list")" \
                "$(shell_quote "${SUBGROUP_DIR}/${check_group}.list")" \
                "$WINDOW_SIZE" \
                "$WINDOW_STEP" \
                "$(shell_quote "$output_prefix")"
        } >> "$script_path"
    done < "$NORMALIZED_COMBINATIONS"

    {
        printf '\n'
        printf '%s\n' 'echo "[fst.sh] Fst 分析完成。"'
    } >> "$script_path"
    chmod +x "$script_path"
}

write_pi_calculator() {
    local calculator_path="${PI_DIR}/calculate_pi_log2_ratio.R"

    cat > "$calculator_path" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7L) {
  stop(
    paste(
      "用法：Rscript calculate_pi_log2_ratio.R",
      "<numerator.pi> <denominator.pi> <numerator_id> <denominator_id>",
      "<window_size> <window_step> <output_file>"
    ),
    call. = FALSE
  )
}

numerator_file <- args[[1]]
denominator_file <- args[[2]]
numerator_id <- args[[3]]
denominator_id <- args[[4]]
window_size <- args[[5]]
window_step <- args[[6]]
output_file <- args[[7]]

read_pi <- function(group_id, input_file) {
  dat <- read.table(
    input_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c("CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "PI")
  missing_columns <- setdiff(required, names(dat))
  if (length(missing_columns) > 0L) {
    stop(
      input_file,
      " 缺少必需列：",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  dat <- dat[, required, drop = FALSE]
  names(dat)[names(dat) == "N_VARIANTS"] <- paste0(group_id, "_N_VARIANTS")
  names(dat)[names(dat) == "PI"] <- paste0(group_id, "_PI")
  dat
}

numerator_data <- read_pi(numerator_id, numerator_file)
denominator_data <- read_pi(denominator_id, denominator_file)
chromosome_order <- unique(numerator_data$CHROM)

merged <- merge(
  numerator_data,
  denominator_data,
  by = c("CHROM", "BIN_START", "BIN_END"),
  all = FALSE,
  sort = FALSE
)

numerator_pi_column <- paste0(numerator_id, "_PI")
denominator_pi_column <- paste0(denominator_id, "_PI")
numerator_pi <- merged[[numerator_pi_column]]
denominator_pi <- merged[[denominator_pi_column]]
valid <- (
  is.finite(numerator_pi) &
  is.finite(denominator_pi) &
  numerator_pi > 0 &
  denominator_pi > 0
)

result <- merged[valid, , drop = FALSE]
result$LOG2_PI_RATIO <- log2(
  result[[numerator_pi_column]] / result[[denominator_pi_column]]
)
result <- result[
  order(
    match(result$CHROM, chromosome_order),
    result$BIN_START,
    result$BIN_END
  ),
  ,
  drop = FALSE
]

write.table(
  result,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  sprintf(
    paste0(
      "%s/%s：共同窗口 %d 个，有效 log2(Pi ratio) 窗口 %d 个，",
      "过滤 %d 个非有限值或非正 Pi 窗口；window=%s，step=%s。"
    ),
    numerator_id,
    denominator_id,
    nrow(merged),
    nrow(result),
    nrow(merged) - nrow(result),
    window_size,
    window_step
  )
)
RSCRIPT
}

write_pi_script() {
    local script_path="${PI_DIR}/pi.sh"
    local group_id ref_group check_group output_prefix ratio_output
    local calculator_path="${PI_DIR}/calculate_pi_log2_ratio.R"

    write_pi_calculator
    write_script_header "$script_path" "Pi"
    {
        printf 'readonly PI_CALCULATOR=%s\n' "$(shell_quote "$calculator_path")"
        printf '%s\n' 'command -v vcftools >/dev/null 2>&1 || { echo "错误：找不到 vcftools。" >&2; exit 127; }'
        printf '%s\n' 'command -v Rscript >/dev/null 2>&1 || { echo "错误：找不到 Rscript。" >&2; exit 127; }'
        printf '\n'
        printf '%s\n' 'echo "[pi.sh] 开始计算各亚群体的 Pi。"'
    } >> "$script_path"

    while IFS= read -r group_id; do
        output_prefix="${PI_DIR}/${group_id}.win${WINDOW_SIZE}.step${WINDOW_STEP}"
        {
            printf '\n'
            printf 'echo %s\n' "$(shell_quote "[pi.sh] 计算 ${group_id} 的 Pi")"
            printf 'vcftools %s %s --window-pi %q --window-pi-step %q --keep %s --out %s\n' \
                "$VCFTOOLS_INPUT_OPTION" \
                "$(shell_quote "$INPUT_VCF")" \
                "$WINDOW_SIZE" \
                "$WINDOW_STEP" \
                "$(shell_quote "${SUBGROUP_DIR}/${group_id}.list")" \
                "$(shell_quote "$output_prefix")"
        } >> "$script_path"
    done < "$GROUP_IDS_FILE"

    {
        printf '\n'
        printf '%s\n' 'echo "[pi.sh] 开始计算各组合的 log2(Pi ratio)。"'
    } >> "$script_path"

    while IFS=$'\t' read -r ref_group check_group; do
        ratio_output="${PI_DIR}/Ref${ref_group}.vs.Check${check_group}.win${WINDOW_SIZE}.step${WINDOW_STEP}.log2_pi_ratio.pi"
        {
            printf '\n'
            printf 'echo %s\n' "$(shell_quote "[pi.sh] 计算 log2(${ref_group} Pi / ${check_group} Pi)")"
            printf 'Rscript %s %s %s %s %s %q %q %s\n' \
                "$(shell_quote "$calculator_path")" \
                "$(shell_quote "${PI_DIR}/${ref_group}.win${WINDOW_SIZE}.step${WINDOW_STEP}.windowed.pi")" \
                "$(shell_quote "${PI_DIR}/${check_group}.win${WINDOW_SIZE}.step${WINDOW_STEP}.windowed.pi")" \
                "$(shell_quote "$ref_group")" \
                "$(shell_quote "$check_group")" \
                "$WINDOW_SIZE" \
                "$WINDOW_STEP" \
                "$(shell_quote "$ratio_output")"
        } >> "$script_path"
    done < "$NORMALIZED_COMBINATIONS"

    {
        printf '\n'
        printf '%s\n' 'echo "[pi.sh] Pi 分析完成。"'
    } >> "$script_path"
    chmod +x "$script_path"
}

write_xpclr_script() {
    local script_path="${XPCLR_DIR}/xpclr.sh"
    local ref_group check_group comparison_prefix comparison_count

    comparison_count="$(awk 'END { print NR }' "$NORMALIZED_COMBINATIONS")"

    write_script_header "$script_path" "XP-CLR"
    {
        printf 'readonly CHROMOSOME_FILE=%s\n' "$(shell_quote "${XPCLR_DIR}/chromosomes.list")"
        printf 'readonly PARALLEL_XPCLR_REQUESTED=%s\n' "$(shell_quote "$PARALLEL_XPCLR")"
        printf 'readonly COMPARISON_COUNT=%q\n' "$comparison_count"
        printf '%s\n' 'command -v xpclr >/dev/null 2>&1 || { echo "错误：找不到 xpclr。" >&2; exit 127; }'
        printf '%s\n' 'if ((BASH_VERSINFO[0] < 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 3))); then echo "错误：并行 XP-CLR 需要 Bash 4.3 或更高版本。" >&2; exit 1; fi'
        if [[ "$INPUT_VCF" == *.gz || "$INPUT_VCF" == *.bgz ]]; then
            printf '%s\n' 'command -v gzip >/dev/null 2>&1 || { echo "错误：找不到 gzip。" >&2; exit 127; }'
        fi
        printf '\n'
        printf '%s\n' 'echo "[xpclr.sh] 按 VCF 中首次出现的顺序提取染色体 ID。"'
        if [[ "$INPUT_VCF" == *.gz || "$INPUT_VCF" == *.bgz ]]; then
            printf '%s\n' 'gzip -dc -- "$INPUT_VCF" | awk '"'"'!/^#/ && !seen[$1]++ { print $1 }'"'"' > "$CHROMOSOME_FILE"'
        else
            printf '%s\n' 'awk '"'"'!/^#/ && !seen[$1]++ { print $1 }'"'"' "$INPUT_VCF" > "$CHROMOSOME_FILE"'
        fi
        printf '%s\n' '[[ -s "$CHROMOSOME_FILE" ]] || { echo "错误：未能从 VCF 提取染色体 ID。" >&2; exit 1; }'
        printf '%s\n' 'chromosome_count="$(awk '"'"'END { print NR }'"'"' "$CHROMOSOME_FILE")"'
        printf '%s\n' 'if [[ "$PARALLEL_XPCLR_REQUESTED" == "auto" ]]; then'
        printf '%s\n' '    PARALLEL_XPCLR="$chromosome_count"'
        printf '%s\n' 'else'
        printf '%s\n' '    PARALLEL_XPCLR="$PARALLEL_XPCLR_REQUESTED"'
        printf '%s\n' 'fi'
        printf '%s\n' 'readonly PARALLEL_XPCLR'
        printf '%s\n' 'total_task_count=$((COMPARISON_COUNT * chromosome_count))'
        printf '%s\n' 'echo "[xpclr.sh] 共 ${COMPARISON_COUNT} 个组合、${chromosome_count} 条染色体、${total_task_count} 个子命令；全局并行数：${PARALLEL_XPCLR}。"'
        printf '\n'
        printf '%s\n' 'run_xpclr_task() {'
        printf '%s\n' '    local ref_group="$1"'
        printf '%s\n' '    local check_group="$2"'
        printf '%s\n' '    local chromosome_id="$3"'
        printf '%s\n' '    local comparison_prefix="$4"'
        printf '%s\n' '    local chromosome_output="${ANALYSIS_DIR}/${comparison_prefix}.${chromosome_id}.xpclr"'
        printf '%s\n' '    echo "[xpclr.sh] 开始：检测 Check=${check_group}（samplesA）相对于 Ref=${ref_group}（samplesB），${chromosome_id}"'
        printf '%s\n' '    if ! xpclr --format vcf --input "$INPUT_VCF" --samplesA "$SUBGROUP_DIR/${check_group}.list" --samplesB "$SUBGROUP_DIR/${ref_group}.list" --chr "$chromosome_id" --ld 0.95 --maxsnps 500 --size "$WINDOW_SIZE" --step "$WINDOW_STEP" --out "$chromosome_output"; then'
        printf '%s\n' '        echo "错误：XP-CLR 运行失败：Check=${check_group} 相对于 Ref=${ref_group}，${chromosome_id}" >&2'
        printf '%s\n' '        return 1'
        printf '%s\n' '    fi'
        printf '%s\n' '    if [[ ! -s "$chromosome_output" ]]; then'
        printf '%s\n' '        echo "错误：XP-CLR 输出不存在或为空：$chromosome_output" >&2'
        printf '%s\n' '        return 1'
        printf '%s\n' '    fi'
        printf '%s\n' '    echo "[xpclr.sh] 完成：Check=${check_group} 相对于 Ref=${ref_group}，${chromosome_id}"'
        printf '%s\n' '}'
        printf '\n'
        printf '%s\n' 'running_tasks=0'
        printf '%s\n' 'failed_tasks=0'
        printf '%s\n' 'launched_tasks=0'
        printf '%s\n' 'wait_for_one_task() {'
        printf '%s\n' '    local task_status=0'
        printf '%s\n' '    if wait -n; then'
        printf '%s\n' '        :'
        printf '%s\n' '    else'
        printf '%s\n' '        task_status=$?'
        printf '%s\n' '        failed_tasks=$((failed_tasks + 1))'
        printf '%s\n' '        echo "错误：一个 XP-CLR 子任务失败，退出状态：${task_status}。" >&2'
        printf '%s\n' '    fi'
        printf '%s\n' '    running_tasks=$((running_tasks - 1))'
        printf '%s\n' '}'
    } >> "$script_path"

    # 先把所有“组合 × 染色体”任务加入同一个全局并发池。
    while IFS=$'\t' read -r ref_group check_group; do
        comparison_prefix="Ref${ref_group}.vs.Check${check_group}.win${WINDOW_SIZE}.step${WINDOW_STEP}"
        {
            printf '\n'
            printf '%s\n' 'while IFS= read -r chromosome_id; do'
            printf '    run_xpclr_task %s %s "$chromosome_id" %s &\n' \
                "$(shell_quote "$ref_group")" \
                "$(shell_quote "$check_group")" \
                "$(shell_quote "$comparison_prefix")"
            printf '%s\n' '    running_tasks=$((running_tasks + 1))'
            printf '%s\n' '    launched_tasks=$((launched_tasks + 1))'
            printf '%s\n' '    if ((running_tasks >= PARALLEL_XPCLR)); then'
            printf '%s\n' '        wait_for_one_task'
            printf '%s\n' '    fi'
            printf '%s\n' 'done < "$CHROMOSOME_FILE"'
        } >> "$script_path"
    done < "$NORMALIZED_COMBINATIONS"

    {
        printf '\n'
        printf '%s\n' 'while ((running_tasks > 0)); do'
        printf '%s\n' '    wait_for_one_task'
        printf '%s\n' 'done'
        printf '%s\n' 'if ((failed_tasks > 0)); then'
        printf '%s\n' '    echo "错误：${failed_tasks} 个 XP-CLR 子任务失败，停止合并。" >&2'
        printf '%s\n' '    exit 1'
        printf '%s\n' 'fi'
        printf '%s\n' '[[ "$launched_tasks" -eq "$total_task_count" ]] || { echo "错误：计划 ${total_task_count} 个任务，实际启动 ${launched_tasks} 个。" >&2; exit 1; }'
        printf '%s\n' 'echo "[xpclr.sh] 全部 ${launched_tasks} 个 XP-CLR 子任务运行成功，开始按 VCF 染色体顺序合并。"'
    } >> "$script_path"

    # 全部并行任务成功后，再按 VCF 中的染色体顺序逐组合合并。
    while IFS=$'\t' read -r ref_group check_group; do
        comparison_prefix="Ref${ref_group}.vs.Check${check_group}.win${WINDOW_SIZE}.step${WINDOW_STEP}"
        {
            printf '\n'
            printf '%s\n' 'chromosome_outputs=()'
            printf '%s\n' 'while IFS= read -r chromosome_id; do'
            printf '    chromosome_output="$ANALYSIS_DIR/%s.${chromosome_id}.xpclr"\n' "$comparison_prefix"
            printf '%s\n' '    [[ -s "$chromosome_output" ]] || { echo "错误：XP-CLR 输出不存在或为空：$chromosome_output" >&2; exit 1; }'
            printf '%s\n' '    chromosome_outputs+=("$chromosome_output")'
            printf '%s\n' 'done < "$CHROMOSOME_FILE"'
            printf 'merged_output="$ANALYSIS_DIR/%s.xpclr"\n' "$comparison_prefix"
            printf '%s\n' 'temporary_output="${merged_output}.tmp.$$"'
            printf '%s\n' 'awk '"'"'FNR == 1 && NR != 1 { next } { print }'"'"' "${chromosome_outputs[@]}" > "$temporary_output"'
            printf '%s\n' 'mv -- "$temporary_output" "$merged_output"'
            printf '%s\n' 'echo "[xpclr.sh] 已合并：$merged_output"'
        } >> "$script_path"
    done < "$NORMALIZED_COMBINATIONS"

    {
        printf '\n'
        printf '%s\n' 'echo "[xpclr.sh] XP-CLR 分析完成。"'
    } >> "$script_path"
    chmod +x "$script_path"
}

INPUT_VCF_ARG=""
OUTPUT_DIR_ARG="$PWD"
ANALYSIS_PROJECT="all"
SAMPLE_GROUP_ARG=""
COMBINATION_ARG=""
WINDOW_SIZE="50000"
WINDOW_STEP="10000"
PARALLEL_XPCLR="auto"
DRY_RUN=false

while (($# > 0)); do
    case "$1" in
        --input-vcf)
            require_value "$1" "${2-}"
            INPUT_VCF_ARG="$2"
            shift 2
            ;;
        --input-vcf=*)
            INPUT_VCF_ARG="${1#*=}"
            [[ -n "$INPUT_VCF_ARG" ]] || die "参数 --input-vcf 需要提供一个值。"
            shift
            ;;
        --output-dir)
            require_value "$1" "${2-}"
            OUTPUT_DIR_ARG="$2"
            shift 2
            ;;
        --output-dir=*)
            OUTPUT_DIR_ARG="${1#*=}"
            [[ -n "$OUTPUT_DIR_ARG" ]] || die "参数 --output-dir 需要提供一个值。"
            shift
            ;;
        --analysis-project)
            require_value "$1" "${2-}"
            ANALYSIS_PROJECT="$2"
            shift 2
            ;;
        --analysis-project=*)
            ANALYSIS_PROJECT="${1#*=}"
            [[ -n "$ANALYSIS_PROJECT" ]] || die "参数 --analysis-project 需要提供一个值。"
            shift
            ;;
        --sample-group)
            require_value "$1" "${2-}"
            SAMPLE_GROUP_ARG="$2"
            shift 2
            ;;
        --sample-group=*)
            SAMPLE_GROUP_ARG="${1#*=}"
            [[ -n "$SAMPLE_GROUP_ARG" ]] || die "参数 --sample-group 需要提供一个值。"
            shift
            ;;
        --combination)
            require_value "$1" "${2-}"
            COMBINATION_ARG="$2"
            shift 2
            ;;
        --combination=*)
            COMBINATION_ARG="${1#*=}"
            [[ -n "$COMBINATION_ARG" ]] || die "参数 --combination 需要提供一个值。"
            shift
            ;;
        --window-size)
            require_value "$1" "${2-}"
            WINDOW_SIZE="$2"
            shift 2
            ;;
        --window-size=*)
            WINDOW_SIZE="${1#*=}"
            [[ -n "$WINDOW_SIZE" ]] || die "参数 --window-size 需要提供一个值。"
            shift
            ;;
        --window-step)
            require_value "$1" "${2-}"
            WINDOW_STEP="$2"
            shift 2
            ;;
        --window-step=*)
            WINDOW_STEP="${1#*=}"
            [[ -n "$WINDOW_STEP" ]] || die "参数 --window-step 需要提供一个值。"
            shift
            ;;
        --parallel-xpclr)
            require_value "$1" "${2-}"
            PARALLEL_XPCLR="$2"
            shift 2
            ;;
        --parallel-xpclr=*)
            PARALLEL_XPCLR="${1#*=}"
            [[ -n "$PARALLEL_XPCLR" ]] || die "参数 --parallel-xpclr 需要提供一个值。"
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            (($# == 0)) || die "不接受位置参数：$*"
            ;;
        -*)
            die "未知参数：$1。请使用 --help 查看帮助。"
            ;;
        *)
            die "不接受位置参数：$1。请使用 --help 查看帮助。"
            ;;
    esac
done

[[ -n "$INPUT_VCF_ARG" ]] || die "缺少必需参数 --input-vcf。"
[[ -n "$SAMPLE_GROUP_ARG" ]] || die "缺少必需参数 --sample-group。"
[[ -n "$COMBINATION_ARG" ]] || die "缺少必需参数 --combination。"
[[ "$ANALYSIS_PROJECT" =~ ^(fst|pi|xpclr|all)$ ]] ||
    die "--analysis-project 必须是 fst、pi、xpclr 或 all。"
[[ "$WINDOW_SIZE" =~ ^[1-9][0-9]*$ ]] ||
    die "--window-size 必须是正整数。"
[[ "$WINDOW_STEP" =~ ^[1-9][0-9]*$ ]] ||
    die "--window-step 必须是正整数。"
[[ "$PARALLEL_XPCLR" == "auto" || "$PARALLEL_XPCLR" =~ ^[1-9][0-9]*$ ]] ||
    die "--parallel-xpclr 必须是正整数。"

OUTPUT_DIR="$(realpath -m "$OUTPUT_DIR_ARG")"

mkdir -p "$OUTPUT_DIR"
LOG_FILE="${OUTPUT_DIR}/SelectiveSweep.log"
command -v tee >/dev/null 2>&1 || die "找不到 tee，无法写入日志文件。"
touch "$LOG_FILE" || die "无法创建或写入日志文件：${LOG_FILE}"
export SELECTIVE_SWEEP_LOG_FILE_ACTIVE="$LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

printf '\n========== SelectiveSweep 运行开始 ==========\n'
printf '  开始时间：%s\n' "$(date '+%Y-%m-%d %H:%M:%S %z')"
printf '  日志文件：%s\n' "$LOG_FILE"

INPUT_VCF="$(absolute_existing_file "$INPUT_VCF_ARG")"
SAMPLE_GROUP_FILE="$(absolute_existing_file "$SAMPLE_GROUP_ARG")"
COMBINATION_FILE="$(absolute_existing_file "$COMBINATION_ARG")"

if [[ "$INPUT_VCF" == *.gz || "$INPUT_VCF" == *.bgz ]]; then
    VCFTOOLS_INPUT_OPTION="--gzvcf"
else
    VCFTOOLS_INPUT_OPTION="--vcf"
fi

SUBGROUP_DIR="${OUTPUT_DIR}/subgroup"
mkdir -p "$SUBGROUP_DIR"

TEMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/selective-sweep.XXXXXXXX")"
trap 'rm -rf -- "$TEMP_DIR"' EXIT
NORMALIZED_SAMPLE_GROUP="${TEMP_DIR}/sample_group.tsv"
GROUP_IDS_FILE="${TEMP_DIR}/group_ids.list"
NORMALIZED_COMBINATIONS="${TEMP_DIR}/combinations.tsv"

awk '
BEGIN {
    FS = "[[:space:]]+"
    OFS = "\t"
}
{
    sub(/\r$/, "")
}
NR == 1 {
    next
}
/^[[:space:]]*$/ || /^[[:space:]]*#/ {
    next
}
NF != 2 {
    printf "错误：sample-group 第 %d 行应恰好有两列。\n", NR > "/dev/stderr"
    failed = 1
    next
}
$1 !~ /^[A-Za-z0-9_.:-]+$/ {
    printf "错误：sample-group 第 %d 行的样本 ID 含不支持字符：%s\n", NR, $1 > "/dev/stderr"
    failed = 1
    next
}
$2 !~ /^[A-Za-z0-9_.-]+$/ {
    printf "错误：sample-group 第 %d 行的亚群体 ID 含不支持字符：%s\n", NR, $2 > "/dev/stderr"
    failed = 1
    next
}
($1 in sample_group) && sample_group[$1] != $2 {
    printf "错误：样本 %s 被分配到多个亚群体：%s、%s。\n", $1, sample_group[$1], $2 > "/dev/stderr"
    failed = 1
    next
}
!seen_line[$1 SUBSEP $2]++ {
    sample_group[$1] = $2
    print $1, $2
    count++
}
END {
    if (count == 0) {
        print "错误：sample-group 没有有效的数据行。" > "/dev/stderr"
        failed = 1
    }
    exit failed
}
' "$SAMPLE_GROUP_FILE" > "$NORMALIZED_SAMPLE_GROUP"

cut -f2 "$NORMALIZED_SAMPLE_GROUP" | awk '!seen[$0]++' > "$GROUP_IDS_FILE"

while IFS= read -r group_id; do
    awk -F '\t' -v group="$group_id" '$2 == group { print $1 }' \
        "$NORMALIZED_SAMPLE_GROUP" > "${SUBGROUP_DIR}/${group_id}.list"
done < "$GROUP_IDS_FILE"

awk -v group_file="$GROUP_IDS_FILE" '
BEGIN {
    FS = "[[:space:]]+"
    OFS = "\t"
    while ((getline group < group_file) > 0) {
        groups[group] = 1
    }
    close(group_file)
}
{
    sub(/\r$/, "")
}
/^[[:space:]]*$/ || /^[[:space:]]*#/ {
    next
}
NF != 2 {
    printf "错误：combination 第 %d 行应恰好有两列。\n", NR > "/dev/stderr"
    failed = 1
    next
}
{
    first = $1
    second = $2
    first_lower = tolower(first)
    second_lower = tolower(second)
    is_header = data_count == 0 && (first_lower ~ /^(ref|reference|reference_group|group1|population1|pop1)$/ || second_lower ~ /^(check|test|test_group|group2|population2|pop2)$/)
    if (is_header) {
        next
    }
    if (!(first in groups)) {
        printf "错误：combination 第 %d 行的参考组不在 sample-group 中：%s\n", NR, first > "/dev/stderr"
        failed = 1
        next
    }
    if (!(second in groups)) {
        printf "错误：combination 第 %d 行的检验组不在 sample-group 中：%s\n", NR, second > "/dev/stderr"
        failed = 1
        next
    }
    if (first == second) {
        printf "错误：combination 第 %d 行不能比较同一个亚群体：%s\n", NR, first > "/dev/stderr"
        failed = 1
        next
    }
    key = first SUBSEP second
    if (!seen[key]++) {
        print first, second
        data_count++
    }
}
END {
    if (data_count == 0) {
        print "错误：combination 没有有效的比较组合。" > "/dev/stderr"
        failed = 1
    }
    exit failed
}
' "$COMBINATION_FILE" > "$NORMALIZED_COMBINATIONS"

FST_DIR="${OUTPUT_DIR}/fst"
PI_DIR="${OUTPUT_DIR}/pi"
XPCLR_DIR="${OUTPUT_DIR}/xpclr"

print_run_summary

case "$ANALYSIS_PROJECT" in
    fst)
        mkdir -p "$FST_DIR"
        write_fst_script
        GENERATED_SCRIPTS=("${FST_DIR}/fst.sh")
        ;;
    pi)
        mkdir -p "$PI_DIR"
        write_pi_script
        GENERATED_SCRIPTS=("${PI_DIR}/pi.sh")
        ;;
    xpclr)
        mkdir -p "$XPCLR_DIR"
        write_xpclr_script
        GENERATED_SCRIPTS=("${XPCLR_DIR}/xpclr.sh")
        ;;
    all)
        mkdir -p "$FST_DIR" "$PI_DIR" "$XPCLR_DIR"
        write_fst_script
        write_pi_script
        write_xpclr_script
        GENERATED_SCRIPTS=(
            "${FST_DIR}/fst.sh"
            "${PI_DIR}/pi.sh"
            "${XPCLR_DIR}/xpclr.sh"
        )
        ;;
esac

info "已生成亚群体样本列表：${SUBGROUP_DIR}"
for generated_script in "${GENERATED_SCRIPTS[@]}"; do
    info "已生成分析脚本：${generated_script}"
done

if [[ "$DRY_RUN" == true ]]; then
    info "dry-run 完成：未执行任何分析命令。"
    exit 0
fi

for generated_script in "${GENERATED_SCRIPTS[@]}"; do
    info "开始执行：${generated_script}"
    bash "$generated_script"
done

info "全部指定分析已完成。"
