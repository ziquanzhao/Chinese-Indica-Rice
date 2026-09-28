#!/usr/bin/env bash
set -euo pipefail

# ggComp.plus：以目标基因的cds、mRNA或链方向扩展gene区间为分析单元，缓存全部样本对N_diff/N_miss并开展群体比较。

SCRIPT_NAME="$(basename "$0")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELPER_SCRIPT_DIR="$SCRIPT_DIR/HelperScript"
COUNTER_SOURCE="$HELPER_SCRIPT_DIR/ggComp.plus.counter.cpp"
GROUP_COMPARE_SCRIPT="$HELPER_SCRIPT_DIR/ggComp.plus.group.compare.R"

INPUT_VCF=""
INPUT_GENE_ID=""
INPUT_GFF3=""
INPUT_SAMPLE_GROUP=""
INPUT_COMPARISON_GROUP=""
OUTPUT_DIR="."
PARALLEL=1
VARIANTION_SOURCE_INPUT="gene"
UPSTREAM=1000
DOWNSTREAM=500
HAPLOTYPE=false
HAPLOTYPE_TOP_NUMBER=5
HAPLOTYPE_TOP_NUMBER_SET=false
QUERY_FORMAT='%CHROM\t%POS[\t%GT]\n'

print_help() {
cat <<'EOF'
脚本功能:
  ggComp.plus以GFF3中目标基因的cds、mRNA或链方向扩展gene区间为分析单元，计算样本对的N_diff/N_miss和DSR_per_kb，并输出Group内均值、指定Group间FC、Welch t检验及BH-FDR。

核心逻辑:
  1. --input-gene-id提供带表头的单列GeneID列表，第一行固定作为表头跳过。GFF3第9列任意位置包含GeneID即匹配，区分大小写，不限制前后缀。
  2. --variantion-source默认gene，支持cds、mRNA、gene及逗号组合。cds使用全部CDS并集；mRNA使用第3列为mRNA的区间并集；gene以mRNA为基础，按--upstream和--downstream扩展。正链在左侧扩展upstream、右侧扩展downstream；负链则在左侧扩展downstream、右侧扩展upstream。
  3. GFF3第4、5列按1-based双端闭合处理；重叠或相邻区间先合并。每种variantion source独立建立GeneCounter/<source>/二进制缓存，不跨类别复用。
  4. 对每个位点和样本对：任一GT缺失或无法识别时N_miss加1；任一样本为杂合GT时该位点不参与DSR差异计数；仅两个样本均为纯合且基因型不同时N_diff加1。
  5. DSR_per_kb=N_diff×1000/(RegionLength-N_miss)。所有输入基因都进入Group比较，DSR_per_kb=0也参与均值。
  6. 分别使用各Group内部的无重复样本对计算Group.Mean_DSRperKb；不输出全样本均值或跨Group样本对均值。
  7. Ref.vs.Check.FC=Ref组内均值/Check组内均值；每个--comparison-group指定的Group对均执行双侧Welch t检验和按列BH-FDR。不输出跨Group样本对对任一Group的比较统计。
  8. 指定--haplotype true时，按每个基因区间的位点次序将纯合参考、纯合突变、杂合或缺失GT编码为0/2/N。杂合与缺失GT均在对应位点编码为N并保留样本。先按全样本频率确定AllSampleNo1~NoN；每个Group的SampleNo1~NoN严格对应相同编号的全样本单倍型，只计算其组内占比，不在Group内重新排序。

程序目录结构:
  <program-dir>/
    ggComp.plus.sh
    HelperScript/
      ggComp.plus.counter.cpp
      ggComp.plus.group.compare.R

用法:
  bash ggComp.plus.sh --input-vcf FILE --input-gene-id FILE --input-gff3 FILE --input-sample-group FILE --comparison-group FILE [options]

主要参数:
  --input-vcf FILE
      必需。单个已索引的VCF/VCF.GZ/BCF文件。

  --input-gene-id FILE
      必需。带表头的单列GeneID列表；第一行固定作为表头跳过，从第二行开始每行一个GeneID且不可重复。每个GeneID必须在GFF3中可以匹配到。

  --input-gff3 FILE
      必需。九列GFF3注释文件。第9列任意位置包含GeneID即匹配，GeneID前后允许任意字符串。

  --input-sample-group FILE
      必需。带表头的制表符分隔两列表，表头必须为SampleID和Group；SampleID必须唯一且必须存在于VCF。

  --comparison-group FILE
      必需。带表头的制表符分隔两列表，表头必须为Ref和Check。每行指定一个定向Group比较，Ref作为Ref.vs.Check.FC的分子；比较两组各自的组内DSR分布。

  --variantion-source VALUE
      可选。可选cds、mRNA、gene，或使用逗号指定多个，例如cds,mRNA,gene。默认gene。

  --upstream INT
      可选。gene模式中mRNA生物学上游扩展长度，bp，允许0。默认1000。对cds和mRNA模式无影响。

  --downstream INT
      可选。gene模式中mRNA生物学下游扩展长度，bp，允许0。默认500。对cds和mRNA模式无影响。

  --haplotype BOOL
      可选。必须指定true或false；true表示先确定全样本Top N单倍型，再计算这些固定单倍型在各Group中的占比；false表示不执行。默认false。

  --haplotype-top-number INT
      可选。单倍型统计的Top数量，必须为正整数，默认5；仅在--haplotype true时生效。

  --output-dir DIR
      可选。输出目录，默认当前目录。

  --parallel INT
      可选。C++基因计数线程数，默认1。

  -h, --help
      显示帮助并退出。

主要输出:
  <output-dir>/GeneCounter/
      按来源建立cds/、mRNA/、gene/子目录。每个子目录包含GeneIndex.tsv、GeneRegions.tsv、SamplePair.list、Ndiff/Nmiss二进制矩阵和完成标记。

  <output-dir>/GeneDSR/
      每种来源生成GroupComparison.Statistics.<source>.tsv和GroupComparison.Statistics.<source>.FDR.tsv。基因信息后给出各Group均值，以及--comparison-group指定的Group间FC和统计检验。

  <output-dir>/GeneHaplotype/
      指定--haplotype true时，每种来源生成Haplotype.Statistics.<source>.tsv；No1~NoN单元格仅输出百分比，同一编号在AllSample和各Group中表示同一种单倍型。

示例:
  bash ggComp.plus.sh --input-vcf cohort.vcf.gz --input-gene-id GeneID.list --input-gff3 OsNIP.longest.gff3 --input-sample-group SampleGroup.list --comparison-group Combination.group.list --variantion-source gene --upstream 1000 --downstream 500 --parallel 8

  bash ggComp.plus.sh --input-vcf cohort.vcf.gz --input-gene-id GeneID.list --input-gff3 OsNIP.longest.gff3 --input-sample-group SampleGroup.list --comparison-group Combination.group.list --variantion-source cds,mRNA,gene --parallel 8

  bash ggComp.plus.sh --input-vcf cohort.vcf.gz --input-gene-id GeneID.list --input-gff3 OsNIP.longest.gff3 --input-sample-group SampleGroup.list --comparison-group Combination.group.list --variantion-source cds --haplotype true --haplotype-top-number 5 --parallel 8

注意事项:
  1. 输入VCF/BCF应预先完成所需QC并包含可用GT；建议仅保留双等位SNP或Kmer。脚本不主动检查这些条件。
  2. VCF/BCF必须有可被bcftools读取的.tbi/.csi索引，并按CHROM、POS排序。
  3. gene模式不再使用另外的FAI，左端扩展坐标低于1时截到1；右端是否超出染色体由VCF索引查询结果自然限制。
  4. 单个基因区间内任一样本对的N_diff或N_miss超过65535时会报错。
  5. 断点记录会改变核心统计结果的解析后参数；未显式指定的可选参数按默认值记录。cds、mRNA、gene按来源独立记录断点，因此来源顺序不影响复用。
  6. Ref.vs.Check.FC=Ref组内均值/Check组内均值。分母为0且分子大于0时为Inf，两者均0时为NA。
  7. 单倍型仅使用--input-sample-group中的样本；杂合、缺失或无法识别的GT均在对应位点编码为N并保留样本。因此每个基因的SampleNumber等于相应分析范围内的输入样本数。AllSampleNo1~NoN按全样本频率排序，各Group的SampleNo1~NoN沿用同一单倍型顺序，不重新排序；无变异基因的No1为100.00%。
EOF
}

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log_info() { printf '[INFO %s] %s\n' "$(timestamp)" "$*" >&2; }
die() { printf '[ERROR %s] %s\n' "$(timestamp)" "$*" >&2; exit 1; }
is_positive_integer() { [[ "$1" =~ ^[1-9][0-9]*$ ]]; }
safe_token() { printf '%s' "$1" | sed 's/[^A-Za-z0-9._-]/_/g'; }
file_size() { stat -c '%s' -- "$1"; }

command -v getopt >/dev/null 2>&1 || { printf '[ERROR %s] 未找到getopt；请安装util-linux。\n' "$(timestamp)" >&2; exit 1; }
ARGS=$(getopt -o h --long help,input-vcf:,input-gene-id:,input-gff3:,input-sample-group:,comparison-group:,variantion-source:,upstream:,downstream:,haplotype:,haplotype-top-number:,output-dir:,parallel: -n "$SCRIPT_NAME" -- "$@") || { print_help >&2; exit 2; }
eval set -- "$ARGS"
while true; do
    case "$1" in
        --input-vcf) INPUT_VCF="$2"; shift 2 ;;
        --input-gene-id) INPUT_GENE_ID="$2"; shift 2 ;;
        --input-gff3) INPUT_GFF3="$2"; shift 2 ;;
        --input-sample-group) INPUT_SAMPLE_GROUP="$2"; shift 2 ;;
        --comparison-group) INPUT_COMPARISON_GROUP="$2"; shift 2 ;;
        --variantion-source) VARIANTION_SOURCE_INPUT="$2"; shift 2 ;;
        --upstream) UPSTREAM="$2"; shift 2 ;;
        --downstream) DOWNSTREAM="$2"; shift 2 ;;
        --haplotype) HAPLOTYPE="$2"; shift 2 ;;
        --haplotype-top-number) HAPLOTYPE_TOP_NUMBER="$2"; HAPLOTYPE_TOP_NUMBER_SET=true; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --parallel) PARALLEL="$2"; shift 2 ;;
        -h|--help) print_help; exit 0 ;;
        --) shift; break ;;
        *) die "内部参数解析错误: $1" ;;
    esac
done


for command_name in bcftools awk sort stat getopt flock g++ tee Rscript mktemp cmp cp mv rm mkdir touch date basename dirname paste sed head tail grep cat cut; do
    command -v "$command_name" >/dev/null 2>&1 || die "未找到必需程序$command_name；请安装后重新运行。"
done
[[ -d "$HELPER_SCRIPT_DIR" ]] || die "缺少辅助脚本目录: $HELPER_SCRIPT_DIR"
[[ -r "$COUNTER_SOURCE" ]] || die "缺少C++辅助程序: $COUNTER_SOURCE"
[[ -r "$GROUP_COMPARE_SCRIPT" ]] || die "缺少Group比较R辅助程序: $GROUP_COMPARE_SCRIPT"

BCFTOOLS_VERSION=$(bcftools --version 2>/dev/null | awk 'NR==1{print;exit}') || die "bcftools无法正常执行。"
GXX_VERSION=$(g++ --version 2>/dev/null | awk 'NR==1{print;exit}') || die "g++无法正常执行。"
PREFLIGHT_TMP=$(mktemp -d "${TMPDIR:-/tmp}/ggComp.plus.preflight.XXXXXX") || die "无法创建依赖预检临时目录。"
if ! printf '#include <omp.h>\nint main(){return omp_get_max_threads()>0?0:1;}\n' | g++ -x c++ -std=c++17 -fopenmp -o "$PREFLIGHT_TMP/openmp.test" - >/dev/null 2>&1 || ! "$PREFLIGHT_TMP/openmp.test"; then
    rm -rf -- "$PREFLIGHT_TMP"
    die "g++无法使用C++17/OpenMP编译并运行测试程序。"
fi
R_VERSION=$(Rscript --version 2>&1 | head -n 1) || { rm -rf -- "$PREFLIGHT_TMP"; die "Rscript无法正常执行。"; }
rm -rf -- "$PREFLIGHT_TMP"

[[ -n "$INPUT_VCF" ]] || die "缺少必需参数--input-vcf。"
[[ -n "$INPUT_GENE_ID" ]] || die "缺少必需参数--input-gene-id。"
[[ -n "$INPUT_GFF3" ]] || die "缺少必需参数--input-gff3。"
[[ -r "$INPUT_VCF" ]] || die "无法读取--input-vcf: $INPUT_VCF"
[[ -r "$INPUT_GENE_ID" ]] || die "无法读取--input-gene-id: $INPUT_GENE_ID"
[[ -r "$INPUT_GFF3" ]] || die "无法读取--input-gff3: $INPUT_GFF3"
[[ -n "$INPUT_SAMPLE_GROUP" ]] || die "缺少必需参数--input-sample-group。"
[[ -r "$INPUT_SAMPLE_GROUP" ]] || die "无法读取--input-sample-group: $INPUT_SAMPLE_GROUP"
[[ -n "$INPUT_COMPARISON_GROUP" ]] || die "缺少必需参数--comparison-group。"
[[ -r "$INPUT_COMPARISON_GROUP" ]] || die "无法读取--comparison-group: $INPUT_COMPARISON_GROUP"
is_positive_integer "$PARALLEL" || die "--parallel必须为正整数。"
is_positive_integer "$HAPLOTYPE_TOP_NUMBER" || die "--haplotype-top-number必须为正整数。"
[[ "$HAPLOTYPE" == true || "$HAPLOTYPE" == false ]] || die "--haplotype必须指定为true或false。"
[[ "$HAPLOTYPE" == true || "$HAPLOTYPE_TOP_NUMBER_SET" == false ]] || die "--haplotype-top-number仅可与--haplotype true同时使用。"
[[ "$UPSTREAM" =~ ^[0-9]+$ ]] || die "--upstream必须为非负整数。"
[[ "$DOWNSTREAM" =~ ^[0-9]+$ ]] || die "--downstream必须为非负整数。"
[[ -n "$VARIANTION_SOURCE_INPUT" ]] || die "--variantion-source不能为空。"
[[ "$VARIANTION_SOURCE_INPUT" != ,* && "$VARIANTION_SOURCE_INPUT" != *, && "$VARIANTION_SOURCE_INPUT" != *,,* ]] || die "--variantion-source逗号分隔格式无效: $VARIANTION_SOURCE_INPUT"
IFS=',' read -r -a VARIANTION_SOURCE_TOKENS <<< "$VARIANTION_SOURCE_INPUT"
VARIANTION_SOURCES=()
declare -A VARIANTION_SOURCE_SEEN=()
for source in "${VARIANTION_SOURCE_TOKENS[@]}"; do
    [[ "$source" == "cds" || "$source" == "mRNA" || "$source" == "gene" ]] || die "--variantion-source仅支持cds、mRNA、gene及其逗号组合: $VARIANTION_SOURCE_INPUT"
    [[ -z "${VARIANTION_SOURCE_SEEN[$source]+x}" ]] || die "--variantion-source中存在重复类别: $source"
    VARIANTION_SOURCE_SEEN["$source"]=1
    VARIANTION_SOURCES+=("$source")
done
mkdir -p -- "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
LOG_FILE="$OUTPUT_DIR/ggComp.plus.log"
touch "$LOG_FILE" || die "无法写入日志文件: $LOG_FILE"
if [[ -z "${GGCOMP_PLUS_LOG_INITIALIZED:-}" ]]; then
    exec 2> >(tee -a "$LOG_FILE" >&2)
    export GGCOMP_PLUS_LOG_INITIALIZED=1
fi

# 全局启动信息只由最外层进程输出一次；多来源子流程仅输出各自的来源提示。
if [[ -z "${GGCOMP_PLUS_SINGLE_SOURCE_CHILD:-}" ]]; then
    log_info "运行脚本: ggComp.plus GeneLevel"
    log_info "软件预检通过：$BCFTOOLS_VERSION；$GXX_VERSION；$R_VERSION。"
    log_info "请确保输入VCF/BCF已完成所需QC、包含可用GT，且建议仅保留双等位SNP或Kmer；脚本不会主动检查这些条件。"
    log_info "C++基因计数并行线程数: $PARALLEL。"
    log_info "DSR计数不使用杂合GT：任一样本为杂合时该位点不参与差异计数，仅比较两个样本均为纯合的位点。"
    if [[ "$HAPLOTYPE" == true ]]; then log_info "已启用单倍型分析：Top N=$HAPLOTYPE_TOP_NUMBER；杂合、缺失或无法识别的GT均在对应位点编码为N并保留样本。先按全样本频率确定No1~NoN，各Group沿用相同单倍型编号且不重新排序。"; fi
    log_info "--input-gene-id中的全部基因均进入分析；结果表只输出各Group组内均值和指定Group间比较统计。"
fi

# 多来源模式按来源依次启动完整的单来源流程，从而使缓存、断点和结果天然隔离。
if (( ${#VARIANTION_SOURCES[@]} > 1 )) && [[ -z "${GGCOMP_PLUS_SINGLE_SOURCE_CHILD:-}" ]]; then
    for source in "${VARIANTION_SOURCES[@]}"; do
        child_optional_args=(--haplotype "$HAPLOTYPE")
        [[ "$HAPLOTYPE" == true ]] && child_optional_args+=(--haplotype-top-number "$HAPLOTYPE_TOP_NUMBER")
        GGCOMP_PLUS_SINGLE_SOURCE_CHILD=1 bash "$SCRIPT_DIR/$SCRIPT_NAME" --input-vcf "$INPUT_VCF" --input-gene-id "$INPUT_GENE_ID" --input-gff3 "$INPUT_GFF3" --input-sample-group "$INPUT_SAMPLE_GROUP" --comparison-group "$INPUT_COMPARISON_GROUP" --variantion-source "$source" --upstream "$UPSTREAM" --downstream "$DOWNSTREAM" --output-dir "$OUTPUT_DIR" --parallel "$PARALLEL" "${child_optional_args[@]}"
    done
    exit 0
fi
VARIANTION_SOURCE="${VARIANTION_SOURCES[0]}"
GENE_EXPANSION_IDENTITY=""
[[ "$VARIANTION_SOURCE" == "gene" ]] && GENE_EXPANSION_IDENTITY=".Upstream.$UPSTREAM.Downstream.$DOWNSTREAM"

log_info "本次--variantion-source=$VARIANTION_SOURCE；每种来源建立独立GeneCounter缓存并分别执行Group比较。"
if [[ "$VARIANTION_SOURCE" == "gene" ]]; then log_info "gene模式区间扩展：upstream=$UPSTREAM bp，downstream=$DOWNSTREAM bp；正链=左侧扩展upstream/右侧扩展downstream，负链相反。"; fi
bcftools view -h "$INPUT_VCF" >/dev/null 2>&1 || die "--input-vcf必须是单个可被bcftools读取的VCF/VCF.GZ/BCF文件；不再支持VCF列表: $INPUT_VCF"
[[ -f "${INPUT_VCF}.csi" || -f "${INPUT_VCF}.tbi" ]] || die "--input-vcf必须已有.csi或.tbi索引: $INPUT_VCF"

TMP_DIR=$(mktemp -d "$OUTPUT_DIR/.ggComp.plus.tmp.XXXXXX")
export TMPDIR="$TMP_DIR"
cleanup() { rm -rf -- "$TMP_DIR"; }
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

CHECKPOINT_FILE="$OUTPUT_DIR/ggComp.plus.checkpoint.tsv"
CHECKPOINT_PREAMBLE="$TMP_DIR/checkpoint.preamble.tsv"
CHECKPOINT_INPUT_VCF=$(basename "$INPUT_VCF")
for checkpoint_parameter_value in "$CHECKPOINT_INPUT_VCF" "$INPUT_GENE_ID" "$INPUT_GFF3" "$INPUT_SAMPLE_GROUP" "$INPUT_COMPARISON_GROUP" "$UPSTREAM" "$DOWNSTREAM" "$HAPLOTYPE" "$HAPLOTYPE_TOP_NUMBER"; do
    [[ "$checkpoint_parameter_value" != *$'\t'* && "$checkpoint_parameter_value" != *$'\n'* && "$checkpoint_parameter_value" != *$'\r'* ]] || die "命令参数值不能包含制表符或换行符，无法写入断点参数快照。"
done
{
    printf '# ggComp.plus checkpoint\n'
    printf '# PARAMETER\t--input-vcf\t%s\n' "$CHECKPOINT_INPUT_VCF"
    printf '# PARAMETER\t--input-gene-id\t%s\n' "$INPUT_GENE_ID"
    printf '# PARAMETER\t--input-gff3\t%s\n' "$INPUT_GFF3"
    printf '# PARAMETER\t--input-sample-group\t%s\n' "$INPUT_SAMPLE_GROUP"
    printf '# PARAMETER\t--comparison-group\t%s\n' "$INPUT_COMPARISON_GROUP"
    printf '# PARAMETER\t--upstream\t%s\n' "$UPSTREAM"
    printf '# PARAMETER\t--downstream\t%s\n' "$DOWNSTREAM"
    printf '# PARAMETER\t--haplotype\t%s\n' "$HAPLOTYPE"
    printf '# PARAMETER\t--haplotype-top-number\t%s\n' "$HAPLOTYPE_TOP_NUMBER"
    printf '# STAGE\tKEY\n'
} > "$CHECKPOINT_PREAMBLE"
CHECKPOINT_PREAMBLE_LINES=$(awk 'END{print NR+0}' "$CHECKPOINT_PREAMBLE")
if [[ ! -f "$CHECKPOINT_FILE" ]]; then
    cp -- "$CHECKPOINT_PREAMBLE" "$CHECKPOINT_FILE"
else
    head -n "$CHECKPOINT_PREAMBLE_LINES" "$CHECKPOINT_FILE" > "$TMP_DIR/checkpoint.existing.preamble.tsv"
    cmp -s "$TMP_DIR/checkpoint.existing.preamble.tsv" "$CHECKPOINT_PREAMBLE" || die "现有断点文件的格式或完整命令参数快照与本次运行不一致: $CHECKPOINT_FILE；请恢复原参数，或使用新的--output-dir。"
fi

BIN_DIR="$OUTPUT_DIR/.ggComp.plus.bin"
COUNTER_BIN="$BIN_DIR/ggComp.plus.counter"
mkdir -p -- "$BIN_DIR"
if [[ ! -x "$COUNTER_BIN" || "$COUNTER_SOURCE" -nt "$COUNTER_BIN" ]]; then
    log_info "编译C++基因区间计数辅助程序。"
    g++ -O3 -std=c++17 -fopenmp -Wall -Wextra -Wpedantic "$COUNTER_SOURCE" -o "$COUNTER_BIN" || die "C++辅助程序编译失败。"
fi

# 第9列任意位置只要包含目标GeneID即视为匹配；按variantion source提取并合并区间。
GENE_BED_NORMALIZED="$TMP_DIR/GeneRegion.$VARIANTION_SOURCE.normalized.tsv"
GENE_CHROMS="$TMP_DIR/gene.chromosome.list"
"$COUNTER_BIN" extract-regions --gene-id "$INPUT_GENE_ID" --gff3 "$INPUT_GFF3" --variantion-source "$VARIANTION_SOURCE" --upstream "$UPSTREAM" --downstream "$DOWNSTREAM" --output "$GENE_BED_NORMALIZED" || die "从GFF3提取目标基因$VARIANTION_SOURCE区间失败。"
GENE_COUNT=$(awk 'END{print NR+0}' "$GENE_BED_NORMALIZED")
awk -F '\t' '!seen[$1]++{print $1}' "$GENE_BED_NORMALIZED" > "$GENE_CHROMS"
CHROM_COUNT=$(awk 'END{print NR+0}' "$GENE_CHROMS")
REGION_INTERVAL_COUNT=$(awk -F '\t' '{n=split($6,a,",");total+=n}END{print total+0}' "$GENE_BED_NORMALIZED")
REGION_TOTAL_LENGTH=$(awk -F '\t' '{total+=$5}END{printf "%.0f",total}' "$GENE_BED_NORMALIZED")
(( GENE_COUNT > 0 && CHROM_COUNT > 0 && REGION_INTERVAL_COUNT > 0 && REGION_TOTAL_LENGTH > 0 )) || die "标准化基因$VARIANTION_SOURCE区间为空。"
log_info "目标基因$VARIANTION_SOURCE区间检查通过：GeneID=$GENE_COUNT，合并后片段=$REGION_INTERVAL_COUNT，总长度=$REGION_TOTAL_LENGTH bp，分布于$CHROM_COUNT条染色体。"

VCF_MAP="$TMP_DIR/vcf.map.tsv"
bcftools index -s "$INPUT_VCF" | awk -F '\t' 'NF{print $1}' > "$TMP_DIR/input.index.chroms" || die "无法读取VCF/BCF索引。"
awk 'NR==FNR{v[$1]=1;next}!($1 in v){print $1}' "$TMP_DIR/input.index.chroms" "$GENE_CHROMS" > "$TMP_DIR/missing.chroms"
[[ ! -s "$TMP_DIR/missing.chroms" ]] || die "--input-vcf未覆盖以下目标基因区间染色体: $(paste -sd, "$TMP_DIR/missing.chroms")"
awk -v path="$INPUT_VCF" -v OFS='\t' '{print $1,path}' "$GENE_CHROMS" > "$VCF_MAP"
log_info "--input-vcf已确认为单个已索引VCF/BCF；仅查询目标基因$VARIANTION_SOURCE区间。"

ALL_SAMPLES="$TMP_DIR/all.samples.list"
bcftools query -l "$INPUT_VCF" > "$ALL_SAMPLES" || die "无法读取VCF样本名称。"
ALL_SAMPLE_COUNT=$(awk 'END{print NR+0}' "$ALL_SAMPLES")
(( ALL_SAMPLE_COUNT >= 2 )) || die "VCF至少需要两个样本。"

REQUESTED_SAMPLES="$TMP_DIR/requested.samples.list"
SAMPLE_GROUP_TABLE="$TMP_DIR/SampleGroup.normalized.tsv"
SAMPLE_GROUP_ORDER="$TMP_DIR/SampleGroup.order.list"
awk -F '\t' -v OFS='\t' 'NR==1{sub(/\r$/, "", $NF);if(NF!=2||$1!="SampleID"||$2!="Group"){print "样本分组表表头必须严格为SampleID和Group两列" > "/dev/stderr";bad=1}next}/^[[:space:]]*$/{next}{sub(/\r$/, "", $NF);if(NF!=2||$1==""||$2==""){printf "样本分组表第%d行必须是非空的SampleID和Group两列\n",NR > "/dev/stderr";bad=1;next}if(seen_sample[$1]++){printf "SampleID重复: %s\n",$1 > "/dev/stderr";bad=1;next}if($1!~/^[A-Za-z0-9._-]+$/){printf "SampleID包含不安全字符: %s\n",$1 > "/dev/stderr";bad=1;next}if($2!~/^[A-Za-z0-9][A-Za-z0-9._-]*$/||$2=="."||$2==".."){printf "Group包含不安全字符: %s\n",$2 > "/dev/stderr";bad=1;next}print $1,$2;rows++}END{if(NR<2||rows==0){print "样本分组表没有数据行" > "/dev/stderr";bad=1}exit bad?1:0}' "$INPUT_SAMPLE_GROUP" > "$SAMPLE_GROUP_TABLE" || die "--input-sample-group格式检查失败。"
cut -f1 "$SAMPLE_GROUP_TABLE" > "$REQUESTED_SAMPLES"
awk -F '\t' '!seen[$2]++{print $2}' "$SAMPLE_GROUP_TABLE" > "$SAMPLE_GROUP_ORDER"
awk -F '\t' '{count[$2]++}END{for(group in count)if(count[group]<2)print group}' "$SAMPLE_GROUP_TABLE" > "$TMP_DIR/groups.with_too_few_samples"
[[ ! -s "$TMP_DIR/groups.with_too_few_samples" ]] || die "以下Group少于2个样本，无法形成样本对: $(paste -sd, "$TMP_DIR/groups.with_too_few_samples")"
awk 'NR==FNR{v[$1]=1;next}!($1 in v){print $1}' "$ALL_SAMPLES" "$REQUESTED_SAMPLES" > "$TMP_DIR/missing.samples"
[[ ! -s "$TMP_DIR/missing.samples" ]] || die "以下样本不在VCF中: $(paste -sd, "$TMP_DIR/missing.samples")"

COMPARISON_GROUP_TABLE="$TMP_DIR/ComparisonGroup.normalized.tsv"
awk -F '\t' -v OFS='\t' 'NR==1{sub(/\r$/, "", $NF);if(NF!=2||$1!="Ref"||$2!="Check"){print "比较组表表头必须严格为Ref和Check两列" > "/dev/stderr";bad=1}next}/^[[:space:]]*$/{next}{sub(/\r$/, "", $NF);if(NF!=2||$1==""||$2==""){printf "比较组表第%d行必须是非空的Ref和Check两列\n",NR > "/dev/stderr";bad=1;next}if($1==$2){printf "Ref与Check不能相同: %s\n",$1 > "/dev/stderr";bad=1;next}if($1!~/^[A-Za-z0-9][A-Za-z0-9._-]*$/||$2!~/^[A-Za-z0-9][A-Za-z0-9._-]*$/){printf "Ref或Check包含不安全字符: %s,%s\n",$1,$2 > "/dev/stderr";bad=1;next}unordered=($1<$2?$1 SUBSEP $2:$2 SUBSEP $1);if(seen_unordered[unordered]++){printf "同一对Group不能重复或反向重复: %s,%s\n",$1,$2 > "/dev/stderr";bad=1;next}label=$1 $2;if(seen_label[label]++){printf "比较组拼接名称重复: %s\n",label > "/dev/stderr";bad=1;next}print $1,$2;rows++}END{if(NR<2||rows==0){print "比较组表没有数据行" > "/dev/stderr";bad=1}exit bad?1:0}' "$INPUT_COMPARISON_GROUP" > "$COMPARISON_GROUP_TABLE" || die "--comparison-group格式检查失败。"
awk 'NR==FNR{group[$1]=1;next}!($1 in group){print "Ref\t" $1}!($2 in group){print "Check\t" $2}' "$SAMPLE_GROUP_ORDER" "$COMPARISON_GROUP_TABLE" > "$TMP_DIR/comparison.groups.missing"
[[ ! -s "$TMP_DIR/comparison.groups.missing" ]] || die "--comparison-group中以下Group未在--input-sample-group中出现: $(awk -F '\t' '{print $1"="$2}' "$TMP_DIR/comparison.groups.missing" | paste -sd, -)"
GROUP_COMPARISON_COUNT=$(awk 'END{print NR+0}' "$COMPARISON_GROUP_TABLE")
log_info "样本分组表检查通过：样本数=$(awk 'END{print NR+0}' "$REQUESTED_SAMPLES")，Group数=$(awk 'END{print NR+0}' "$SAMPLE_GROUP_ORDER")。"
log_info "比较组表检查通过：定向比较数=$GROUP_COMPARISON_COUNT；只计算指定Ref组内对Check组内的FC/Welch t检验。"

EXISTING_CACHE_CONFIG="$OUTPUT_DIR/GeneCounter/$VARIANTION_SOURCE/cache.config.tsv"
CACHE_SAMPLES="$TMP_DIR/cache.samples.list"
if [[ -f "$EXISTING_CACHE_CONFIG" ]]; then
    awk -F '\t' '$1=="SAMPLES"{n=split($2,a,",");for(i=1;i<=n;i++)print a[i];found++}END{exit found==1?0:1}' "$EXISTING_CACHE_CONFIG" > "$CACHE_SAMPLES" || die "现有GeneCounter缓存缺少唯一SAMPLES顺序记录: $EXISTING_CACHE_CONFIG"
    awk 'NF!=1||seen[$1]++{bad=1}END{exit bad?1:0}' "$CACHE_SAMPLES" || die "现有GeneCounter缓存的SAMPLES顺序记录无效。"
    awk 'NR==FNR{v[$1]=1;next}!($1 in v){print $1}' "$CACHE_SAMPLES" "$REQUESTED_SAMPLES" > "$TMP_DIR/samples.not_in_background_cache"
    [[ ! -s "$TMP_DIR/samples.not_in_background_cache" ]] || die "以下样本未在现有全群体GeneCounter数据库中出现，无法继续: $(paste -sd, "$TMP_DIR/samples.not_in_background_cache")"
    awk 'NR==FNR{v[$1]=1;next}!($1 in v){print $1}' "$CACHE_SAMPLES" "$ALL_SAMPLES" > "$TMP_DIR/cache.missing.samples"
    awk 'NR==FNR{v[$1]=1;next}!($1 in v){print $1}' "$ALL_SAMPLES" "$CACHE_SAMPLES" > "$TMP_DIR/cache.extra.samples"
    [[ ! -s "$TMP_DIR/cache.missing.samples" && ! -s "$TMP_DIR/cache.extra.samples" ]] || die "现有GeneCounter不是当前VCF全部样本形成的全群体缓存；请检查--output-dir。"
    log_info "读取$VARIANTION_SOURCE GeneCounter缓存中的持久化样本顺序：$(awk 'END{print NR+0}' "$CACHE_SAMPLES")个样本。"
else
    cp -- "$ALL_SAMPLES" "$CACHE_SAMPLES"
    log_info "尚无全群体缓存；将以VCF全部样本顺序建立GeneCounter缓存主顺序。"
fi

awk -F '\t' -v OFS='\t' 'NR==FNR{group[$1]=$2;next}$1 in group{print $1,group[$1]}' "$SAMPLE_GROUP_TABLE" "$CACHE_SAMPLES" > "$SAMPLE_GROUP_TABLE.cache_sorted"
[[ $(awk 'END{print NR+0}' "$SAMPLE_GROUP_TABLE.cache_sorted") -eq $(awk 'END{print NR+0}' "$SAMPLE_GROUP_TABLE") ]] || die "--input-sample-group无法按照GeneCounter样本顺序完整重排。"
mv -f -- "$SAMPLE_GROUP_TABLE.cache_sorted" "$SAMPLE_GROUP_TABLE"
log_info "分组表中的全部样本已按照GeneCounter持久化样本顺序重排；Group输出顺序仍按--input-sample-group中Group首次出现顺序。"

CACHE_SAMPLE_COUNT=$(awk 'END{print NR+0}' "$CACHE_SAMPLES")
CACHE_CSV=$(paste -sd, "$CACHE_SAMPLES")
CACHE_PAIR_COUNT=$(( CACHE_SAMPLE_COUNT * (CACHE_SAMPLE_COUNT - 1) / 2 ))
log_info "GeneCounter全群体缓存范围：VCF全部$CACHE_SAMPLE_COUNT个样本，共$CACHE_PAIR_COUNT个无重复样本对。"

GENE_DIR="$OUTPUT_DIR/GeneCounter/$VARIANTION_SOURCE"
PAIR_GROUP_DIR="$GENE_DIR/PairGroups"
SHARD_DIR="$GENE_DIR/.ChromosomeShards"
PAIR_LIST="$GENE_DIR/SamplePair.list"
GROUP_INDEX="$GENE_DIR/PairGroup.index.tsv"
GENE_INDEX="$GENE_DIR/GeneIndex.tsv"
GENE_CDS_CACHE="$GENE_DIR/GeneRegions.tsv"
CACHE_CONFIG="$GENE_DIR/cache.config.tsv"
CACHE_COMPLETE="$GENE_DIR/cache.complete.tsv"
mkdir -p -- "$GENE_DIR" "$PAIR_GROUP_DIR" "$BIN_DIR"

if [[ -f "$GENE_INDEX" || -f "$GENE_CDS_CACHE" ]]; then
    [[ -f "$GENE_INDEX" && -f "$GENE_CDS_CACHE" ]] || die "现有$VARIANTION_SOURCE GeneCounter缺少GeneIndex.tsv或GeneRegions.tsv；请使用新的--output-dir。"
    head -n 1 "$GENE_INDEX" | grep -qx $'GeneIndex\tChr\tStart\tEnd\tGeneID\tRegionLength\tRegionIntervals\tVariantionSource' || die "现有GeneIndex.tsv不是当前区间缓存格式；请使用新的--output-dir。"
    [[ $(awk 'END{print NR+0}' "$GENE_INDEX") -eq $((GENE_COUNT+1)) ]] || die "现有GeneIndex.tsv基因数与本次GeneID列表不一致；请使用新的--output-dir。"
    cmp -s "$GENE_CDS_CACHE" "$GENE_BED_NORMALIZED" || die "现有$VARIANTION_SOURCE GeneCounter的GeneID/GFF3区间定义与本次输入不一致；请使用新的--output-dir。"
    CACHE_GENE_COUNT="$GENE_COUNT"
    log_info "本次提取的$GENE_COUNT个基因$VARIANTION_SOURCE区间与持久化GeneRegions定义完全一致。"
else
    CACHE_GENE_COUNT="$GENE_COUNT"
    cp -- "$GENE_BED_NORMALIZED" "$GENE_CDS_CACHE.tmp"
    mv -f -- "$GENE_CDS_CACHE.tmp" "$GENE_CDS_CACHE"
    { printf 'GeneIndex\tChr\tStart\tEnd\tGeneID\tRegionLength\tRegionIntervals\tVariantionSource\n'; awk -F '\t' -v OFS='\t' -v source="$VARIANTION_SOURCE" '{print NR,$1,$2,$3,$4,$5,$6,source}' "$GENE_CDS_CACHE"; } > "$GENE_INDEX.tmp"
    mv -f -- "$GENE_INDEX.tmp" "$GENE_INDEX"
    log_info "以本次输入的$CACHE_GENE_COUNT个GeneID及其$VARIANTION_SOURCE区间建立全群体GeneCounter主顺序。"
fi

CACHE_GENE_BED="$GENE_CDS_CACHE"
awk '!seen[$1]++{print $1}' "$CACHE_GENE_BED" > "$GENE_CHROMS"
CHROM_COUNT=$(awk 'END{print NR+0}' "$GENE_CHROMS")

REQUESTED_CACHE_CONFIG="$TMP_DIR/cache.config.requested.tsv"
{
    printf 'FORMAT\tggComp.plus.gene-region-cache\n'
    printf 'VARIANTION_SOURCE\t%s\n' "$VARIANTION_SOURCE"
    if [[ -n "$GENE_EXPANSION_IDENTITY" ]]; then printf 'UPSTREAM\t%s\nDOWNSTREAM\t%s\n' "$UPSTREAM" "$DOWNSTREAM"; fi
    printf 'GENOTYPE_RULE\tmerged_region_1based_closed_missing_priority_then_ignore_heterozygous_then_homozygous_difference\n'
    printf 'GFF3_MATCH_RULE\tcolumn9_contains_GeneID_anywhere_case_sensitive\n'
    printf 'GENE_ID\t%s\n' "$(basename "$INPUT_GENE_ID")"
    printf 'GFF3\t%s\n' "$(basename "$INPUT_GFF3")"
    printf 'GENE_COUNT\t%s\n' "$CACHE_GENE_COUNT"
    printf 'SAMPLES\t%s\n' "$(paste -sd, "$CACHE_SAMPLES")"
    while IFS=$'\t' read -r chrom vcf_path; do printf 'VCF\t%s\t%s\n' "$chrom" "$(basename "$vcf_path")"; done < "$VCF_MAP"
} > "$REQUESTED_CACHE_CONFIG"
if [[ -f "$CACHE_CONFIG" ]]; then
    cmp -s "$CACHE_CONFIG" "$REQUESTED_CACHE_CONFIG" || die "现有GeneCounter缓存文件名或分析配置与本次输入不一致；请使用新--output-dir，或备份后删除$GENE_DIR。"
else
    cp -- "$REQUESTED_CACHE_CONFIG" "$CACHE_CONFIG"
fi
awk '{sample[NR]=$1}END{for(i=1;i<NR;i++)for(j=i+1;j<=NR;j++)print sample[i]","sample[j]}' "$CACHE_SAMPLES" > "$PAIR_LIST.tmp"
awk -F ',' '{key=$1"."$2;if(seen[key]++)print key}' "$PAIR_LIST.tmp" > "$TMP_DIR/pair.filename.collisions"
[[ ! -s "$TMP_DIR/pair.filename.collisions" ]] || die "样本ID组合产生重复结果文件名。"
mv -f -- "$PAIR_LIST.tmp" "$PAIR_LIST"
: > "$GROUP_INDEX.tmp"
while IFS= read -r first_sample; do
    group_pairs="$PAIR_GROUP_DIR/${first_sample}.SamplePair.list"
    awk -F ',' -v sample="$first_sample" '$1==sample' "$PAIR_LIST" > "$group_pairs.tmp"
    if [[ -s "$group_pairs.tmp" ]]; then
        mv -f -- "$group_pairs.tmp" "$group_pairs"
        printf '%s\t%s\t%s\n' "$first_sample" "$(awk 'END{print NR+0}' "$group_pairs")" "$group_pairs" >> "$GROUP_INDEX.tmp"
    else
        rm -f -- "$group_pairs.tmp"
    fi
done < "$CACHE_SAMPLES"
mv -f -- "$GROUP_INDEX.tmp" "$GROUP_INDEX"
log_info "全群体样本对列表: $PAIR_LIST；划分为$((CACHE_SAMPLE_COUNT-1))个基因计数组。"

declare -A CHECKPOINT_INDEX=()
reload_checkpoint_index() {
    local checkpoint_stage checkpoint_key
    CHECKPOINT_INDEX=()
    while IFS=$'\t' read -r checkpoint_stage checkpoint_key _; do
        [[ -n "$checkpoint_stage" && "$checkpoint_stage" != \#* && -n "$checkpoint_key" ]] || continue
        CHECKPOINT_INDEX["$checkpoint_stage"$'\t'"$checkpoint_key"]=1
    done < "$CHECKPOINT_FILE"
}
reload_checkpoint_index

checkpoint_matches() {
    local stage="$1" key="$2" output_file="$3" index_key
    [[ -f "$output_file" ]] || return 1
    index_key="$stage"$'\t'"$key"
    [[ -n "${CHECKPOINT_INDEX[$index_key]+x}" ]]
}
record_checkpoint() {
    local stage="$1" key="$2" output_file="$3" index_key
    [[ -f "$output_file" ]] || return 1
    exec 9>> "$CHECKPOINT_FILE" || return 1
    flock -x 9 || { exec 9>&-; return 1; }
    printf '%s\t%s\n' "$stage" "$key" >&9
    flock -u 9
    exec 9>&-
    if declare -p CHECKPOINT_INDEX >/dev/null 2>&1; then index_key="$stage"$'\t'"$key";CHECKPOINT_INDEX["$index_key"]=1; fi
}
deduplicate_checkpoints() {
    local compact_tmp="$CHECKPOINT_FILE.compact.tmp"
    exec 9>> "$CHECKPOINT_FILE" || return 1
    flock -x 9 || { exec 9>&-; return 1; }
    awk -F '\t' -v OFS='\t' '/^#/{print;next}$1!=""&&$2!=""&&!seen[$1 SUBSEP $2]++{print $1,$2}' "$CHECKPOINT_FILE" > "$compact_tmp" || { flock -u 9;exec 9>&-;return 1; }
    mv -f -- "$compact_tmp" "$CHECKPOINT_FILE"
    flock -u 9
    exec 9>&-
    reload_checkpoint_index
}
export CHECKPOINT_FILE
export -f record_checkpoint timestamp

GENE_CACHE_STAGE="GENE_REGION_CACHE:$VARIANTION_SOURCE"
CACHE_GENE_ID_TAG=$(awk -F '\t' '$1=="GENE_ID"{print $2;exit}' "$CACHE_CONFIG")
CACHE_GFF3_TAG=$(awk -F '\t' '$1=="GFF3"{print $2;exit}' "$CACHE_CONFIG")
[[ -n "$CACHE_GENE_ID_TAG" && -n "$CACHE_GFF3_TAG" ]] || die "GeneCounter缓存配置缺少GENE_ID或GFF3记录。"
GENE_CACHE_KEY="$(safe_token "$CACHE_GENE_ID_TAG.$CACHE_GFF3_TAG.$VARIANTION_SOURCE$GENE_EXPANSION_IDENTITY.HeterozygousExcluded")"
MEM_AVAILABLE_BYTES=$(awk '$1=="MemAvailable:"{print $2*1024;exit}' /proc/meminfo)
[[ "$MEM_AVAILABLE_BYTES" =~ ^[0-9]+$ ]] || die "无法读取当前可用内存。"
if [[ -r /sys/fs/cgroup/memory.max && -r /sys/fs/cgroup/memory.current ]]; then
    CGROUP_MEMORY_MAX=$(< /sys/fs/cgroup/memory.max)
    CGROUP_MEMORY_CURRENT=$(< /sys/fs/cgroup/memory.current)
    if [[ "$CGROUP_MEMORY_MAX" =~ ^[0-9]+$ && "$CGROUP_MEMORY_CURRENT" =~ ^[0-9]+$ && $CGROUP_MEMORY_MAX -gt $CGROUP_MEMORY_CURRENT ]]; then
        CGROUP_MEMORY_AVAILABLE=$(( CGROUP_MEMORY_MAX - CGROUP_MEMORY_CURRENT ))
        (( MEM_AVAILABLE_BYTES <= CGROUP_MEMORY_AVAILABLE )) || MEM_AVAILABLE_BYTES=$CGROUP_MEMORY_AVAILABLE
    fi
fi
MEMORY_BUDGET=$(( MEM_AVAILABLE_BYTES / 2 ))
MAX_MEMORY_BUDGET=$(( 10 * 1024 * 1024 * 1024 ))
(( MEMORY_BUDGET <= MAX_MEMORY_BUDGET )) || MEMORY_BUDGET=$MAX_MEMORY_BUDGET
(( MEMORY_BUDGET > 0 )) || die "基因计数内存预算为0。"
log_info "基因计数内存预算=min(当前MemAvailable的一半,10 GiB)=$MEMORY_BUDGET bytes。"

validate_gene_matrices() {
    local first_sample ndiff nmiss
    [[ -f "$GENE_INDEX" ]] || return 1
    while IFS=$'\t' read -r first_sample _; do
        ndiff="$PAIR_GROUP_DIR/$first_sample/$first_sample.Ndiff.u16le.bin"
        nmiss="$PAIR_GROUP_DIR/$first_sample/$first_sample.Nmiss.u16le.bin"
        [[ -f "$ndiff" && -f "$nmiss" ]] || return 1
    done < "$GROUP_INDEX"
}
validate_gene_cache_artifacts() {
    [[ -f "$CACHE_COMPLETE" ]] || return 1
    awk -F '\t' -v source="$VARIANTION_SOURCE" -v genes="$CACHE_GENE_COUNT" -v pairs="$CACHE_PAIR_COUNT" -v groups="$((CACHE_SAMPLE_COUNT-1))" '$1=="FORMAT"&&$2=="ggComp.plus.gene-region-cache.complete"{f++}$1=="VARIANTION_SOURCE"&&$2==source{s++}$1=="GENE_COUNT"&&$2==genes{g++}$1=="PAIR_COUNT"&&$2==pairs{p++}$1=="GROUP_COUNT"&&$2==groups{r++}END{exit !(NR==5&&f==1&&s==1&&g==1&&p==1&&r==1)}' "$CACHE_COMPLETE" || return 1
}
validate_gene_cache() { validate_gene_cache_artifacts && checkpoint_matches "$GENE_CACHE_STAGE" "$GENE_CACHE_KEY" "$CACHE_COMPLETE"; }

if validate_gene_cache; then
    log_info "基因二进制缓存和断点记录检查通过，跳过VCF重新读取。"
else
    mkdir -p -- "$SHARD_DIR"
    chrom_index=0
    while IFS= read -r chrom; do
        ((chrom_index+=1))
        vcf_path=$(awk -F '\t' -v chrom="$chrom" '$1==chrom{print $2;exit}' "$VCF_MAP")
        chrom_gene_bed="$TMP_DIR/$(printf '%06d' "$chrom_index").genes.bed"
        awk -F '\t' -v chrom="$chrom" '$1==chrom' "$CACHE_GENE_BED" > "$chrom_gene_bed"
        chrom_gene_count=$(awk 'END{print NR+0}' "$chrom_gene_bed")
        (( chrom_gene_count > 0 )) || die "染色体$chrom没有基因。"
        chrom_cds_regions="$TMP_DIR/$(printf '%06d' "$chrom_index").$VARIANTION_SOURCE.regions.tsv"
        awk -F '\t' -v OFS='\t' '{n=split($6,intervals,",");for(i=1;i<=n;i++){split(intervals[i],position,"-");print $1,position[1],position[2]}}' "$chrom_gene_bed" | sort -t $'\t' -k2,2n -k3,3n | awk -F '\t' -v OFS='\t' 'NR==1{chr=$1;start=$2;end=$3;next}$2<=end+1{if($3>end)end=$3;next}{print chr,start,end;chr=$1;start=$2;end=$3}END{if(NR)print chr,start,end}' > "$chrom_cds_regions" || die "染色体$chrom $VARIANTION_SOURCE查询区间生成失败。"
        [[ -s "$chrom_cds_regions" ]] || die "染色体$chrom $VARIANTION_SOURCE查询区间为空。"
        shard_name=$(printf '%06d' "$chrom_index")
        shard="$SHARD_DIR/$shard_name"
        if [[ -f "$shard/chromosome.meta.tsv" ]]; then
            log_info "断点续跑：跳过已完成的染色体基因计数$chrom。"
            continue
        fi
        shard_tmp="$SHARD_DIR/.${shard_name}.tmp"
        if [[ -f "$shard_tmp/chromosome.meta.tsv" ]]; then
            shard_tmp_valid=true
            while IFS=$'\t' read -r first_sample columns _; do
                expected=$(( chrom_gene_count * columns * 2 ))
                if [[ ! -f "$shard_tmp/$first_sample.Ndiff.u16le.bin" || ! -f "$shard_tmp/$first_sample.Nmiss.u16le.bin" ]] ||
                   [[ $(file_size "$shard_tmp/$first_sample.Ndiff.u16le.bin") -ne $expected || $(file_size "$shard_tmp/$first_sample.Nmiss.u16le.bin") -ne $expected ]]; then
                    shard_tmp_valid=false
                    break
                fi
            done < "$GROUP_INDEX"
            if [[ "$shard_tmp_valid" == true ]]; then
                rm -rf -- "$shard"
                mv -- "$shard_tmp" "$shard"
                log_info "断点续跑：恢复已完成但尚未正式提交的染色体基因计数$chrom。"
                continue
            fi
            log_info "检测到不完整的染色体临时分片$chrom，将清理后重新计算。"
        fi
        rm -rf -- "$shard_tmp"
        mkdir -p -- "$shard_tmp"
        log_info "仅查询染色体$chrom的合并$VARIANTION_SOURCE区间，并累计$chrom_gene_count个基因的N_diff/N_miss。"
        bcftools query -R "$chrom_cds_regions" -s "$CACHE_CSV" -f "$QUERY_FORMAT" "$vcf_path" | "$COUNTER_BIN" count --chrom "$chrom" --gene-bed "$chrom_gene_bed" --memory-bytes "$MEMORY_BUDGET" --threads "$PARALLEL" --samples "$CACHE_SAMPLES" --output-dir "$shard_tmp" || die "染色体$chrom基因$VARIANTION_SOURCE计数失败。"
        while IFS=$'\t' read -r first_sample columns _; do
            expected=$(( chrom_gene_count * columns * 2 ))
            [[ $(file_size "$shard_tmp/$first_sample.Ndiff.u16le.bin") -eq $expected && $(file_size "$shard_tmp/$first_sample.Nmiss.u16le.bin") -eq $expected ]] || die "染色体$chrom二进制分片大小验证失败: $first_sample"
        done < "$GROUP_INDEX"
        rm -rf -- "$shard"
        mv -- "$shard_tmp" "$shard"
    done < "$GENE_CHROMS"

    log_info "按目标基因$VARIANTION_SOURCE染色体顺序合并基因计数分片。"
    while IFS=$'\t' read -r first_sample columns group_pairs; do
        final_group="$PAIR_GROUP_DIR/$first_sample"
        mkdir -p -- "$final_group"
        ndiff_tmp="$final_group/.${first_sample}.Ndiff.tmp"
        nmiss_tmp="$final_group/.${first_sample}.Nmiss.tmp"
        : > "$ndiff_tmp"
        : > "$nmiss_tmp"
        for ((index=1;index<=CHROM_COUNT;index++)); do
            shard="$SHARD_DIR/$(printf '%06d' "$index")"
            cat "$shard/$first_sample.Ndiff.u16le.bin" >> "$ndiff_tmp"
            cat "$shard/$first_sample.Nmiss.u16le.bin" >> "$nmiss_tmp"
        done
        expected=$(( CACHE_GENE_COUNT * columns * 2 ))
        [[ $(file_size "$ndiff_tmp") -eq $expected && $(file_size "$nmiss_tmp") -eq $expected ]] || die "最终基因矩阵大小验证失败: $first_sample"
        mv -f -- "$ndiff_tmp" "$final_group/$first_sample.Ndiff.u16le.bin"
        mv -f -- "$nmiss_tmp" "$final_group/$first_sample.Nmiss.u16le.bin"
        cp -- "$group_pairs" "$final_group/$first_sample.SamplePair.list"
    done < "$GROUP_INDEX"
    validate_gene_matrices || die "最终基因二进制矩阵验证失败。"
    {
        printf 'FORMAT\tggComp.plus.gene-region-cache.complete\n'
        printf 'VARIANTION_SOURCE\t%s\n' "$VARIANTION_SOURCE"
        printf 'GENE_COUNT\t%s\n' "$CACHE_GENE_COUNT"
        printf 'PAIR_COUNT\t%s\n' "$CACHE_PAIR_COUNT"
        printf 'GROUP_COUNT\t%s\n' "$((CACHE_SAMPLE_COUNT-1))"
    } > "$CACHE_COMPLETE.tmp"
    mv -f -- "$CACHE_COMPLETE.tmp" "$CACHE_COMPLETE"
    record_checkpoint "$GENE_CACHE_STAGE" "$GENE_CACHE_KEY" "$CACHE_COMPLETE" || die "无法写入GENE_CACHE断点。"
    validate_gene_cache || die "基因缓存与断点联合验证失败。"
    rm -rf -- "$SHARD_DIR"
    log_info "基因缓存完成且验证通过；已删除成功合并的染色体临时分片。"
fi

# 一次流式汇总全部Group组内统计。
run_normal_group_comparison() {
NORMAL_DIR="$OUTPUT_DIR/GeneDSR"
mkdir -p -- "$NORMAL_DIR"
NORMAL_OUTPUT="$NORMAL_DIR/GroupComparison.Statistics.$VARIANTION_SOURCE.tsv"
NORMAL_FDR_OUTPUT="$NORMAL_DIR/GroupComparison.Statistics.$VARIANTION_SOURCE.FDR.tsv"
NORMAL_PAIR_COUNTS="$NORMAL_DIR/Group.SamplePair.Count.tsv"
GROUP_COUNT=$(awk 'END{print NR+0}' "$SAMPLE_GROUP_ORDER")
(( GROUP_COUNT >= 2 )) || die "至少需要2个Group才能执行组间Welch t检验。"
GROUP_COMPARISON_COUNT=$(awk 'END{print NR+0}' "$COMPARISON_GROUP_TABLE")
EXPECTED_NORMAL_COLUMNS=$(( 5 + GROUP_COUNT + 2 * GROUP_COMPARISON_COUNT ))

printf 'Group\tSampleCount\tSamplePairCount\n' > "$NORMAL_PAIR_COUNTS.tmp"
while IFS= read -r group; do
    sample_count=$(awk -F '\t' -v target="$group" '$2==target{n++}END{print n+0}' "$SAMPLE_GROUP_TABLE")
    pair_count=$(( sample_count * (sample_count - 1) / 2 ))
    printf '%s\t%s\t%s\n' "$group" "$sample_count" "$pair_count" >> "$NORMAL_PAIR_COUNTS.tmp"
done < "$SAMPLE_GROUP_ORDER"
mv -f -- "$NORMAL_PAIR_COUNTS.tmp" "$NORMAL_PAIR_COUNTS"

# 第一层断点只对应原始均值、FC和Welch t检验表；完成后立即落盘并记录。
NORMAL_RAW_CONFIG="$NORMAL_DIR/GroupComparison.Statistics.$VARIANTION_SOURCE.config.tsv"
NORMAL_RAW_MARKER="$NORMAL_DIR/GroupComparison.Statistics.$VARIANTION_SOURCE.complete.tsv"
NORMAL_RAW_STAGE="GROUP_RAW_STATISTICS:$VARIANTION_SOURCE"
NORMAL_RAW_KEY="$(safe_token "AllInputGenes.$VARIANTION_SOURCE$GENE_EXPANSION_IDENTITY.HeterozygousExcluded.GeneID.$(basename "$INPUT_GENE_ID").GFF3.$(basename "$INPUT_GFF3").SampleGroup.$(basename "$INPUT_SAMPLE_GROUP").ComparisonGroup.$(basename "$INPUT_COMPARISON_GROUP").WithinGroupMean.RefVsCheck.FC.WelchT.OutputSchemaV2")"
NORMAL_RAW_REQUESTED_CONFIG="$TMP_DIR/GroupComparison.$VARIANTION_SOURCE.raw.requested.tsv"
{
    printf 'FORMAT\tggComp.plus.group-region-raw-statistics\n'
    printf 'VARIANTION_SOURCE\t%s\n' "$VARIANTION_SOURCE"
    if [[ -n "$GENE_EXPANSION_IDENTITY" ]]; then printf 'UPSTREAM\t%s\nDOWNSTREAM\t%s\n' "$UPSTREAM" "$DOWNSTREAM"; fi
    printf 'INPUT_GENE_ID\t%s\n' "$(basename "$INPUT_GENE_ID")"
    printf 'INPUT_GFF3\t%s\n' "$(basename "$INPUT_GFF3")"
    printf 'INPUT_SAMPLE_GROUP\t%s\n' "$(basename "$INPUT_SAMPLE_GROUP")"
    printf 'INPUT_COMPARISON_GROUP\t%s\n' "$(basename "$INPUT_COMPARISON_GROUP")"
    printf 'GENE_COUNT\t%s\n' "$GENE_COUNT"
    printf 'GROUP_COUNT\t%s\n' "$GROUP_COUNT"
    printf 'OUTPUT_SCHEMA\tgroup_means_and_ref_vs_check_only\n'
    while IFS=$'\t' read -r sample group; do printf 'SAMPLE_GROUP\t%s\t%s\n' "$sample" "$group"; done < "$SAMPLE_GROUP_TABLE"
    while IFS=$'\t' read -r ref check; do printf 'COMPARISON_GROUP\t%s\t%s\n' "$ref" "$check"; done < "$COMPARISON_GROUP_TABLE"
} > "$NORMAL_RAW_REQUESTED_CONFIG"

if checkpoint_matches "$NORMAL_RAW_STAGE" "$NORMAL_RAW_KEY" "$NORMAL_RAW_MARKER" && [[ -f "$NORMAL_OUTPUT" && -f "$NORMAL_RAW_CONFIG" ]] && [[ $(awk 'END{print NR+0}' "$NORMAL_OUTPUT") -eq $((GENE_COUNT+1)) ]] && awk -F '\t' -v expected="$EXPECTED_NORMAL_COLUMNS" 'NR==1{exit NF==expected?0:1}' "$NORMAL_OUTPUT" && head -n 1 "$NORMAL_OUTPUT" | grep -q '\.FC' && head -n 1 "$NORMAL_OUTPUT" | grep -q '\.Welch\.t' && cmp -s "$NORMAL_RAW_CONFIG" "$NORMAL_RAW_REQUESTED_CONFIG"; then
    log_info "断点续跑：复用已完成的Group内均值、Ref对Check FC和原始Welch t检验表: $NORMAL_OUTPUT"
else
    NORMAL_SUMMARY_TMP="$TMP_DIR/Normal.GroupComparison.sufficient-statistics.tsv"
    log_info "开始一次流式读取$VARIANTION_SOURCE GeneCounter缓存，计算$GENE_COUNT个基因的$GROUP_COUNT个Group内均值；DSR_per_kb=0参与统计，不再汇总全样本或跨Group样本对。"
    "$COUNTER_BIN" compare-groups --gene-bed "$CACHE_GENE_BED" --selected-gene-bed "$GENE_BED_NORMALIZED" --group-index "$GROUP_INDEX" --pair-group-dir "$PAIR_GROUP_DIR" --sample-group "$SAMPLE_GROUP_TABLE" --group-order "$SAMPLE_GROUP_ORDER" --output "$NORMAL_SUMMARY_TMP" || die "从GeneCounter缓存汇总Group内统计失败。"
    log_info "开始计算指定Group组合的Ref组内对Check组内FC和双侧Welch t检验p-value。"
    Rscript "$GROUP_COMPARE_SCRIPT" --raw "$NORMAL_SUMMARY_TMP" "$SAMPLE_GROUP_ORDER" "$COMPARISON_GROUP_TABLE" "$NORMAL_OUTPUT.tmp" || die "计算指定Group比较的逐基因FC及Welch t检验失败。"
    [[ $(awk 'END{print NR+0}' "$NORMAL_OUTPUT.tmp") -eq $((GENE_COUNT+1)) ]] || die "Group比较原始结果表行数检查失败。"
    awk -F '\t' -v expected="$EXPECTED_NORMAL_COLUMNS" 'NR==1{exit NF==expected?0:1}' "$NORMAL_OUTPUT.tmp" || die "Group比较原始结果表列数检查失败。"
    head -n 1 "$NORMAL_OUTPUT.tmp" | cut -f1-5 | grep -qx $'GeneID\tChr\tStart\tEnd\tRegionLength' || die "Group比较原始结果表基因信息表头异常。"
    head -n 1 "$NORMAL_OUTPUT.tmp" | grep -q '\.FC' || die "Group比较原始结果表缺少.FC列。"
    head -n 1 "$NORMAL_OUTPUT.tmp" | grep -qv 'AllSample\.Mean_DSRperKb' || die "Group比较原始结果表仍包含AllSample.Mean_DSRperKb列。"
    head -n 1 "$NORMAL_OUTPUT.tmp" | grep -q '\.Welch\.t' || die "Group比较原始结果表缺少.Welch.t列。"
    mv -f -- "$NORMAL_OUTPUT.tmp" "$NORMAL_OUTPUT"
    cp -- "$NORMAL_RAW_REQUESTED_CONFIG" "$NORMAL_RAW_CONFIG.tmp";mv -f -- "$NORMAL_RAW_CONFIG.tmp" "$NORMAL_RAW_CONFIG"
    {
        printf 'FORMAT\tggComp.plus.group-region-raw-statistics.complete\n'
        printf 'VARIANTION_SOURCE\t%s\n' "$VARIANTION_SOURCE"
        printf 'STAGE\t%s\n' "$NORMAL_RAW_STAGE"
        printf 'KEY\t%s\n' "$NORMAL_RAW_KEY"
        printf 'GENE_COUNT\t%s\n' "$GENE_COUNT"
        printf 'GROUP_COUNT\t%s\n' "$GROUP_COUNT"
        printf 'COMPARISON_COUNT\t%s\n' "$GROUP_COMPARISON_COUNT"
        printf 'OUTPUT\t%s\n' "$(basename "$NORMAL_OUTPUT")"
    } > "$NORMAL_RAW_MARKER.tmp"
    mv -f -- "$NORMAL_RAW_MARKER.tmp" "$NORMAL_RAW_MARKER"
    record_checkpoint "$NORMAL_RAW_STAGE" "$NORMAL_RAW_KEY" "$NORMAL_RAW_MARKER" || die "无法写入Group原始统计表断点。"
    log_info "Group原始统计表已完成并写入断点: $NORMAL_RAW_MARKER"
fi

# 第二层断点只对应固定的Benjamini-Hochberg校正，可直接读取原始统计表。
NORMAL_FDR_CONFIG="$NORMAL_DIR/GroupComparison.Statistics.$VARIANTION_SOURCE.FDR.config.tsv"
NORMAL_FDR_MARKER="$NORMAL_DIR/GroupComparison.Statistics.$VARIANTION_SOURCE.FDR.complete.tsv"
NORMAL_FDR_STAGE="GROUP_BH_FDR:$VARIANTION_SOURCE"
NORMAL_FDR_KEY="$(safe_token "$NORMAL_RAW_KEY.BenjaminiHochberg")"
NORMAL_FDR_REQUESTED_CONFIG="$TMP_DIR/GroupComparison.$VARIANTION_SOURCE.fdr.requested.tsv"
{
    printf 'FORMAT\tggComp.plus.group-region-bh-fdr\n'
    printf 'VARIANTION_SOURCE\t%s\n' "$VARIANTION_SOURCE"
    printf 'RAW_STAGE\t%s\n' "$NORMAL_RAW_STAGE"
    printf 'RAW_KEY\t%s\n' "$NORMAL_RAW_KEY"
    printf 'GENE_COUNT\t%s\n' "$GENE_COUNT"
    printf 'GROUP_COUNT\t%s\n' "$GROUP_COUNT"
    printf 'METHOD\tBenjamini-Hochberg\n'
} > "$NORMAL_FDR_REQUESTED_CONFIG"

if checkpoint_matches "$NORMAL_FDR_STAGE" "$NORMAL_FDR_KEY" "$NORMAL_FDR_MARKER" && [[ -f "$NORMAL_FDR_OUTPUT" && -f "$NORMAL_FDR_CONFIG" ]] && [[ $(awk 'END{print NR+0}' "$NORMAL_FDR_OUTPUT") -eq $((GENE_COUNT+1)) ]] && awk -F '\t' -v expected="$EXPECTED_NORMAL_COLUMNS" 'NR==1{exit NF==expected?0:1}' "$NORMAL_FDR_OUTPUT" && head -n 1 "$NORMAL_FDR_OUTPUT" | grep -q '\.FC' && head -n 1 "$NORMAL_FDR_OUTPUT" | grep -q '\.Welch\.t\.FDR' && cmp -s "$NORMAL_FDR_CONFIG" "$NORMAL_FDR_REQUESTED_CONFIG"; then
    log_info "断点续跑：复用已完成的BH-FDR校正表: $NORMAL_FDR_OUTPUT"
else
    log_info "开始从原始统计表执行BH-FDR校正；每个Group比较均在全部$GENE_COUNT个基因范围内独立校正。"
    Rscript "$GROUP_COMPARE_SCRIPT" --fdr "$NORMAL_OUTPUT" "$NORMAL_FDR_OUTPUT.tmp" || die "从原始p-value表生成BH-FDR校正表失败。"
    [[ $(awk 'END{print NR+0}' "$NORMAL_FDR_OUTPUT.tmp") -eq $((GENE_COUNT+1)) ]] || die "BH-FDR校正结果表行数检查失败。"
    awk -F '\t' -v expected="$EXPECTED_NORMAL_COLUMNS" 'NR==1{exit NF==expected?0:1}' "$NORMAL_FDR_OUTPUT.tmp" || die "BH-FDR校正结果表列数检查失败。"
    head -n 1 "$NORMAL_FDR_OUTPUT.tmp" | grep -q '\.FC' || die "BH-FDR校正结果表缺少.FC列。"
    head -n 1 "$NORMAL_FDR_OUTPUT.tmp" | grep -qv 'AllSample\.Mean_DSRperKb' || die "BH-FDR校正结果表仍包含AllSample.Mean_DSRperKb列。"
    head -n 1 "$NORMAL_FDR_OUTPUT.tmp" | grep -q '\.Welch\.t\.FDR' || die "BH-FDR校正结果表缺少.Welch.t.FDR列。"
    mv -f -- "$NORMAL_FDR_OUTPUT.tmp" "$NORMAL_FDR_OUTPUT"
    cp -- "$NORMAL_FDR_REQUESTED_CONFIG" "$NORMAL_FDR_CONFIG.tmp";mv -f -- "$NORMAL_FDR_CONFIG.tmp" "$NORMAL_FDR_CONFIG"
    {
        printf 'FORMAT\tggComp.plus.group-region-bh-fdr.complete\n'
        printf 'VARIANTION_SOURCE\t%s\n' "$VARIANTION_SOURCE"
        printf 'STAGE\t%s\n' "$NORMAL_FDR_STAGE"
        printf 'KEY\t%s\n' "$NORMAL_FDR_KEY"
        printf 'GENE_COUNT\t%s\n' "$GENE_COUNT"
        printf 'GROUP_COUNT\t%s\n' "$GROUP_COUNT"
        printf 'OUTPUT\t%s\n' "$(basename "$NORMAL_FDR_OUTPUT")"
    } > "$NORMAL_FDR_MARKER.tmp"
    mv -f -- "$NORMAL_FDR_MARKER.tmp" "$NORMAL_FDR_MARKER"
    record_checkpoint "$NORMAL_FDR_STAGE" "$NORMAL_FDR_KEY" "$NORMAL_FDR_MARKER" || die "无法写入BH-FDR校正断点。"
fi
deduplicate_checkpoints || die "整理Group比较断点失败。"
log_info "$VARIANTION_SOURCE来源Group比较完成：每个基因包含$GROUP_COUNT个Group内均值，以及$GROUP_COMPARISON_COUNT组Ref对Check的FC/双侧Welch t检验和对应BH-FDR。"
log_info "原始p-value结果表: $NORMAL_OUTPUT"
log_info "BH-FDR校正结果表: $NORMAL_FDR_OUTPUT"
log_info "各Group样本数与样本对数: $NORMAL_PAIR_COUNTS"
}

run_haplotype_analysis() {
HAPLOTYPE_DIR="$OUTPUT_DIR/GeneHaplotype"
HAPLOTYPE_OUTPUT="$HAPLOTYPE_DIR/Haplotype.Statistics.$VARIANTION_SOURCE.tsv"
HAPLOTYPE_CONFIG="$HAPLOTYPE_DIR/Haplotype.Statistics.$VARIANTION_SOURCE.config.tsv"
HAPLOTYPE_MARKER="$HAPLOTYPE_DIR/Haplotype.Statistics.$VARIANTION_SOURCE.complete.tsv"
HAPLOTYPE_SHARD_DIR="$HAPLOTYPE_DIR/.Haplotype.$VARIANTION_SOURCE.ChromosomeShards"
HAPLOTYPE_STAGE="HAPLOTYPE_STATISTICS:$VARIANTION_SOURCE"
HAPLOTYPE_KEY="$(safe_token "AllInputSamples.$VARIANTION_SOURCE$GENE_EXPANSION_IDENTITY.GeneID.$(basename "$INPUT_GENE_ID").GFF3.$(basename "$INPUT_GFF3").SampleGroup.$(basename "$INPUT_SAMPLE_GROUP").Top.$HAPLOTYPE_TOP_NUMBER.FixedGlobalNoOrder.EncodeHeterozygousAsN")"
HAPLOTYPE_REQUESTED_CONFIG="$TMP_DIR/Haplotype.Statistics.$VARIANTION_SOURCE.requested.tsv"
GROUP_COUNT=$(awk 'END{print NR+0}' "$SAMPLE_GROUP_ORDER")
EXPECTED_HAPLOTYPE_COLUMNS=$(( 5 + (1 + HAPLOTYPE_TOP_NUMBER) * (1 + GROUP_COUNT) ))
mkdir -p -- "$HAPLOTYPE_DIR"
{
    printf 'FORMAT\tggComp.plus.haplotype-statistics\n'
    printf 'VARIANTION_SOURCE\t%s\n' "$VARIANTION_SOURCE"
    if [[ -n "$GENE_EXPANSION_IDENTITY" ]]; then printf 'UPSTREAM\t%s\nDOWNSTREAM\t%s\n' "$UPSTREAM" "$DOWNSTREAM"; fi
    printf 'INPUT_VCF\t%s\n' "$(basename "$INPUT_VCF")"
    printf 'INPUT_GENE_ID\t%s\n' "$(basename "$INPUT_GENE_ID")"
    printf 'INPUT_GFF3\t%s\n' "$(basename "$INPUT_GFF3")"
    printf 'INPUT_SAMPLE_GROUP\t%s\n' "$(basename "$INPUT_SAMPLE_GROUP")"
    printf 'HAPLOTYPE_TOP_NUMBER\t%s\n' "$HAPLOTYPE_TOP_NUMBER"
    printf 'MISSING_RULE\tencode_missing_or_undefined_GT_as_N_and_keep_sample\n'
    printf 'HETEROZYGOUS_RULE\tencode_heterozygous_GT_as_N_and_keep_sample\n'
    printf 'GT_ENCODING\thom_ref=0,hom_alt=2,heterozygous=N,missing=N\n'
    printf 'RANKING_RULE\tAllSample_ranks_all_haplotypes_by_frequency_then_code;each_Group_keeps_the_same_AllSample_No_order_without_reranking\n'
    while IFS=$'\t' read -r sample group; do printf 'SAMPLE_GROUP\t%s\t%s\n' "$sample" "$group"; done < "$SAMPLE_GROUP_TABLE"
} > "$HAPLOTYPE_REQUESTED_CONFIG"

if checkpoint_matches "$HAPLOTYPE_STAGE" "$HAPLOTYPE_KEY" "$HAPLOTYPE_MARKER" &&
   [[ -f "$HAPLOTYPE_OUTPUT" && -f "$HAPLOTYPE_CONFIG" ]] &&
   [[ $(awk 'END{print NR+0}' "$HAPLOTYPE_OUTPUT") -eq $((GENE_COUNT+1)) ]] &&
   awk -F '\t' -v expected="$EXPECTED_HAPLOTYPE_COLUMNS" 'NR==1{exit NF==expected?0:1}' "$HAPLOTYPE_OUTPUT" &&
   cmp -s "$HAPLOTYPE_CONFIG" "$HAPLOTYPE_REQUESTED_CONFIG"; then
    log_info "断点续跑：复用已完成的$VARIANTION_SOURCE单倍型统计表: $HAPLOTYPE_OUTPUT"
    return
fi

if [[ -d "$HAPLOTYPE_SHARD_DIR" ]]; then
    if [[ ! -f "$HAPLOTYPE_SHARD_DIR/request.config.tsv" ]] || ! cmp -s "$HAPLOTYPE_SHARD_DIR/request.config.tsv" "$HAPLOTYPE_REQUESTED_CONFIG"; then
        log_info "单倍型参数或样本分组与当前临时分片不一致，清理后重新计算$VARIANTION_SOURCE染色体分片。"
        rm -rf -- "$HAPLOTYPE_SHARD_DIR"
    fi
fi
mkdir -p -- "$HAPLOTYPE_SHARD_DIR"
cp -- "$HAPLOTYPE_REQUESTED_CONFIG" "$HAPLOTYPE_SHARD_DIR/request.config.tsv"
HAPLOTYPE_SAMPLES="$TMP_DIR/Haplotype.samples.list"
cut -f1 "$SAMPLE_GROUP_TABLE" > "$HAPLOTYPE_SAMPLES"
HAPLOTYPE_CSV=$(paste -sd, "$HAPLOTYPE_SAMPLES")
log_info "开始$VARIANTION_SOURCE单倍型分析：输入样本=$(awk 'END{print NR+0}' "$HAPLOTYPE_SAMPLES")，Group=$GROUP_COUNT，Top N=$HAPLOTYPE_TOP_NUMBER；杂合、缺失或无法识别的GT均在对应位点编码为N并保留样本；先按全样本频率确定No1~NoN，各Group严格沿用相同单倍型编号且不重新排序。"
chrom_index=0
while IFS= read -r chrom; do
    ((chrom_index+=1))
    shard_name=$(printf '%06d' "$chrom_index")
    haplotype_shard="$HAPLOTYPE_SHARD_DIR/$shard_name.tsv"
    chrom_gene_bed="$TMP_DIR/Haplotype.$shard_name.genes.tsv"
    chrom_regions="$TMP_DIR/Haplotype.$shard_name.regions.tsv"
    awk -F '\t' -v chrom="$chrom" '$1==chrom' "$CACHE_GENE_BED" > "$chrom_gene_bed"
    chrom_gene_count=$(awk 'END{print NR+0}' "$chrom_gene_bed")
    awk -F '\t' -v OFS='\t' '{n=split($6,intervals,",");for(i=1;i<=n;i++){split(intervals[i],position,"-");print $1,position[1],position[2]}}' "$chrom_gene_bed" | sort -t $'\t' -k2,2n -k3,3n | awk -F '\t' -v OFS='\t' 'NR==1{chr=$1;start=$2;end=$3;next}$2<=end+1{if($3>end)end=$3;next}{print chr,start,end;chr=$1;start=$2;end=$3}END{if(NR)print chr,start,end}' > "$chrom_regions" || die "单倍型染色体$chrom查询区间生成失败。"
    if [[ -f "$haplotype_shard" ]] &&
       [[ $(awk 'END{print NR+0}' "$haplotype_shard") -eq $((chrom_gene_count+1)) ]] &&
       awk -F '\t' -v expected="$EXPECTED_HAPLOTYPE_COLUMNS" 'NR==1{exit NF==expected?0:1}' "$haplotype_shard"; then
        log_info "断点续跑：跳过已完成的单倍型染色体$chrom。"
        continue
    fi
    haplotype_tmp="$HAPLOTYPE_SHARD_DIR/.$shard_name.tmp.tsv"
    log_info "统计染色体$chrom的$chrom_gene_count个$VARIANTION_SOURCE基因单倍型。"
    bcftools query -R "$chrom_regions" -s "$HAPLOTYPE_CSV" -f "$QUERY_FORMAT" "$INPUT_VCF" | "$COUNTER_BIN" haplotype --chrom "$chrom" --gene-bed "$chrom_gene_bed" --samples "$HAPLOTYPE_SAMPLES" --sample-group "$SAMPLE_GROUP_TABLE" --group-order "$SAMPLE_GROUP_ORDER" --top-number "$HAPLOTYPE_TOP_NUMBER" --output "$haplotype_tmp" || die "染色体$chrom单倍型统计失败。"
    [[ $(awk 'END{print NR+0}' "$haplotype_tmp") -eq $((chrom_gene_count+1)) ]] || die "染色体$chrom单倍型结果行数异常。"
    awk -F '\t' -v expected="$EXPECTED_HAPLOTYPE_COLUMNS" 'NR==1{exit NF==expected?0:1}' "$haplotype_tmp" || die "染色体$chrom单倍型结果列数异常。"
    mv -f -- "$haplotype_tmp" "$haplotype_shard"
done < "$GENE_CHROMS"

HAPLOTYPE_OUTPUT_TMP="$HAPLOTYPE_OUTPUT.tmp"
: > "$HAPLOTYPE_OUTPUT_TMP"
for ((index=1;index<=CHROM_COUNT;index++)); do
    haplotype_shard="$HAPLOTYPE_SHARD_DIR/$(printf '%06d' "$index").tsv"
    if (( index == 1 )); then cat "$haplotype_shard" >> "$HAPLOTYPE_OUTPUT_TMP"; else tail -n +2 "$haplotype_shard" >> "$HAPLOTYPE_OUTPUT_TMP"; fi
done
[[ $(awk 'END{print NR+0}' "$HAPLOTYPE_OUTPUT_TMP") -eq $((GENE_COUNT+1)) ]] || die "合并后单倍型统计表行数异常。"
awk -F '\t' -v expected="$EXPECTED_HAPLOTYPE_COLUMNS" 'NR==1{exit NF==expected?0:1}' "$HAPLOTYPE_OUTPUT_TMP" || die "合并后单倍型统计表列数异常。"
mv -f -- "$HAPLOTYPE_OUTPUT_TMP" "$HAPLOTYPE_OUTPUT"
cp -- "$HAPLOTYPE_REQUESTED_CONFIG" "$HAPLOTYPE_CONFIG.tmp";mv -f -- "$HAPLOTYPE_CONFIG.tmp" "$HAPLOTYPE_CONFIG"
{
    printf 'FORMAT\tggComp.plus.haplotype-statistics.complete\n'
    printf 'VARIANTION_SOURCE\t%s\n' "$VARIANTION_SOURCE"
    printf 'STAGE\t%s\n' "$HAPLOTYPE_STAGE"
    printf 'KEY\t%s\n' "$HAPLOTYPE_KEY"
    printf 'GENE_COUNT\t%s\n' "$GENE_COUNT"
    printf 'TOP_NUMBER\t%s\n' "$HAPLOTYPE_TOP_NUMBER"
    printf 'OUTPUT\t%s\n' "$(basename "$HAPLOTYPE_OUTPUT")"
} > "$HAPLOTYPE_MARKER.tmp"
mv -f -- "$HAPLOTYPE_MARKER.tmp" "$HAPLOTYPE_MARKER"
record_checkpoint "$HAPLOTYPE_STAGE" "$HAPLOTYPE_KEY" "$HAPLOTYPE_MARKER" || die "无法写入单倍型断点记录。"
deduplicate_checkpoints || die "整理单倍型断点失败。"
rm -rf -- "$HAPLOTYPE_SHARD_DIR"
log_info "$VARIANTION_SOURCE单倍型统计完成: $HAPLOTYPE_OUTPUT"
}

run_normal_group_comparison
if [[ "$HAPLOTYPE" == true ]]; then run_haplotype_analysis; fi

log_info "$VARIANTION_SOURCE来源全部分析完成。"
log_info "$VARIANTION_SOURCE基因二进制缓存: $GENE_DIR"
log_info "Group比较结果目录: $OUTPUT_DIR/GeneDSR"
if [[ "$HAPLOTYPE" == true ]]; then log_info "单倍型结果目录: $OUTPUT_DIR/GeneHaplotype"; fi
log_info "日志文件: $LOG_FILE"
