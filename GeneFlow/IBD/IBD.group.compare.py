#!/usr/bin/env python3

import os
import re
import csv
import gzip
import argparse
import statistics
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "根据样本分组信息，对 hap-ibd 结果进行群体间 IBD 统计。\n"
            "同一样本对、同一染色体上重叠的 IBD 区间会先进行合并，\n"
            "然后计算每个样本对的 IBD 数量和共享长度。\n"
            "没有检测到 IBD 的样本对也会输出，并以 0 计入群体统计。"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--input-ibd",
        required=True,
        help="hap-ibd 输出的 .ibd 或 .ibd.gz 文件"
    )

    parser.add_argument(
        "--input-sample-group",
        required=True,
        help=(
            "样本分组文件，第一行为表头，必须包含 SampleID 和 Group 两列。\n"
            "支持逗号、Tab或空白字符分隔。"
        )
    )

    parser.add_argument(
        "--input-compartion",
        "--input-comparison",
        dest="input_comparison",
        required=True,
        help=(
            "需要比较的群体组合文件，无表头，每行两列，例如：\n"
            "XianLandrace1 XianLandrace2\n"
            "XianLandrace1 XianCultivar1"
        )
    )

    parser.add_argument(
        "--output-dir",
        default=".",
        help="输出目录，默认当前目录 [default: .]"
    )

    parser.add_argument(
        "--bp-per-cm",
        type=float,
        default=1000000,
        help=(
            "每 1 cM 对应的物理距离(bp)。\n"
            "你当前使用 1 cM/Mb 的人工遗传图谱，因此默认 1000000。\n"
            "[default: 1000000]"
        )
    )

    return parser.parse_args()


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def split_line(line):
    line = line.strip()

    if "," in line:
        return [x.strip() for x in line.split(",")]

    return re.split(r"\s+", line)


def read_sample_group(path):
    sample_to_group = {}
    group_to_samples = defaultdict(list)

    with open_text(path) as f:

        header = None

        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            header = split_line(line)
            break

        if header is None:
            raise ValueError("样本分组文件为空。")

        if "SampleID" not in header:
            raise ValueError("样本分组文件中没有找到 SampleID 列。")

        if "Group" not in header:
            raise ValueError("样本分组文件中没有找到 Group 列。")

        sample_idx = header.index("SampleID")
        group_idx = header.index("Group")

        for line_number, line in enumerate(f, start=2):

            if not line.strip() or line.startswith("#"):
                continue

            fields = split_line(line)

            if len(fields) <= max(sample_idx, group_idx):
                raise ValueError(
                    f"样本分组文件第 {line_number} 行列数不足："
                    f"{line.strip()}"
                )

            sample = fields[sample_idx]
            group = fields[group_idx]

            if sample in sample_to_group:
                raise ValueError(
                    f"样本 {sample} 在分组文件中重复出现。"
                )

            sample_to_group[sample] = group
            group_to_samples[group].append(sample)

    return sample_to_group, group_to_samples


def read_comparisons(path):
    comparisons = []

    with open_text(path) as f:

        for line_number, line in enumerate(f, start=1):

            if not line.strip() or line.startswith("#"):
                continue

            fields = split_line(line)

            if len(fields) < 2:
                raise ValueError(
                    f"比较文件第 {line_number} 行不足两列："
                    f"{line.strip()}"
                )

            comparisons.append((fields[0], fields[1]))

    if not comparisons:
        raise ValueError("比较文件中没有找到有效的群体组合。")

    return comparisons


def canonical_pair(sample1, sample2):
    if sample1 <= sample2:
        return sample1, sample2
    return sample2, sample1


def build_wanted_pairs(comparisons, group_to_samples):
    wanted_pairs = set()

    for group1, group2 in comparisons:

        if group1 not in group_to_samples:
            raise ValueError(
                f"群体 {group1} 不存在于样本分组文件。"
            )

        if group2 not in group_to_samples:
            raise ValueError(
                f"群体 {group2} 不存在于样本分组文件。"
            )

        samples1 = group_to_samples[group1]
        samples2 = group_to_samples[group2]

        if group1 == group2:

            for i in range(len(samples1)):
                for j in range(i + 1, len(samples1)):
                    wanted_pairs.add(
                        canonical_pair(
                            samples1[i],
                            samples1[j]
                        )
                    )

        else:

            for sample1 in samples1:
                for sample2 in samples2:
                    wanted_pairs.add(
                        canonical_pair(
                            sample1,
                            sample2
                        )
                    )

    return wanted_pairs


def read_ibd(path, wanted_pairs):
    """
    保存格式：
        pair -> chromosome -> [(start, end), ...]

    hap-ibd:
        col1 Sample1
        col3 Sample2
        col5 Chr
        col6 Start
        col7 End
        col8 cM

    这里不直接使用第8列进行累加，
    而是在区间合并后依据 bp-per-cm 重新计算长度。
    """

    ibd_intervals = defaultdict(
        lambda: defaultdict(list)
    )

    total_lines = 0
    used_lines = 0
    malformed_lines = 0

    with open_text(path) as f:

        for line in f:

            if not line.strip() or line.startswith("#"):
                continue

            total_lines += 1

            fields = re.split(r"\s+", line.strip())

            if len(fields) < 8:
                malformed_lines += 1
                continue

            sample1 = fields[0]
            sample2 = fields[2]
            chromosome = fields[4]

            try:
                start = int(fields[5])
                end = int(fields[6])
            except ValueError:
                malformed_lines += 1
                continue

            if end < start:
                start, end = end, start

            pair = canonical_pair(
                sample1,
                sample2
            )

            if pair not in wanted_pairs:
                continue

            ibd_intervals[pair][chromosome].append(
                (start, end)
            )

            used_lines += 1

    return (
        ibd_intervals,
        total_lines,
        used_lines,
        malformed_lines
    )


def merge_intervals(intervals):
    """
    合并重叠区间。

    例如：
        100-500
        300-700
        800-1000

    合并为：
        100-700
        800-1000
    """

    if not intervals:
        return []

    intervals = sorted(
        intervals,
        key=lambda x: (x[0], x[1])
    )

    merged = []

    current_start, current_end = intervals[0]

    for start, end in intervals[1:]:

        if start <= current_end:

            if end > current_end:
                current_end = end

        else:

            merged.append(
                (current_start, current_end)
            )

            current_start = start
            current_end = end

    merged.append(
        (current_start, current_end)
    )

    return merged


def get_merged_lengths(pair_intervals, bp_per_cm):
    """
    对一个 sample pair：
      1. 每条染色体分别合并
      2. 将 merged physical interval 转换为 cM

    当前：
        1 cM = bp_per_cm bp
    """

    lengths_cm = []

    for chromosome in pair_intervals:

        intervals = pair_intervals[chromosome]

        merged = merge_intervals(intervals)

        for start, end in merged:

            length_bp = end - start

            length_cm = (
                length_bp / bp_per_cm
            )

            lengths_cm.append(length_cm)

    return lengths_cm


def calculate_stats(
    pair_intervals,
    bp_per_cm
):

    if not pair_intervals:
        return (
            0,
            0.0,
            0.0,
            0.0,
            0.0
        )

    lengths = get_merged_lengths(
        pair_intervals,
        bp_per_cm
    )

    if not lengths:
        return (
            0,
            0.0,
            0.0,
            0.0,
            0.0
        )

    n_ibd = len(lengths)

    total_ibd = sum(lengths)

    mean_ibd = (
        total_ibd / n_ibd
    )

    median_ibd = statistics.median(
        lengths
    )

    max_ibd = max(lengths)

    return (
        n_ibd,
        total_ibd,
        mean_ibd,
        median_ibd,
        max_ibd
    )


def safe_filename(text):
    return re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        text
    )


def get_sample_pairs(
    group1,
    group2,
    group_to_samples
):

    samples1 = sorted(
        group_to_samples[group1]
    )

    samples2 = sorted(
        group_to_samples[group2]
    )

    if group1 == group2:

        for i in range(len(samples1)):
            for j in range(i + 1, len(samples1)):

                yield (
                    samples1[i],
                    samples1[j]
                )

    else:

        for sample1 in samples1:
            for sample2 in samples2:

                yield (
                    sample1,
                    sample2
                )


def write_group_comparison(
    group1,
    group2,
    group_to_samples,
    ibd_intervals,
    output_dir,
    bp_per_cm
):

    filename = (
        f"{safe_filename(group1)}."
        f"vs."
        f"{safe_filename(group2)}."
        f"stats.csv"
    )

    output_path = os.path.join(
        output_dir,
        filename
    )

    total_pairs = 0
    pairs_with_ibd = 0

    total_ibd_values = []
    n_ibd_values = []

    with open(
        output_path,
        "w",
        newline="",
        encoding="utf-8"
    ) as out:

        writer = csv.writer(out)

        writer.writerow([
            "Sample1",
            "Sample2",
            "Group1",
            "Group2",
            "N_IBD",
            "Total_IBD_cM",
            "Mean_IBD_cM",
            "Median_IBD_cM",
            "Max_IBD_cM"
        ])

        for sample1, sample2 in get_sample_pairs(
            group1,
            group2,
            group_to_samples
        ):

            pair = canonical_pair(
                sample1,
                sample2
            )

            pair_intervals = ibd_intervals.get(
                pair,
                {}
            )

            stats = calculate_stats(
                pair_intervals,
                bp_per_cm
            )

            (
                n_ibd,
                total_ibd,
                mean_ibd,
                median_ibd,
                max_ibd
            ) = stats

            if n_ibd > 0:
                pairs_with_ibd += 1

            total_pairs += 1

            total_ibd_values.append(
                total_ibd
            )

            n_ibd_values.append(
                n_ibd
            )

            writer.writerow([
                sample1,
                sample2,
                group1,
                group2,
                n_ibd,
                f"{total_ibd:.6f}",
                f"{mean_ibd:.6f}",
                f"{median_ibd:.6f}",
                f"{max_ibd:.6f}"
            ])

    if total_pairs > 0:

        fraction_with_ibd = (
            pairs_with_ibd / total_pairs
        )

        mean_total_ibd = statistics.mean(
            total_ibd_values
        )

        median_total_ibd = statistics.median(
            total_ibd_values
        )

        mean_n_ibd = statistics.mean(
            n_ibd_values
        )

    else:

        fraction_with_ibd = 0.0
        mean_total_ibd = 0.0
        median_total_ibd = 0.0
        mean_n_ibd = 0.0

    summary = {
        "Group1": group1,
        "Group2": group2,
        "N_Group1": len(
            group_to_samples[group1]
        ),
        "N_Group2": len(
            group_to_samples[group2]
        ),
        "N_SamplePairs": total_pairs,
        "N_Pairs_with_IBD": pairs_with_ibd,
        "Fraction_Pairs_with_IBD": fraction_with_ibd,
        "Mean_N_IBD_per_sample_pair": mean_n_ibd,
        "Mean_Total_IBD_cM_per_sample_pair": mean_total_ibd,
        "Median_Total_IBD_cM_per_sample_pair": median_total_ibd
    }

    return (
        output_path,
        summary
    )


def write_group_summary(
    summaries,
    output_dir
):

    output_path = os.path.join(
        output_dir,
        "IBD.group.comparison.summary.csv"
    )

    fields = [
        "Group1",
        "Group2",
        "N_Group1",
        "N_Group2",
        "N_SamplePairs",
        "N_Pairs_with_IBD",
        "Fraction_Pairs_with_IBD",
        "Mean_N_IBD_per_sample_pair",
        "Mean_Total_IBD_cM_per_sample_pair",
        "Median_Total_IBD_cM_per_sample_pair"
    ]

    with open(
        output_path,
        "w",
        newline="",
        encoding="utf-8"
    ) as out:

        writer = csv.DictWriter(
            out,
            fieldnames=fields
        )

        writer.writeheader()

        for summary in summaries:

            row = summary.copy()

            row["Fraction_Pairs_with_IBD"] = (
                f"{row['Fraction_Pairs_with_IBD']:.6f}"
            )

            row["Mean_N_IBD_per_sample_pair"] = (
                f"{row['Mean_N_IBD_per_sample_pair']:.6f}"
            )

            row["Mean_Total_IBD_cM_per_sample_pair"] = (
                f"{row['Mean_Total_IBD_cM_per_sample_pair']:.6f}"
            )

            row["Median_Total_IBD_cM_per_sample_pair"] = (
                f"{row['Median_Total_IBD_cM_per_sample_pair']:.6f}"
            )

            writer.writerow(row)

    return output_path


def main():

    args = parse_args()

    os.makedirs(
        args.output_dir,
        exist_ok=True
    )

    print("[INFO] 读取样本分组文件……")

    (
        sample_to_group,
        group_to_samples
    ) = read_sample_group(
        args.input_sample_group
    )

    print(
        f"[INFO] 共读取 "
        f"{len(sample_to_group)} 个样本，"
        f"{len(group_to_samples)} 个群体。"
    )

    for group in sorted(group_to_samples):

        print(
            f"[INFO] Group {group}: "
            f"{len(group_to_samples[group])} samples"
        )

    print("[INFO] 读取群体比较文件……")

    comparisons = read_comparisons(
        args.input_comparison
    )

    print(
        f"[INFO] 共读取 "
        f"{len(comparisons)} 个群体比较。"
    )

    print(
        f"[INFO] 遗传距离转换："
        f"1 cM = {args.bp_per_cm:g} bp"
    )

    print("[INFO] 建立目标样本对……")

    wanted_pairs = build_wanted_pairs(
        comparisons,
        group_to_samples
    )

    print(
        f"[INFO] 共需要分析 "
        f"{len(wanted_pairs)} 个唯一样本对。"
    )

    print("[INFO] 读取 hap-ibd 文件……")

    (
        ibd_intervals,
        total_lines,
        used_lines,
        malformed_lines
    ) = read_ibd(
        args.input_ibd,
        wanted_pairs
    )

    print(
        f"[INFO] IBD 文件数据行："
        f"{total_lines}"
    )

    print(
        f"[INFO] 属于目标样本对的 IBD 记录："
        f"{used_lines}"
    )

    print(
        f"[INFO] 至少具有一条 IBD 的样本对："
        f"{len(ibd_intervals)}"
    )

    if malformed_lines > 0:

        print(
            f"[WARNING] 格式异常记录："
            f"{malformed_lines}"
        )

    summaries = []

    print(
        "[INFO] 合并重叠 IBD 区间并输出结果……"
    )

    for group1, group2 in comparisons:

        output_path, summary = (
            write_group_comparison(
                group1,
                group2,
                group_to_samples,
                ibd_intervals,
                args.output_dir,
                args.bp_per_cm
            )
        )

        summaries.append(summary)

        print(
            f"[INFO] {group1} vs {group2}: "
            f"{summary['N_SamplePairs']} pairs; "
            f"{summary['N_Pairs_with_IBD']} "
            f"pairs with IBD; "
            f"Mean total IBD = "
            f"{summary['Mean_Total_IBD_cM_per_sample_pair']:.6f} cM/pair"
        )

        print(
            f"[INFO] 输出："
            f"{output_path}"
        )

    summary_path = write_group_summary(
        summaries,
        args.output_dir
    )

    print(
        f"[INFO] 群体比较汇总："
        f"{summary_path}"
    )

    print("[INFO] 全部处理完成。")


if __name__ == "__main__":
    main()