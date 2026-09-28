#!/usr/bin/env python3

import argparse
import errno
import multiprocessing
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass, fields
from itertools import zip_longest

try:
    import pysam
except ImportError:
    sys.stderr.write(
        "[ERROR] 未检测到 pysam 模块。\n"
        "请先安装：conda install -c bioconda pysam\n"
    )
    sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(
        prog="restore_missing.py",
        description=(
            "保留 Beagle 的定相 GT，同时完整恢复原始 VCF 的缺失 GT、Header、"
            "QUAL、FILTER、INFO 以及其他 FORMAT 信息。\n\n"
            "处理原则：\n"
            "  1. 输出 VCF 完全以原始 VCF 为基础；\n"
            "  2. 所有 ## Header 注释均来自原始 VCF；\n"
            "  3. 样本 ID 及样本顺序均使用原始 VCF；\n"
            "  4. CHROM、POS、ID、REF、ALT、QUAL、FILTER、INFO 均来自原始 VCF；\n"
            "  5. 原始 VCF 中 GT 缺失的样本保持原始缺失状态；\n"
            "  6. 原始 VCF 中 GT 非缺失的样本使用 Beagle 定相后的 GT；\n"
            "  7. DP、AD、GQ 等其他 FORMAT 信息保持原始 VCF 内容；\n"
            "  8. 不保留 Beagle 新增的 INFO 或 FORMAT 信息。"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--input-original-vcf",
        required=True,
        metavar="FILE",
        help=(
            "原始 VCF 文件。\n"
            "支持 .vcf、.vcf.gz 和 .bcf。"
        )
    )

    parser.add_argument(
        "--input-beagle-vcf",
        required=True,
        metavar="FILE",
        help=(
            "Beagle 定相后的 VCF 文件。\n"
            "仅使用其中的 GT 和 phase 信息。"
        )
    )

    parser.add_argument(
        "--out-keep-missing",
        required=True,
        metavar="FILE",
        help=(
            "输出文件。\n"
            "例如：sample.beagle.keepMissing.vcf.gz"
        )
    )

    parser.add_argument(
        "--parallel", type=str.lower, choices=("true", "false"), default="true",
        help=(
            "是否按染色体并行，默认 true。\n"
            "true：自动为每条有记录的染色体启动一个进程，每个文件使用 1 个 I/O 线程。\n"
            "输入须为已建 .tbi/.csi 索引的 BGZF VCF 或 BCF；输出保持原始染色体顺序。\n"
            "false：顺序读取，无需输入索引。\n"
            "染色体临时文件存放在当前工作目录下，结束后自动清理。"
        )
    )
    args = parser.parse_args()
    args.parallel = args.parallel == "true"
    return args


def variant_key(record):
    """
    返回用于检查两个 VCF 位点是否一致的键。
    """
    return (
        record.contig,
        record.pos,
        record.ref,
        tuple(record.alts or ())
    )


def is_missing_gt(gt):
    """
    判断 GT 是否包含缺失 allele。

    以下均认为存在缺失：
        ./.
        0/.
        ./1
        .
    """
    if gt is None:
        return True

    return any(allele is None for allele in gt)


def output_mode(filename):
    """
    根据文件后缀选择 pysam 输出模式。
    """
    if filename.endswith(".vcf.gz") or filename.endswith(".vcf.bgz"):
        return "wz"

    if filename.endswith(".bcf"):
        return "wb"

    return "w"


@dataclass
class Counts:
    variants: int = 0
    sample_gts: int = 0
    missing: int = 0
    phased: int = 0
    original_hets: int = 0
    phased_hets: int = 0

    def add(self, other):
        for field in fields(self):
            setattr(self, field.name,
                    getattr(self, field.name) + getattr(other, field.name))


def log(message):
    sys.stderr.write(f"[INFO] {message}\n")
    sys.stderr.flush()


def check_samples(original_vcf, beagle_vcf):
    original_samples = list(original_vcf.header.samples)
    beagle_samples = list(beagle_vcf.header.samples)
    log(f"原始 VCF 样本数：{len(original_samples):,}")
    log(f"Beagle VCF 样本数：{len(beagle_samples):,}")
    if set(original_samples) != set(beagle_samples):
        details = ["原始 VCF 与 Beagle VCF 的样本集合不一致。"]
        for label, samples in (
            ("原始 VCF", set(original_samples) - set(beagle_samples)),
            ("Beagle VCF", set(beagle_samples) - set(original_samples)),
        ):
            if samples:
                details.append(f"仅存在于{label}的样本："
                               + ",".join(sorted(samples)[:20]))
        raise ValueError("\n".join(details))
    if original_samples != beagle_samples:
        log("样本顺序不同，将按 Sample ID 匹配 GT，输出使用原始样本顺序。")
    return original_samples


def process_records(original_records, beagle_records, output_vcf, samples,
                    label):
    """顺序和并行共用同一套 GT 恢复逻辑，只修改原始记录中的非缺失 GT。"""
    counts = Counts()
    for number, (original, beagle) in enumerate(
        zip_longest(original_records, beagle_records), start=1
    ):
        if original is None or beagle is None:
            relation = "多于" if original is None else "少于"
            raise RuntimeError(f"{label}：Beagle VCF 的位点数量{relation}原始 VCF。")
        if variant_key(original) != variant_key(beagle):
            raise RuntimeError(
                f"{label}：检测到两个 VCF 的位点不一致；记录编号：{number:,}\n"
                f"原始 VCF：{variant_key(original)}\n"
                f"Beagle VCF：{variant_key(beagle)}"
            )

        # 记录用完即写出，无需复制整个记录；缓存样本映射减少属性访问。
        original_calls = original.samples
        beagle_calls = beagle.samples
        for sample in samples:
            original_call = original_calls[sample]
            original_gt = original_call.get("GT")
            if is_missing_gt(original_gt):
                counts.missing += 1
                continue

            # 原始缺失 GT 无需读取对应的 Beagle GT。
            beagle_call = beagle_calls[sample]
            beagle_gt = beagle_call.get("GT")
            if is_missing_gt(beagle_gt):
                raise RuntimeError(
                    "检测到原始 GT 非缺失，但 Beagle GT 缺失：\n"
                    f"位点：{original.contig}:{original.pos}；样本：{sample}\n"
                    f"原始 GT：{original_gt}；Beagle GT：{beagle_gt}"
                )
            if len(original_gt) == 2 and original_gt[0] != original_gt[1]:
                counts.original_hets += 1
            phased = beagle_call.phased
            original_call["GT"] = beagle_gt
            original_call.phased = phased
            if phased:
                counts.phased += 1
                if len(beagle_gt) == 2 and beagle_gt[0] != beagle_gt[1]:
                    counts.phased_hets += 1

        output_vcf.write(original)
        counts.variants += 1
        if counts.variants % 100000 == 0:
            log(f"{label}：已处理 {counts.variants:,} 个位点；"
                f"保留 {counts.missing:,} 个原始缺失 GT。")
    counts.sample_gts = counts.variants * len(samples)
    return counts


def indexed_contigs(vcf, filename):
    """仅探测每条染色体的首条记录，按文件偏移确定原始排列顺序。"""
    if vcf.index is None:
        raise ValueError(
            f"并行模式需要可用索引：{filename}\n"
            "请先用 bgzip 压缩普通 VCF，并用 bcftools index 建立索引；"
            "或者设置 --parallel false 顺序处理。"
        )
    locations = []
    # BCF 索引可能包含 Header 中没有实际记录的染色体，需要跳过。
    for contig in vcf.index:
        if next(vcf.fetch(contig), None) is not None:
            locations.append((vcf.tell(), contig))
    return [contig for _, contig in sorted(locations)]


def process_chromosome(task):
    """子进程自行打开文件，避免共享 HTSlib 句柄及传输大批记录。"""
    original_file, beagle_file, part_file, contig = task
    try:
        with pysam.VariantFile(original_file, threads=1) as original, \
             pysam.VariantFile(beagle_file, threads=1) as beagle, \
             pysam.VariantFile(part_file, "wb", header=original.header,
                               threads=1) as output:
            counts = process_records(
                original.fetch(contig), beagle.fetch(contig), output,
                list(original.header.samples), contig
            )
        return contig, counts
    except Exception as exc:
        raise RuntimeError(f"染色体 {contig} 处理失败：{exc}") from exc


def run_parallel(args, contigs, header, stage_file):
    counts = Counts()
    if not contigs:
        with pysam.VariantFile(stage_file, output_mode(args.out_keep_missing),
                               header=header, threads=1):
            pass
        return counts

    workers = len(contigs)
    with tempfile.TemporaryDirectory(prefix="restore_missing.parts.",
                                     dir=os.getcwd()) as tmp_dir:
        tasks = [
            (args.input_original_vcf, args.input_beagle_vcf,
             os.path.join(tmp_dir, f"{index:06d}.bcf"), contig)
            for index, contig in enumerate(contigs)
        ]
        log(f"按染色体并行：{len(contigs)} 条染色体，{workers} 个进程；"
            "每个文件 1 个 I/O 线程。")
        log(f"临时 BCF 目录：{tmp_dir}")
        # spawn 避免继承父进程的 HTSlib 文件句柄和压缩线程。
        with multiprocessing.get_context("spawn").Pool(workers) as pool:
            for completed, (contig, result) in enumerate(
                pool.imap_unordered(process_chromosome, tasks, chunksize=1), 1
            ):
                counts.add(result)
                log(f"染色体 {contig} 完成：{result.variants:,} 个位点；"
                    f"进度 {completed}/{len(contigs)}。")

        log("开始按原始文件的染色体顺序合并结果...")
        merged = 0
        with pysam.VariantFile(stage_file, output_mode(args.out_keep_missing),
                               header=header, threads=1) as output:
            for task in tasks:
                # 临时 BCF 仅顺序读取，不需要索引。
                # 使用文件对象，避免 HTSlib 尝试加载不需要的索引并打印错误信息。
                with open(task[2], "rb") as stream, \
                     pysam.VariantFile(stream, threads=1) as part:
                    for record in part:
                        record.translate(output.header)
                        output.write(record)
                        merged += 1
                os.unlink(task[2])
        if merged != counts.variants:
            raise RuntimeError("合并位点数与分染色体统计不一致。")
    return counts


def validate_paths(args):
    for filename in (args.input_original_vcf, args.input_beagle_vcf):
        if not os.path.isfile(filename):
            raise ValueError(f"输入文件不存在：{filename}")
        if (os.path.realpath(filename) == os.path.realpath(args.out_keep_missing)
                or (os.path.exists(args.out_keep_missing)
                    and os.path.samefile(filename, args.out_keep_missing))):
            raise ValueError("输出文件不能与任一输入文件相同。")
    for suffix in (".tbi", ".csi"):
        if os.path.exists(args.out_keep_missing + suffix):
            raise ValueError(
                f"输出路径已有索引：{args.out_keep_missing + suffix}；"
                "请使用新的输出路径，避免旧索引与新结果不匹配。"
            )


def publish_output(stage_file, output_file):
    """同一文件系统直接替换；跨文件系统先复制到目标目录再原子替换。"""
    try:
        os.replace(stage_file, output_file)
        return
    except OSError as exc:
        if exc.errno != errno.EXDEV:
            raise

    # 跨文件系统无法原子移动；仅发布时在目标目录短暂暂存，避免覆盖半成品。
    fd, destination_stage = tempfile.mkstemp(
        prefix=".restore_missing.publish.", suffix=".tmp",
        dir=os.path.dirname(os.path.abspath(output_file))
    )
    os.close(fd)
    try:
        shutil.copyfile(stage_file, destination_stage)
        os.replace(destination_stage, output_file)
    finally:
        if os.path.exists(destination_stage):
            os.unlink(destination_stage)
    os.unlink(stage_file)


def main():
    args = parse_args()
    parallel = args.parallel
    stage_file = None
    try:
        validate_paths(args)
        log(f"原始 VCF：{args.input_original_vcf}")
        log(f"Beagle VCF：{args.input_beagle_vcf}")
        log(f"输出文件：{args.out_keep_missing}")
        with pysam.VariantFile(args.input_original_vcf,
                               threads=1) as original, \
             pysam.VariantFile(args.input_beagle_vcf,
                               threads=1) as beagle:
            samples = check_samples(original, beagle)
            header = original.header.copy()
            contigs = None
            if parallel:
                contigs = indexed_contigs(original, args.input_original_vcf)
                beagle_contigs = indexed_contigs(beagle, args.input_beagle_vcf)
                if set(contigs) != set(beagle_contigs):
                    raise ValueError(
                        "两个 VCF 中有记录的染色体集合不一致。\n"
                        f"仅原始 VCF 有：{sorted(set(contigs) - set(beagle_contigs))}\n"
                        f"仅 Beagle VCF 有：{sorted(set(beagle_contigs) - set(contigs))}"
                    )

            # 在当前工作目录暂存结果，全部处理成功后再发布到输出路径。
            fd, stage_file = tempfile.mkstemp(
                prefix=".restore_missing.", suffix=".tmp",
                dir=os.getcwd()
            )
            os.close(fd)
            if not parallel:
                log("开始顺序处理 VCF...")
                with pysam.VariantFile(
                    stage_file, output_mode(args.out_keep_missing),
                    header=header, threads=1
                ) as output:
                    counts = process_records(original, beagle, output, samples,
                                             "全部位点")

        if parallel:
            counts = run_parallel(args, contigs, header, stage_file)
        publish_output(stage_file, args.out_keep_missing)
        stage_file = None
        log("处理完成。")
        for label, value in (
            ("总位点数", counts.variants),
            ("总样本 GT 数", counts.sample_gts),
            ("保留的原始缺失 GT 数", counts.missing),
            ("使用 Beagle phased GT 数", counts.phased),
            ("原始非缺失杂合 GT 数", counts.original_hets),
            ("Beagle 已定相杂合 GT 数", counts.phased_hets),
        ):
            log(f"{label}：{value:,}")
        log(f"输出文件：{args.out_keep_missing}")
        return 0
    except KeyboardInterrupt:
        sys.stderr.write("[ERROR] 用户中断，未生成最终输出。\n")
        return 130
    except Exception as exc:
        sys.stderr.write(f"[ERROR] VCF 处理失败：\n{exc}\n")
        return 1
    finally:
        if stage_file is not None:
            os.unlink(stage_file)


if __name__ == "__main__":
    sys.exit(main())
