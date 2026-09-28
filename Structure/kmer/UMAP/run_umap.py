#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
使用 GCTA PCA 结果中的主成分运行 UMAP，并绘制 PCA 碎石图。

输入前缀对应以下两个 GCTA 输出文件：
    <prefix>.eigenvec
    <prefix>.eigenval

UMAP 输出为制表符分隔文件：
    FID    IID    UMAP1    UMAP2    ...
"""

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


class CustomHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    """保留描述文本中的换行，并显示参数默认值。"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "1.使用 GCTA PCA 结果的前 N 个主成分运行 UMAP，并通过 R 绘制 PCA 碎石图;\n"
            "2.输入前缀对应以下两个 GCTA 输出文件：<prefix>.eigenvec和<prefix>.eigenval;\n"
            "3.UMAP 坐标输出为制表符分隔文件，并根据 Excel 样本分组绘制 UMAP 图。\n"
            "\n使用示例:\n"
            "python3 run_umap.py --pca-result XianSample.OnlySNP.LDfilter.PCA "
            "--input-sample-group SampleGroup.xlsx\n"
            "python3 run_umap.py --pca-result XianSample.OnlySNP.LDfilter.PCA "
            "--input-sample-group SampleGroup.xlsx --n-pc 10 "
            "--n-neighbors 30 --min-dist 0.3\n"
        ),
        formatter_class=CustomHelpFormatter,
        add_help=False,
    )

    parser.add_argument(
        "--pca-result",
        required=True,
        metavar="PREFIX",
        help="GCTA PCA 结果文件前缀；脚本将读取 PREFIX.eigenvec 和 PREFIX.eigenval。",
    )

    parser.add_argument(
        "--input-sample-group",
        required=True,
        metavar="XLSX",
        help=(
            "Excel 样本分组文件（.xlsx 或 .xlsm）；第一行表头必须为 "
            "SampleID 和 Group，SampleID 与 eigenvec 的 IID 匹配。"
        ),
    )

    parser.add_argument(
        "--output-umap",
        default=argparse.SUPPRESS,
        metavar="FILE",
        help=(
            "UMAP 坐标输出文件（TSV 格式）；省略时按 PCA 参数自动命名为 "
            "PCA_PREFIX.pcN.nneighborsN.mindistN.UMAP.list。"
        ),
    )

    parser.add_argument(
        "--output-umap-plot",
        default=argparse.SUPPRESS,
        metavar="PREFIX",
        help=(
            "UMAP 图输出前缀；脚本同时生成 PREFIX.pdf 和 PREFIX.png。"
            "省略时按 PCA 参数自动命名为 PCA_PREFIX.pcN.nneighborsN."
            "mindistN.UMAP.plot。"
        ),
    )

    parser.add_argument(
        "--color-mapping-group",
        metavar="FILE",
        help=(
            "可选的群体颜色映射文件，无表头、两列，依次为群体 ID 和 "
            "#RRGGBB 颜色代码；省略时按群体首次出现顺序使用内置颜色列表。"
        ),
    )

    parser.add_argument(
        "--n-pc",
        type=int,
        default=10,
        help="用于 UMAP 输入的前置主成分数量。",
    )

    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=30,
        help="构建 UMAP 最近邻图时使用的邻居样本数。",
    )

    parser.add_argument(
        "--min-dist",
        type=float,
        default=0.2,
        help="UMAP 嵌入中样本点之间的最小距离。",
    )

    parser.add_argument(
        "--metric",
        default="euclidean",
        help="在 PCA 空间中使用的距离度量。",
    )

    parser.add_argument(
        "--n-components",
        dest="n_components",
        type=int,
        default=2,
        help="输出的 UMAP 轴（维度）数量。",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="用于复现 UMAP 结果的随机种子。",
    )

    parser.add_argument(
        "--n-epochs",
        type=int,
        default=None,
        help="UMAP 优化轮数；省略时由 UMAP 自动选择。",
    )

    parser.add_argument(
        "--init",
        choices=["spectral", "random", "pca", "tswspectral"],
        default="spectral",
        help="UMAP 嵌入的初始化方法。",
    )

    parser.add_argument(
        "--scree-plot",
        metavar="FILE",
        help=(
            "PCA 碎石图输出路径，支持 .pdf 或 .png；省略时输出为 "
            "PCA_PREFIX.scree_plot.pdf。图中固定展示全部 PC。"
        ),
    )

    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="显示此帮助文档并退出。",
    )

    args = parser.parse_args()
    output_prefix = (
        f"{args.pca_result}.pc{args.n_pc}"
        f".nneighbors{args.n_neighbors}"
        f".mindist{args.min_dist}"
    )
    if not hasattr(args, "output_umap"):
        args.output_umap = f"{output_prefix}.UMAP.list"
    if not hasattr(args, "output_umap_plot"):
        args.output_umap_plot = f"{output_prefix}.UMAP.plot"
    return args


def get_gcta_paths(prefix: str) -> tuple[Path, Path]:
    """根据 GCTA 结果前缀生成 eigenvec 和 eigenval 文件路径。"""
    return Path(f"{prefix}.eigenvec"), Path(f"{prefix}.eigenval")


def get_umap_plot_paths(output_prefix: str) -> tuple[Path, Path]:
    """根据 UMAP 图输出前缀生成 PDF 和 PNG 路径。"""
    plot_prefix = Path(output_prefix)
    if plot_prefix.suffix.lower() in {".pdf", ".png"}:
        plot_prefix = plot_prefix.with_suffix("")
    return Path(f"{plot_prefix}.pdf"), Path(f"{plot_prefix}.png")


def run_r_code(
    rscript: str,
    r_code: str,
    arguments: list[str],
) -> subprocess.CompletedProcess[str]:
    """通过临时 R 脚本运行代码，避免命令行长度限制。"""
    script_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            suffix=".R",
            delete=False,
        ) as handle:
            handle.write(r_code)
            script_path = Path(handle.name)

        return subprocess.run(
            [rscript, str(script_path), *arguments],
            check=False,
            capture_output=True,
            text=True,
        )
    finally:
        if script_path is not None:
            script_path.unlink(missing_ok=True)


def read_gcta_eigenvec(
    filename: Path,
) -> tuple[list[tuple[str, str]], list[str], np.ndarray]:
    """
    读取无表头的 GCTA eigenvec 文件。

    GCTA 格式：
        FID IID PC1 PC2 PC3 ...
    """
    if not filename.is_file():
        raise FileNotFoundError(f"GCTA eigenvec 文件不存在：{filename}")

    rows: list[list[str]] = []
    with filename.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            if rows and len(fields) != len(rows[0]):
                raise ValueError(
                    f"GCTA eigenvec 第 {line_number} 行有 "
                    f"{len(fields)} 列，预期为 {len(rows[0])} 列。"
                )
            rows.append(fields)

    if not rows:
        raise ValueError(f"GCTA eigenvec 文件为空：{filename}")

    if len(rows[0]) < 3:
        raise ValueError(
            "GCTA eigenvec 文件至少需要包含 FID、IID 和一个主成分。"
        )

    number_of_pcs = len(rows[0]) - 2
    pc_columns = [f"PC{i}" for i in range(1, number_of_pcs + 1)]
    sample_ids = [(row[0], row[1]) for row in rows]

    seen_ids: set[tuple[str, str]] = set()
    duplicated_ids: list[str] = []
    for sample_id in sample_ids:
        if sample_id in seen_ids and len(duplicated_ids) < 10:
            duplicated_ids.append(":".join(sample_id))
        seen_ids.add(sample_id)

    if duplicated_ids:
        raise ValueError(
            "检测到重复的 FID/IID 样本标识："
            + ", ".join(duplicated_ids)
        )

    try:
        pca_matrix = np.asarray(
            [row[2:] for row in rows],
            dtype=np.float64,
        )
    except ValueError as error:
        raise ValueError(
            "GCTA eigenvec 的主成分列包含非数值内容。"
        ) from error

    return sample_ids, pc_columns, pca_matrix


def read_gcta_eigenval(filename: Path) -> np.ndarray:
    """读取并校验 GCTA 单列 eigenval 文件。"""
    if not filename.is_file():
        raise FileNotFoundError(f"GCTA eigenval 文件不存在：{filename}")

    values: list[float] = []
    with filename.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 1:
                raise ValueError(
                    f"GCTA eigenval 第 {line_number} 行包含 "
                    f"{len(fields)} 列，预期为 1 列。"
                )
            try:
                values.append(float(fields[0]))
            except ValueError as error:
                raise ValueError(
                    f"GCTA eigenval 第 {line_number} 行包含非数值内容。"
                ) from error

    if not values:
        raise ValueError(f"GCTA eigenval 文件为空：{filename}")

    eigenvalues = np.asarray(values, dtype=np.float64)

    if not np.isfinite(eigenvalues).all():
        raise ValueError(
            "GCTA eigenval 文件包含缺失值、NaN 或无穷值。"
        )

    negative_count = int((eigenvalues < 0).sum())
    if negative_count:
        print(
            f"[WARNING] eigenval 中有 {negative_count} 个负特征值；"
            "计算解释方差时将其按 0 处理。",
            file=sys.stderr,
        )
        eigenvalues = np.maximum(eigenvalues, 0)

    if eigenvalues.sum() <= 0:
        raise ValueError(
            "GCTA eigenval 文件中的特征值总和必须大于 0。"
        )

    return eigenvalues


def draw_scree_plot(
    eigenval_path: Path,
    output_path: Path,
    number_of_pcs: int,
) -> None:
    """调用 WSL2 中的 R 和 ggplot2 绘制 PCA 碎石图。"""
    if number_of_pcs < 1:
        raise ValueError("PCA 碎石图至少需要一个特征值。")

    suffix = output_path.suffix.lower()
    if suffix not in {".pdf", ".png"}:
        raise ValueError(
            "--scree-plot 仅支持扩展名 .pdf 或 .png。"
        )

    rscript = shutil.which("Rscript")
    if rscript is None:
        raise RuntimeError(
            "未找到 Rscript；绘制碎石图需要 WSL2 中的 R。"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # R 代码用两个共享 PC 横轴的面板展示单个及累计解释方差。
    r_code = r"""
args <- commandArgs(trailingOnly = TRUE)
eigenval_file <- args[[1]]
output_file <- args[[2]]
n_requested <- as.integer(args[[3]])

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("R package 'ggplot2' is required.")
}

eigenvalues <- scan(eigenval_file, what = double(), quiet = TRUE)
eigenvalues <- pmax(eigenvalues, 0)
n_total <- length(eigenvalues)
n_show <- min(n_requested, n_total)
explained <- eigenvalues / sum(eigenvalues) * 100
cumulative <- cumsum(explained)

plot_data <- rbind(
    data.frame(
        PC = seq_len(n_show),
        Value = explained[seq_len(n_show)],
        Panel = "Individual explained variance (%)"
    ),
    data.frame(
        PC = seq_len(n_show),
        Value = cumulative[seq_len(n_show)],
        Panel = "Cumulative explained variance (%)"
    )
)

plot_data$Panel <- factor(
    plot_data$Panel,
    levels = c(
        "Individual explained variance (%)",
        "Cumulative explained variance (%)"
    )
)

individual_data <- plot_data[
    plot_data$Panel == "Individual explained variance (%)",
]
cumulative_data <- plot_data[
    plot_data$Panel == "Cumulative explained variance (%)",
]

p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = PC, y = Value)
    ) +
    ggplot2::geom_col(
        data = individual_data,
        width = 0.82,
        fill = "#3B6FB6"
    ) +
    ggplot2::geom_line(
        data = cumulative_data,
        linewidth = 0.8,
        color = "#D98B2B"
    )

if (n_show <= 100) {
    p <- p + ggplot2::geom_point(
        data = cumulative_data,
        size = 1.5,
        color = "#D98B2B"
    )
}

p <- p +
    ggplot2::facet_wrap(
        ggplot2::vars(Panel),
        ncol = 1,
        scales = "free_y"
    ) +
    ggplot2::scale_x_continuous(
        breaks = scales::breaks_pretty(n = 10),
        expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::scale_y_continuous(
        breaks = scales::breaks_pretty(n = 6),
        expand = ggplot2::expansion(mult = c(0, 0.06))
    ) +
    ggplot2::labs(
        title = "PCA Scree Plot",
        subtitle = sprintf(
            "First %d of %d PCs; percentages use all eigenvalues",
            n_show,
            n_total
        ),
        x = "Principal component",
        y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(
            face = "bold",
            color = "#252525"
        ),
        plot.subtitle = ggplot2::element_text(color = "#555555"),
        plot.caption = ggplot2::element_text(
            color = "#666666",
            hjust = 0
        ),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(
            color = "#E5E5E5",
            linewidth = 0.35
        ),
        panel.spacing.y = grid::unit(0.7, "lines"),
        strip.background = ggplot2::element_rect(
            fill = "#F2F2F2",
            color = "#D8D8D8"
        ),
        strip.text = ggplot2::element_text(
            face = "bold",
            color = "#333333"
        ),
        axis.title = ggplot2::element_text(color = "#333333"),
        axis.text = ggplot2::element_text(color = "#4D4D4D")
    )

ggplot2::ggsave(
    filename = output_file,
    plot = p,
    width = 10,
    height = 6.5,
    units = "in",
    dpi = 600,
    bg = "white"
)
"""

    completed = run_r_code(
        rscript=rscript,
        r_code=r_code,
        arguments=[
            str(eigenval_path),
            str(output_path),
            str(number_of_pcs),
        ],
    )

    if completed.returncode != 0:
        error_message = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(
            "R 绘制 PCA 碎石图失败："
            + (error_message or "未知错误")
        )

    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(
            f"R 未生成有效的 PCA 碎石图文件：{output_path}"
        )


def draw_umap_plot(
    umap_path: Path,
    sample_group_path: Path,
    color_mapping_path: Path | None,
    pdf_path: Path,
    png_path: Path,
) -> None:
    """调用 WSL2 中的 R，根据 Excel 样本分组绘制 UMAP 图。"""
    if not sample_group_path.is_file():
        raise FileNotFoundError(
            f"样本分组文件不存在：{sample_group_path}"
        )

    if sample_group_path.suffix.lower() not in {".xlsx", ".xlsm"}:
        raise ValueError(
            "--input-sample-group 必须是 .xlsx 或 .xlsm 格式的 Excel 文件。"
        )

    if color_mapping_path is not None and not color_mapping_path.is_file():
        raise FileNotFoundError(
            f"群体颜色映射文件不存在：{color_mapping_path}"
        )

    rscript = shutil.which("Rscript")
    if rscript is None:
        raise RuntimeError(
            "未找到 Rscript；绘制 UMAP 图需要 WSL2 中的 R。"
        )

    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.parent.mkdir(parents=True, exist_ok=True)

    # 绘图参数参考 PCA.R，并对未知分组自动补充分辨度较高的颜色。
    r_code = r"""
args <- commandArgs(trailingOnly = TRUE)
umap_file <- args[[1]]
sample_group_file <- args[[2]]
color_mapping_file <- args[[3]]
pdf_file <- args[[4]]
png_file <- args[[5]]

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("R package 'ggplot2' is required.")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("R package 'openxlsx' is required.")
}

umap_data <- utils::read.delim(
    umap_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
)

required_umap_columns <- c("IID", "UMAP1", "UMAP2")
missing_umap_columns <- setdiff(required_umap_columns, names(umap_data))
if (length(missing_umap_columns) > 0) {
    stop(sprintf(
        "UMAP output is missing required columns: %s",
        paste(missing_umap_columns, collapse = ", ")
    ))
}

umap_data$IID <- trimws(umap_data$IID)
umap_data$UMAP1 <- suppressWarnings(as.numeric(umap_data$UMAP1))
umap_data$UMAP2 <- suppressWarnings(as.numeric(umap_data$UMAP2))
if (any(!is.finite(umap_data$UMAP1)) ||
        any(!is.finite(umap_data$UMAP2))) {
    stop("UMAP1 or UMAP2 contains non-numeric or non-finite values.")
}

sample_groups <- openxlsx::read.xlsx(
    sample_group_file,
    sheet = 1,
    colNames = TRUE,
    check.names = FALSE
)

expected_headers <- c("SampleID", "Group")
if (ncol(sample_groups) < 2 ||
        !identical(names(sample_groups)[1:2], expected_headers)) {
    stop(
        "The first two Excel headers must be exactly SampleID and Group."
    )
}

sample_groups <- sample_groups[, expected_headers, drop = FALSE]
sample_groups$SampleID <- trimws(as.character(sample_groups$SampleID))
sample_groups$Group <- trimws(as.character(sample_groups$Group))

invalid_sample_ids <- is.na(sample_groups$SampleID) |
    sample_groups$SampleID == ""
if (any(invalid_sample_ids)) {
    stop("The SampleID column contains empty values.")
}

duplicated_sample_ids <- unique(
    sample_groups$SampleID[duplicated(sample_groups$SampleID)]
)
if (length(duplicated_sample_ids) > 0) {
    stop(sprintf(
        "The SampleID column contains duplicates: %s",
        paste(utils::head(duplicated_sample_ids, 10), collapse = ", ")
    ))
}

match_index <- match(umap_data$IID, sample_groups$SampleID)
matched_rows <- !is.na(match_index)
unmatched_umap_count <- sum(!matched_rows)
unused_group_count <- sum(
    !sample_groups$SampleID %in% umap_data$IID
)

plot_data <- umap_data[matched_rows, , drop = FALSE]
plot_data$Group <- sample_groups$Group[match_index[matched_rows]]
valid_group_rows <- !is.na(plot_data$Group) & plot_data$Group != ""
empty_group_count <- sum(!valid_group_rows)
plot_data <- plot_data[valid_group_rows, , drop = FALSE]

if (nrow(plot_data) == 0) {
    stop(
        "No UMAP samples matched non-empty groups by IID and SampleID."
    )
}

group_levels <- unique(
    sample_groups$Group[
        !is.na(sample_groups$Group) &
        sample_groups$Group != "" &
        sample_groups$Group %in% plot_data$Group
    ]
)
plot_data$Group <- factor(plot_data$Group, levels = group_levels)

default_colors <- c(
    "#C5000B",
    "#004586",
    "#579D1C",
    "#FF950E",
    "#7E0021",
    "#83CAFF",
    "#314004",
    "#AECF00",
    "#4B1F6F",
    "#DC3912",
    "#0084D1"
)

if (nzchar(color_mapping_file)) {
    color_mapping <- utils::read.table(
        color_mapping_file,
        header = FALSE,
        sep = "\t",
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE,
        colClasses = "character"
    )
    if (ncol(color_mapping) != 2) {
        stop(paste0(
            "The color mapping file must contain exactly two columns ",
            "without a header."
        ))
    }
    names(color_mapping) <- c("Group", "Color")
    color_mapping$Group <- trimws(color_mapping$Group)
    color_mapping$Color <- trimws(color_mapping$Color)

    if (any(color_mapping$Group == "") ||
            any(color_mapping$Color == "")) {
        stop("The color mapping file contains empty values.")
    }
    duplicated_mapping_groups <- unique(
        color_mapping$Group[duplicated(color_mapping$Group)]
    )
    if (length(duplicated_mapping_groups) > 0) {
        stop(sprintf(
            "The color mapping file contains duplicate groups: %s",
            paste(duplicated_mapping_groups, collapse = ", ")
        ))
    }
    invalid_colors <- color_mapping$Color[
        !grepl("^#[0-9A-Fa-f]{6}$", color_mapping$Color)
    ]
    if (length(invalid_colors) > 0) {
        stop(sprintf(
            "Invalid #RRGGBB color codes: %s",
            paste(unique(invalid_colors), collapse = ", ")
        ))
    }

    missing_color_groups <- setdiff(
        group_levels,
        color_mapping$Group
    )
    if (length(missing_color_groups) > 0) {
        message(sprintf(
            paste0(
                "[WARNING] 颜色映射文件缺少以下群体：%s。",
                "程序退出。"
            ),
            paste(missing_color_groups, collapse = ", ")
        ))
        quit(save = "no", status = 1)
    }
    group_colors <- stats::setNames(
        color_mapping$Color,
        color_mapping$Group
    )[group_levels]
} else {
    if (length(group_levels) > length(default_colors)) {
        message(sprintf(
            paste0(
                "[WARNING] 检测到 %d 个群体，但内置颜色仅有 %d 个。",
                "请使用 --color-mapping-group 提供完整映射；程序退出。"
            ),
            length(group_levels),
            length(default_colors)
        ))
        quit(save = "no", status = 1)
    }
    group_colors <- stats::setNames(
        default_colors[seq_along(group_levels)],
        group_levels
    )
}

p <- ggplot2::ggplot(plot_data,ggplot2::aes(x = UMAP1, y = UMAP2)) +
    ggplot2::geom_point(ggplot2::aes(color = Group),size = 0.6) +
    ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
    ggplot2::labs(x = "UMAP1", y = "UMAP2") +
    ggplot2::theme_bw() +
    ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(family = "sans",size = 12,colour = "black",face = "plain",angle = 0),
        axis.title.x = ggplot2::element_text(size = 12,margin = ggplot2::margin(t = 3)),
        axis.text.y = ggplot2::element_text(family = "sans",size = 12,colour = "black",face = "plain",angle = 0),
        axis.title.y = ggplot2::element_text(size = 12,margin = ggplot2::margin(r = 3)),
        legend.position = "right",
        legend.background = ggplot2::element_rect(fill = NA,color = NA),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 6),
        legend.key.height = grid::unit(10, "pt")
    )

ggplot2::ggsave(
    filename = pdf_file,plot = p,width = 10,height = 6,units = "cm",bg = "white")
ggplot2::ggsave(filename = png_file,plot = p,width = 10,height = 6,units = "cm",dpi = 600,bg = "white")

message(sprintf("[INFO] UMAP 图使用样本数：%d", nrow(plot_data)))
if (unmatched_umap_count > 0) {
    message(sprintf(
        "[WARNING] UMAP 中有 %d 个 IID 未在分组文件中找到。",
        unmatched_umap_count
    ))
}
if (empty_group_count > 0) {
    message(sprintf(
        "[WARNING] 有 %d 个已匹配样本的 Group 为空。",
        empty_group_count
    ))
}
if (unused_group_count > 0) {
    message(sprintf(
        "[WARNING] 分组文件中有 %d 个 SampleID 不在 UMAP 结果中。",
        unused_group_count
    ))
}
"""

    completed = run_r_code(
        rscript=rscript,
        r_code=r_code,
        arguments=[
            str(umap_path),
            str(sample_group_path),
            str(color_mapping_path) if color_mapping_path else "",
            str(pdf_path),
            str(png_path),
        ],
    )

    if completed.returncode != 0:
        error_message = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(
            "R 绘制 UMAP 图失败："
            + (error_message or "未知错误")
        )

    missing_plot_paths = [
        path
        for path in (pdf_path, png_path)
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing_plot_paths:
        raise RuntimeError(
            "R 未生成以下有效 UMAP 图文件："
            + ", ".join(str(path) for path in missing_plot_paths)
        )

    r_messages = completed.stderr.strip() or completed.stdout.strip()
    if r_messages:
        print(r_messages, file=sys.stderr)


def load_umap_module():
    """延迟导入 umap，并在缺少依赖时提供清晰提示。"""
    try:
        import umap
    except ModuleNotFoundError as error:
        raise RuntimeError(
            "未安装 Python 包 'umap-learn'。请在目标虚拟环境中安装后"
            "再运行 UMAP。"
        ) from error
    return umap


def main() -> None:
    args = parse_args()
    sample_group_path = Path(args.input_sample_group)
    color_mapping_path = (
        Path(args.color_mapping_group)
        if args.color_mapping_group
        else None
    )
    if not sample_group_path.is_file():
        raise FileNotFoundError(
            f"样本分组文件不存在：{sample_group_path}"
        )
    if sample_group_path.suffix.lower() not in {".xlsx", ".xlsm"}:
        raise ValueError(
            "--input-sample-group 必须是 .xlsx 或 .xlsm 格式的 "
            "Excel 文件。"
        )
    if color_mapping_path is not None and not color_mapping_path.is_file():
        raise FileNotFoundError(
            f"群体颜色映射文件不存在：{color_mapping_path}"
        )

    if args.n_pc < 1:
        raise ValueError("--n-pc 必须至少为 1。")

    if args.n_components < 1:
        raise ValueError("--n-components 必须至少为 1。")

    if args.n_components < 2:
        raise ValueError(
            "绘制 UMAP 图至少需要两个输出维度；"
            "请将 --n-components 设置为不小于 2。"
        )

    if args.n_neighbors < 2:
        raise ValueError("--n-neighbors 必须至少为 2。")

    if not np.isfinite(args.min_dist) or args.min_dist < 0:
        raise ValueError("--min-dist 必须是大于或等于 0 的有限数值。")

    eigenvec_path, eigenval_path = get_gcta_paths(args.pca_result)
    eigenvalues = read_gcta_eigenval(eigenval_path)

    scree_plot_path = Path(
        args.scree_plot
        if args.scree_plot
        else f"{args.pca_result}.scree_plot.pdf"
    )

    draw_scree_plot(
        eigenval_path=eigenval_path,
        output_path=scree_plot_path,
        number_of_pcs=len(eigenvalues),
    )
    print(
        f"[INFO] PCA 碎石图（全部 PC）：{scree_plot_path}",
        file=sys.stderr,
    )

    sample_ids, available_pc_columns, pca_matrix = (
        read_gcta_eigenvec(eigenvec_path)
    )
    number_of_samples = len(sample_ids)

    if len(eigenvalues) < len(available_pc_columns):
        raise ValueError(
            f"eigenval 仅包含 {len(eigenvalues)} 个特征值，但 "
            f"eigenvec 包含 {len(available_pc_columns)} 个 PC。"
        )

    if args.n_pc > len(available_pc_columns):
        raise ValueError(
            f"--n-pc 设置为 {args.n_pc}，但 GCTA eigenvec "
            f"仅包含 {len(available_pc_columns)} 个 PC。"
        )

    if args.n_neighbors >= number_of_samples:
        raise ValueError(
            "--n-neighbors 必须小于样本数。"
            f"当前 n_neighbors={args.n_neighbors}，"
            f"n_samples={number_of_samples}。"
        )

    selected_pc_columns = available_pc_columns[:args.n_pc]
    pca_matrix = pca_matrix[:, :args.n_pc]

    if not np.isfinite(pca_matrix).all():
        raise ValueError(
            "PCA 矩阵包含缺失值、NaN 或无穷值。"
        )

    umap = load_umap_module()
    reducer = umap.UMAP(
        n_neighbors=args.n_neighbors,
        n_components=args.n_components,
        metric=args.metric,
        min_dist=args.min_dist,
        init=args.init,
        n_epochs=args.n_epochs,
        random_state=args.seed,
        transform_seed=args.seed,
        n_jobs=1,
    )

    embedding = reducer.fit_transform(pca_matrix)

    output_path = Path(args.output_umap)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(
            [
                "FID",
                "IID",
                *[
                    f"UMAP{index + 1}"
                    for index in range(args.n_components)
                ],
            ]
        )
        for (fid, iid), coordinates in zip(sample_ids, embedding):
            writer.writerow(
                [
                    fid,
                    iid,
                    *[
                        f"{coordinate:.10f}"
                        for coordinate in coordinates
                    ],
                ]
            )

    umap_plot_pdf_path, umap_plot_png_path = get_umap_plot_paths(
        args.output_umap_plot
    )
    draw_umap_plot(
        umap_path=output_path,
        sample_group_path=sample_group_path,
        color_mapping_path=color_mapping_path,
        pdf_path=umap_plot_pdf_path,
        png_path=umap_plot_png_path,
    )

    print(f"[INFO] 样本数：{number_of_samples}", file=sys.stderr)
    print(
        f"[INFO] 输入 PC：{selected_pc_columns[0]}-"
        f"{selected_pc_columns[-1]}",
        file=sys.stderr,
    )
    print(f"[INFO] n_neighbors：{args.n_neighbors}", file=sys.stderr)
    print(f"[INFO] min_dist：{args.min_dist}", file=sys.stderr)
    print(f"[INFO] metric：{args.metric}", file=sys.stderr)
    print(f"[INFO] 随机种子：{args.seed}", file=sys.stderr)
    print(f"[INFO] UMAP 输出轴数：{args.n_components}", file=sys.stderr)
    print(f"[INFO] UMAP 输出文件：{output_path}", file=sys.stderr)
    print(f"[INFO] UMAP PDF 图：{umap_plot_pdf_path}", file=sys.stderr)
    print(f"[INFO] UMAP PNG 图：{umap_plot_png_path}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        sys.exit(1)
