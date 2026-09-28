#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

struct Arguments {
    std::unordered_map<std::string, std::string> values;
};

static std::string timestamp() {
    std::time_t now = std::time(nullptr);
    std::tm local{};
    localtime_r(&now, &local);
    char buffer[32];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", &local);
    return buffer;
}

static void log_info(const std::string& message) {
    const std::string line = "[INFO " + timestamp() + "] " + message + "\n";
    (void)::write(STDERR_FILENO, line.data(), line.size());
}

[[noreturn]] static void fail(const std::string& message) {
    throw std::runtime_error(message);
}

static Arguments parse_arguments(int argc, char** argv, int start) {
    Arguments parsed;
    for (int index = start; index < argc; index += 2) {
        if (index + 1 >= argc || std::string(argv[index]).rfind("--", 0) != 0) {
            fail("C++辅助程序参数格式无效: " + std::string(argv[index]));
        }
        parsed.values[argv[index]] = argv[index + 1];
    }
    return parsed;
}

static std::string required(const Arguments& args, const std::string& key) {
    const auto iterator = args.values.find(key);
    if (iterator == args.values.end() || iterator->second.empty()) {
        fail("C++辅助程序缺少参数: " + key);
    }
    return iterator->second;
}

static uint64_t positive_u64(const std::string& value, const std::string& key) {
    std::size_t parsed = 0;
    uint64_t number = 0;
    try {
        number = std::stoull(value, &parsed);
    } catch (...) {
        fail(key + "必须为正整数: " + value);
    }
    if (parsed != value.size() || number == 0) fail(key + "必须为正整数: " + value);
    return number;
}

static uint64_t nonnegative_u64(const std::string& value, const std::string& key) {
    if (value.empty() || !std::all_of(value.begin(), value.end(), [](const unsigned char character) {
            return character >= '0' && character <= '9';
        })) {
        fail(key + "必须为非负整数: " + value);
    }
    std::size_t parsed = 0;
    uint64_t number = 0;
    try {
        number = std::stoull(value, &parsed);
    } catch (...) {
        fail(key + "必须为非负整数: " + value);
    }
    if (parsed != value.size()) fail(key + "必须为非负整数: " + value);
    return number;
}

static std::vector<std::string> read_nonempty_lines(const std::string& path) {
    std::ifstream input(path);
    if (!input) fail("无法读取文件: " + path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) lines.push_back(line);
    }
    return lines;
}

static std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        if (end == std::string::npos) {
            fields.emplace_back(line.substr(begin));
            break;
        }
        fields.emplace_back(line.substr(begin, end - begin));
        begin = end + 1;
    }
    return fields;
}

struct Gene {
    std::string chrom;
    uint64_t start = 0;
    uint64_t end = 0;
    std::string id;
    uint64_t length = 0;
    std::vector<std::pair<uint64_t, uint64_t>> cds_intervals;
};

static std::vector<Gene> read_gene_bed(const std::string& path) {
    std::ifstream input(path);
    if (!input) fail("无法读取标准化基因区间文件: " + path);
    std::vector<Gene> genes;
    std::string line;
    std::string previous_chrom;
    uint64_t previous_start = 0;
    bool have_previous = false;
    std::unordered_map<std::string, bool> seen_gene_ids;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() != 6) {
            fail("标准化基因区间文件必须恰好六列Chr/Start/End/GeneID/RegionLength/RegionIntervals: " + line);
        }
        Gene gene;
        gene.chrom = fields[0];
        gene.id = fields[3];
        try {
            std::size_t parsed_start = 0;
            std::size_t parsed_end = 0;
            std::size_t parsed_length = 0;
            gene.start = std::stoull(fields[1], &parsed_start);
            gene.end = std::stoull(fields[2], &parsed_end);
            gene.length = std::stoull(fields[4], &parsed_length);
            if (parsed_start != fields[1].size() || parsed_end != fields[2].size() ||
                parsed_length != fields[4].size()) fail("坐标或RegionLength不是整数");
        } catch (...) {
            fail("标准化基因区间坐标无效: " + line);
        }
        if (gene.chrom.empty() || gene.id.empty() || gene.start == 0 || gene.end < gene.start ||
            gene.length == 0 || fields[5].empty()) {
            fail("标准化基因区间、RegionLength或GeneID无效: " + line);
        }
        uint64_t calculated_length = 0;
        std::size_t interval_begin = 0;
        while (interval_begin < fields[5].size()) {
            const std::size_t interval_end = fields[5].find(',', interval_begin);
            const std::string interval = fields[5].substr(
                interval_begin, interval_end == std::string::npos ? std::string::npos : interval_end - interval_begin);
            const std::size_t dash = interval.find('-');
            if (dash == std::string::npos || dash == 0 || dash + 1 >= interval.size() ||
                interval.find('-', dash + 1) != std::string::npos) {
                fail("RegionIntervals格式无效: " + interval);
            }
            uint64_t cds_start = 0;
            uint64_t cds_end = 0;
            try {
                std::size_t parsed_cds_start = 0;
                std::size_t parsed_cds_end = 0;
                cds_start = std::stoull(interval.substr(0, dash), &parsed_cds_start);
                cds_end = std::stoull(interval.substr(dash + 1), &parsed_cds_end);
                if (parsed_cds_start != dash || parsed_cds_end != interval.size() - dash - 1) {
                    fail("基因区间不是整数");
                }
            } catch (...) {
                fail("RegionIntervals坐标无效: " + interval);
            }
            if (cds_start == 0 || cds_end < cds_start || cds_start < gene.start || cds_end > gene.end) {
                fail("RegionIntervals超出基因区间外边界: " + interval);
            }
            if (!gene.cds_intervals.empty() && cds_start <= gene.cds_intervals.back().second) {
                fail("RegionIntervals必须按坐标排序且互不重叠: " + gene.id);
            }
            gene.cds_intervals.emplace_back(cds_start, cds_end);
            const uint64_t interval_length = cds_end - cds_start + 1;
            if (calculated_length > std::numeric_limits<uint64_t>::max() - interval_length) {
                fail("RegionLength求和溢出: " + gene.id);
            }
            calculated_length += interval_length;
            if (interval_end == std::string::npos) break;
            interval_begin = interval_end + 1;
        }
        if (gene.cds_intervals.empty() || gene.cds_intervals.front().first != gene.start ||
            gene.cds_intervals.back().second != gene.end || calculated_length != gene.length) {
            fail("RegionLength或外边界与RegionIntervals不一致: " + gene.id);
        }
        if (seen_gene_ids[gene.id]) fail("标准化基因区间文件包含重复GeneID: " + gene.id);
        seen_gene_ids[gene.id] = true;
        if (have_previous && gene.chrom == previous_chrom && gene.start < previous_start) {
            fail("标准化基因区间文件未按染色体和Start排序: " + line);
        }
        previous_chrom = gene.chrom;
        previous_start = gene.start;
        have_previous = true;
        genes.push_back(std::move(gene));
    }
    if (genes.empty()) fail("标准化基因区间文件中没有基因: " + path);
    return genes;
}

static bool gene_contains_position(const Gene& gene, uint64_t position) {
    const auto iterator = std::upper_bound(
        gene.cds_intervals.begin(), gene.cds_intervals.end(), position,
        [](uint64_t value, const std::pair<uint64_t, uint64_t>& interval) {
            return value < interval.first;
        });
    if (iterator == gene.cds_intervals.begin()) return false;
    return std::prev(iterator)->second >= position;
}

static bool same_gene_definition(const Gene& first, const Gene& second) {
    return first.chrom == second.chrom && first.start == second.start && first.end == second.end &&
           first.length == second.length && first.cds_intervals == second.cds_intervals;
}

struct PatternNode {
    std::unordered_map<unsigned char, std::size_t> next;
    std::size_t failure = 0;
    std::vector<std::size_t> outputs;
};

static std::vector<std::string> read_gene_ids(const std::string& path) {
    std::ifstream input(path);
    if (!input) fail("无法读取GeneID列表: " + path);
    std::vector<std::string> ids;
    std::unordered_set<std::string> seen;
    std::string line;
    if (!std::getline(input, line)) fail("GeneID列表为空，第一行必须是表头: " + path);
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) fail("GeneID列表第一行表头不能为空: " + path);
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line.find_first_of("\t ") != std::string::npos) {
            fail("GeneID列表必须每行恰好一个不含空白字符的GeneID: " + line);
        }
        if (!seen.insert(line).second) fail("GeneID列表包含重复ID: " + line);
        ids.push_back(line);
    }
    if (ids.empty()) fail("GeneID列表除表头外没有有效GeneID: " + path);
    return ids;
}

static std::vector<PatternNode> build_pattern_automaton(const std::vector<std::string>& patterns) {
    std::vector<PatternNode> nodes(1);
    for (std::size_t pattern_index = 0; pattern_index < patterns.size(); ++pattern_index) {
        std::size_t state = 0;
        for (const unsigned char character : patterns[pattern_index]) {
            auto iterator = nodes[state].next.find(character);
            if (iterator == nodes[state].next.end()) {
                const std::size_t new_state = nodes.size();
                nodes[state].next.emplace(character, new_state);
                nodes.emplace_back();
                state = new_state;
            } else {
                state = iterator->second;
            }
        }
        nodes[state].outputs.push_back(pattern_index);
    }
    std::queue<std::size_t> pending;
    for (const auto& edge : nodes[0].next) {
        nodes[edge.second].failure = 0;
        pending.push(edge.second);
    }
    while (!pending.empty()) {
        const std::size_t state = pending.front();
        pending.pop();
        for (const auto& edge : nodes[state].next) {
            const unsigned char character = edge.first;
            const std::size_t next_state = edge.second;
            std::size_t failure = nodes[state].failure;
            while (failure != 0 && nodes[failure].next.find(character) == nodes[failure].next.end()) {
                failure = nodes[failure].failure;
            }
            const auto failure_edge = nodes[failure].next.find(character);
            if (failure_edge != nodes[failure].next.end() && failure_edge->second != next_state) {
                nodes[next_state].failure = failure_edge->second;
            } else {
                nodes[next_state].failure = 0;
            }
            const auto& inherited = nodes[nodes[next_state].failure].outputs;
            nodes[next_state].outputs.insert(nodes[next_state].outputs.end(), inherited.begin(), inherited.end());
            pending.push(next_state);
        }
    }
    return nodes;
}

static int extract_regions_mode(const Arguments& args) {
    const std::string gene_id_path = required(args, "--gene-id");
    const std::string gff3_path = required(args, "--gff3");
    const std::string variantion_source = required(args, "--variantion-source");
    const uint64_t upstream_bp = nonnegative_u64(required(args, "--upstream"), "--upstream");
    const uint64_t downstream_bp = nonnegative_u64(required(args, "--downstream"), "--downstream");
    const fs::path output_path = required(args, "--output");
    if (variantion_source != "cds" && variantion_source != "mRNA" && variantion_source != "gene") {
        fail("--variantion-source必须是cds、mRNA或gene: " + variantion_source);
    }
    const std::string target_feature = variantion_source == "cds" ? "CDS" : "mRNA";
    const auto gene_ids = read_gene_ids(gene_id_path);
    const auto automaton = build_pattern_automaton(gene_ids);
    std::ifstream gff3(gff3_path);
    if (!gff3) fail("无法读取GFF3: " + gff3_path);

    struct ExtractedGene {
        std::string chrom;
        std::string strand;
        std::size_t chrom_rank = 0;
        std::vector<std::pair<uint64_t, uint64_t>> intervals;
    };
    std::vector<ExtractedGene> extracted(gene_ids.size());
    std::unordered_map<std::string, std::size_t> chromosome_ranks;
    std::size_t next_chromosome_rank = 0;
    std::vector<uint64_t> last_seen(gene_ids.size(), 0);
    uint64_t feature_row_number = 0;
    uint64_t matched_feature_rows = 0;
    uint64_t gene_assignments = 0;
    std::string line;
    while (std::getline(gff3, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tab(line);
        if (fields.size() != 9) fail("GFF3数据行必须恰好九列: " + line);
        if (fields[2] != target_feature) continue;
        ++feature_row_number;
        uint64_t start = 0;
        uint64_t end = 0;
        try {
            std::size_t parsed_start = 0;
            std::size_t parsed_end = 0;
            start = std::stoull(fields[3], &parsed_start);
            end = std::stoull(fields[4], &parsed_end);
            if (parsed_start != fields[3].size() || parsed_end != fields[4].size()) fail("GFF3坐标不是整数");
        } catch (...) {
            fail("GFF3 " + target_feature + "坐标无效: " + line);
        }
        if (fields[0].empty() || start == 0 || end < start) fail("GFF3 " + target_feature + "区间无效: " + line);
        if (variantion_source == "gene") {
            if (fields[6] != "+" && fields[6] != "-") {
                fail("gene模式要求mRNA记录第7列链方向必须是+或-: " + line);
            }
            const uint64_t left_extension = fields[6] == "+" ? upstream_bp : downstream_bp;
            const uint64_t right_extension = fields[6] == "+" ? downstream_bp : upstream_bp;
            start = start > left_extension ? start - left_extension : 1;
            if (end > std::numeric_limits<uint64_t>::max() - right_extension) {
                fail("gene模式扩展后终止坐标溢出: " + line);
            }
            end += right_extension;
        }
        auto rank_iterator = chromosome_ranks.find(fields[0]);
        if (rank_iterator == chromosome_ranks.end()) {
            rank_iterator = chromosome_ranks.emplace(fields[0], next_chromosome_rank++).first;
        }
        std::size_t state = 0;
        std::vector<std::size_t> row_matches;
        for (const unsigned char character : fields[8]) {
            while (state != 0 && automaton[state].next.find(character) == automaton[state].next.end()) {
                state = automaton[state].failure;
            }
            const auto edge = automaton[state].next.find(character);
            if (edge != automaton[state].next.end()) state = edge->second;
            for (const std::size_t pattern_index : automaton[state].outputs) {
                if (last_seen[pattern_index] != feature_row_number) {
                    last_seen[pattern_index] = feature_row_number;
                    row_matches.push_back(pattern_index);
                }
            }
        }
        if (row_matches.empty()) continue;
        ++matched_feature_rows;
        for (const std::size_t pattern_index : row_matches) {
            auto& gene = extracted[pattern_index];
            if (gene.chrom.empty()) {
                gene.chrom = fields[0];
                gene.chrom_rank = rank_iterator->second;
            } else if (gene.chrom != fields[0]) {
                fail("同一GeneID匹配到多条染色体: " + gene_ids[pattern_index]);
            }
            if (variantion_source == "gene") {
                if (gene.strand.empty()) gene.strand = fields[6];
                else if (gene.strand != fields[6]) fail("同一GeneID的mRNA记录链方向不一致: " + gene_ids[pattern_index]);
            }
            gene.intervals.emplace_back(start, end);
            ++gene_assignments;
        }
    }

    std::vector<std::size_t> missing;
    std::vector<Gene> genes;
    genes.reserve(gene_ids.size());
    std::vector<std::size_t> gene_chromosome_ranks;
    gene_chromosome_ranks.reserve(gene_ids.size());
    for (std::size_t index = 0; index < gene_ids.size(); ++index) {
        auto& source = extracted[index];
        if (source.intervals.empty()) {
            missing.push_back(index);
            continue;
        }
        std::sort(source.intervals.begin(), source.intervals.end());
        std::vector<std::pair<uint64_t, uint64_t>> merged;
        for (const auto& interval : source.intervals) {
            if (merged.empty() || interval.first > merged.back().second + 1) {
                merged.push_back(interval);
            } else if (interval.second > merged.back().second) {
                merged.back().second = interval.second;
            }
        }
        Gene gene;
        gene.chrom = source.chrom;
        gene.start = merged.front().first;
        gene.end = merged.back().second;
        gene.id = gene_ids[index];
        gene.cds_intervals = std::move(merged);
        for (const auto& interval : gene.cds_intervals) gene.length += interval.second - interval.first + 1;
        genes.push_back(std::move(gene));
        gene_chromosome_ranks.push_back(source.chrom_rank);
    }
    if (!missing.empty()) {
        std::string message = "以下GeneID未在GFF3第9列匹配到" + target_feature + "记录";
        const std::size_t display = std::min<std::size_t>(missing.size(), 20);
        for (std::size_t index = 0; index < display; ++index) message += (index == 0 ? ": " : ",") + gene_ids[missing[index]];
        if (missing.size() > display) message += "（另有" + std::to_string(missing.size() - display) + "个未显示）";
        fail(message);
    }
    std::vector<std::size_t> order(genes.size());
    for (std::size_t index = 0; index < order.size(); ++index) order[index] = index;
    std::sort(order.begin(), order.end(), [&](std::size_t first, std::size_t second) {
        if (gene_chromosome_ranks[first] != gene_chromosome_ranks[second]) {
            return gene_chromosome_ranks[first] < gene_chromosome_ranks[second];
        }
        if (genes[first].start != genes[second].start) return genes[first].start < genes[second].start;
        if (genes[first].end != genes[second].end) return genes[first].end < genes[second].end;
        return genes[first].id < genes[second].id;
    });
    std::ofstream output(output_path, std::ios::trunc);
    if (!output) fail("无法创建标准化基因区间文件: " + output_path.string());
    for (const std::size_t index : order) {
        const auto& gene = genes[index];
        output << gene.chrom << '\t' << gene.start << '\t' << gene.end << '\t' << gene.id << '\t'
               << gene.length << '\t';
        for (std::size_t interval = 0; interval < gene.cds_intervals.size(); ++interval) {
            if (interval != 0) output << ',';
            output << gene.cds_intervals[interval].first << '-' << gene.cds_intervals[interval].second;
        }
        output << '\n';
    }
    output.close();
    if (!output) fail("写入标准化基因区间文件失败");
    log_info("完成GFF3 " + variantion_source + "区间提取: requested_genes=" + std::to_string(gene_ids.size()) +
             ", matched_feature_rows=" + std::to_string(matched_feature_rows) +
             ", gene_assignments=" + std::to_string(gene_assignments));
    return 0;
}

enum class GenotypeState : uint8_t { HomRef, HomAlt, Heterozygous, Missing };

static GenotypeState classify_gt(const std::string& genotype) {
    const std::size_t separator = genotype.find_first_of("/|");
    if (separator == std::string::npos || separator == 0 || separator + 1 >= genotype.size()) {
        return GenotypeState::Missing;
    }
    if (genotype.find_first_of("/|", separator + 1) != std::string::npos) {
        return GenotypeState::Missing;
    }
    const std::string allele1 = genotype.substr(0, separator);
    const std::string allele2 = genotype.substr(separator + 1);
    if ((allele1 != "0" && allele1 != "1") || (allele2 != "0" && allele2 != "1")) {
        return GenotypeState::Missing;
    }
    if (allele1 != allele2) return GenotypeState::Heterozygous;
    return allele1 == "0" ? GenotypeState::HomRef : GenotypeState::HomAlt;
}

struct SampleMasks {
    std::vector<uint64_t> hom_ref;
    std::vector<uint64_t> hom_alt;
    std::vector<uint64_t> missing;
};

struct GroupBuffer {
    std::string first_sample;
    std::size_t columns = 0;
    std::vector<uint16_t> ndiff;
    std::vector<uint16_t> nmiss;
    std::ofstream ndiff_output;
    std::ofstream nmiss_output;
};

static uint64_t pair_count_for_samples(uint64_t sample_count) {
    return sample_count * (sample_count - 1) / 2;
}

static void allocate_group_buffers(std::vector<GroupBuffer>& groups, std::size_t rows) {
    for (auto& group : groups) {
        if (group.columns && rows > std::numeric_limits<std::size_t>::max() / group.columns) {
            fail("基因计数缓冲区大小溢出");
        }
        group.ndiff.assign(rows * group.columns, 0);
        group.nmiss.assign(rows * group.columns, 0);
    }
}

static int count_mode(const Arguments& args) {
    const std::string chromosome = required(args, "--chrom");
    const std::string gene_bed_path = required(args, "--gene-bed");
    const uint64_t memory_budget = positive_u64(required(args, "--memory-bytes"), "--memory-bytes");
    const uint64_t threads = positive_u64(required(args, "--threads"), "--threads");
    const std::string samples_path = required(args, "--samples");
    const fs::path output_dir = required(args, "--output-dir");
    const std::vector<std::string> samples = read_nonempty_lines(samples_path);
    const std::vector<Gene> genes = read_gene_bed(gene_bed_path);
    if (samples.size() < 2) fail("基因计数至少需要两个样本");
    for (const auto& gene : genes) {
        if (gene.chrom != chromosome) fail("染色体基因区间文件包含非目标染色体: " + gene.chrom);
    }

#ifdef _OPENMP
    omp_set_num_threads(static_cast<int>(std::min<uint64_t>(threads, std::numeric_limits<int>::max())));
#else
    if (threads > 1) log_info("C++辅助程序未启用OpenMP，将使用单线程计数。");
#endif

    fs::create_directories(output_dir);
    const uint64_t total_genes = genes.size();
    const uint64_t pair_count = pair_count_for_samples(samples.size());
    if (pair_count > std::numeric_limits<uint64_t>::max() / total_genes / 4) {
        fail("染色体基因计数矩阵大小溢出");
    }
    const uint64_t matrix_bytes = total_genes * pair_count * sizeof(uint16_t) * 2;
    if (matrix_bytes > memory_budget * 9 / 10) {
        fail("为支持重叠基因，单条染色体需要同时保存全部基因计数；预计需要" +
             std::to_string(matrix_bytes) + " bytes，超过本次内存预算的90%=" +
             std::to_string(memory_budget * 9 / 10) + " bytes");
    }
    std::vector<GroupBuffer> groups(samples.size() - 1);
    for (std::size_t first = 0; first + 1 < samples.size(); ++first) {
        groups[first].first_sample = samples[first];
        groups[first].columns = samples.size() - first - 1;
    }
    try {
        allocate_group_buffers(groups, static_cast<std::size_t>(total_genes));
    } catch (const std::bad_alloc&) {
        fail("无法分配单条染色体的基因计数矩阵内存，预计需要" +
             std::to_string(matrix_bytes) + " bytes");
    }
    log_info("基因计数 " + chromosome + ": samples=" + std::to_string(samples.size()) +
             ", pairs=" + std::to_string(pair_count) + ", genes=" + std::to_string(total_genes) +
             ", matrix_memory=" + std::to_string(matrix_bytes) + " bytes" +
             ", memory_budget=" + std::to_string(memory_budget) + " bytes");

    struct ActiveGene {
        std::size_t index = 0;
        uint64_t variant_count = 0;
        std::vector<SampleMasks> masks;
    };
    std::vector<ActiveGene> active_genes;
    std::size_t gene_cursor = 0;
    uint64_t record_count = 0;
    uint64_t in_gene_record_count = 0;
    uint64_t gene_assignment_count = 0;
    uint64_t genes_with_variants = 0;
    uint64_t duplicate_position_count = 0;
    uint64_t last_position = 0;
    bool have_last_position = false;
    std::atomic<bool> overflow(false);
    std::vector<GenotypeState> states(samples.size());

    auto finalize_gene = [&](ActiveGene& active) {
        const std::size_t row = active.index;
        const std::size_t words = active.masks.empty() ? 0 : active.masks.front().hom_ref.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int64_t first_signed = 0; first_signed < static_cast<int64_t>(samples.size() - 1);
             ++first_signed) {
            const std::size_t first = static_cast<std::size_t>(first_signed);
            auto& group = groups[first];
            for (std::size_t second = first + 1; second < samples.size(); ++second) {
                uint64_t ndiff = 0;
                uint64_t nmiss = 0;
                for (std::size_t word = 0; word < words; ++word) {
                    const uint64_t missing_bits = active.masks[first].missing[word] |
                                                  active.masks[second].missing[word];
                    uint64_t diff_bits =
                        (active.masks[first].hom_ref[word] & active.masks[second].hom_alt[word]) |
                        (active.masks[first].hom_alt[word] & active.masks[second].hom_ref[word]);
                    nmiss += static_cast<uint64_t>(__builtin_popcountll(missing_bits));
                    ndiff += static_cast<uint64_t>(__builtin_popcountll(diff_bits));
                }
                if (ndiff > std::numeric_limits<uint16_t>::max() ||
                    nmiss > std::numeric_limits<uint16_t>::max()) {
                    overflow.store(true, std::memory_order_relaxed);
                    continue;
                }
                const std::size_t column = second - first - 1;
                group.ndiff[row * group.columns + column] = static_cast<uint16_t>(ndiff);
                group.nmiss[row * group.columns + column] = static_cast<uint16_t>(nmiss);
            }
        }
        if (overflow.load(std::memory_order_relaxed)) {
            fail("单个基因的N_diff或N_miss超过uint16上限65535: " +
                 genes[active.index].id);
        }
        if (active.variant_count > 0) ++genes_with_variants;
    };

    std::string line;
    while (std::getline(std::cin, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_tab(line);
        if (fields.size() != samples.size() + 2) {
            fail("bcftools GT数据列数异常：预期" + std::to_string(samples.size() + 2) +
                 "列，实际" + std::to_string(fields.size()) + "列");
        }
        if (fields[0] != chromosome) {
            fail("染色体VCF包含非目标染色体记录：预期" + chromosome + "，实际" + fields[0]);
        }
        uint64_t position = 0;
        try {
            std::size_t parsed = 0;
            position = std::stoull(fields[1], &parsed);
            if (parsed != fields[1].size()) fail("POS不是整数");
        } catch (...) {
            fail("VCF POS无效: " + fields[1]);
        }
        if (position == 0) fail("VCF POS必须大于0: " + fields[1]);
        if (have_last_position && position < last_position) {
            fail("VCF位置未排序: " + chromosome + ":" + std::to_string(position) +
                 " < " + std::to_string(last_position));
        }
        if (have_last_position && position == last_position) ++duplicate_position_count;
        last_position = position;
        have_last_position = true;
        ++record_count;
        auto active_iterator = active_genes.begin();
        while (active_iterator != active_genes.end()) {
            if (genes[active_iterator->index].end < position) {
                finalize_gene(*active_iterator);
                active_iterator = active_genes.erase(active_iterator);
            } else {
                ++active_iterator;
            }
        }
        while (gene_cursor < genes.size() && genes[gene_cursor].start <= position) {
            if (genes[gene_cursor].end >= position) {
                ActiveGene active;
                active.index = gene_cursor;
                active.masks.resize(samples.size());
                active_genes.push_back(std::move(active));
            }
            ++gene_cursor;
        }
        if (active_genes.empty()) continue;
        bool position_in_cds = false;
        for (const auto& active : active_genes) {
            if (gene_contains_position(genes[active.index], position)) {
                position_in_cds = true;
                break;
            }
        }
        if (!position_in_cds) continue;
        for (std::size_t sample = 0; sample < samples.size(); ++sample) {
            states[sample] = classify_gt(fields[sample + 2]);
        }
        uint64_t assignments_this_position = 0;
        for (auto& active : active_genes) {
            if (!gene_contains_position(genes[active.index], position)) continue;
            const std::size_t word = static_cast<std::size_t>(active.variant_count / 64);
            const uint64_t bit = uint64_t{1} << (active.variant_count % 64);
            if (active.variant_count % 64 == 0) {
                for (auto& sample_mask : active.masks) {
                    sample_mask.hom_ref.push_back(0);
                    sample_mask.hom_alt.push_back(0);
                    sample_mask.missing.push_back(0);
                }
            }
            for (std::size_t sample = 0; sample < samples.size(); ++sample) {
                switch (states[sample]) {
                    case GenotypeState::HomRef: active.masks[sample].hom_ref[word] |= bit; break;
                    case GenotypeState::HomAlt: active.masks[sample].hom_alt[word] |= bit; break;
                    case GenotypeState::Missing: active.masks[sample].missing[word] |= bit; break;
                    case GenotypeState::Heterozygous: break;
                }
            }
            ++active.variant_count;
            ++assignments_this_position;
        }
        ++in_gene_record_count;
        gene_assignment_count += assignments_this_position;
    }

    for (auto& active : active_genes) finalize_gene(active);
    for (auto& group : groups) {
        group.ndiff_output.open(output_dir / (group.first_sample + ".Ndiff.u16le.bin"),
                                std::ios::binary | std::ios::trunc);
        group.nmiss_output.open(output_dir / (group.first_sample + ".Nmiss.u16le.bin"),
                                std::ios::binary | std::ios::trunc);
        if (!group.ndiff_output || !group.nmiss_output) {
            fail("无法创建染色体基因计数分片: " + group.first_sample);
        }
        group.ndiff_output.write(reinterpret_cast<const char*>(group.ndiff.data()),
                                 static_cast<std::streamsize>(group.ndiff.size() * sizeof(uint16_t)));
        group.nmiss_output.write(reinterpret_cast<const char*>(group.nmiss.data()),
                                 static_cast<std::streamsize>(group.nmiss.size() * sizeof(uint16_t)));
        group.ndiff_output.close();
        group.nmiss_output.close();
        if (!group.ndiff_output || !group.nmiss_output) {
            fail("关闭染色体基因计数分片失败: " + group.first_sample);
        }
    }

    std::ofstream metadata(output_dir / "chromosome.meta.tsv", std::ios::trunc);
    if (!metadata) fail("无法创建染色体基因计数元数据");
    metadata << "CHROM\tGENE_COUNT\tSAMPLE_COUNT\tPAIR_COUNT\tVCF_RECORD_COUNT\tIN_GENE_RECORD_COUNT"
                "\tGENE_ASSIGNMENT_COUNT\tGENES_WITH_VARIANTS\tDUPLICATE_POSITION_COUNT\tMATRIX_MEMORY_BYTES\n";
    metadata << chromosome << '\t' << total_genes << '\t' << samples.size() << '\t' << pair_count
             << '\t' << record_count << '\t' << in_gene_record_count << '\t' << gene_assignment_count
             << '\t' << genes_with_variants << '\t' << duplicate_position_count << '\t' << matrix_bytes << '\n';
    metadata.close();
    log_info("完成基因计数 " + chromosome + ": VCF记录=" + std::to_string(record_count) +
             ", 基因内记录=" + std::to_string(in_gene_record_count) +
             ", 基因归属次数=" + std::to_string(gene_assignment_count) +
             ", 含变异基因=" + std::to_string(genes_with_variants) +
             ", 重复CHROM+POS记录=" + std::to_string(duplicate_position_count));
    return 0;
}

static std::string haplotype_frequency_cell(uint64_t count, uint64_t denominator) {
    if (denominator == 0) return "NA";
    std::ostringstream value;
    value << std::fixed << std::setprecision(2)
          << (100.0 * static_cast<double>(count) / static_cast<double>(denominator)) << '%';
    return value.str();
}

static int haplotype_mode(const Arguments& args) {
    const std::string chromosome = required(args, "--chrom");
    const auto genes = read_gene_bed(required(args, "--gene-bed"));
    const auto samples = read_nonempty_lines(required(args, "--samples"));
    const auto group_names = read_nonempty_lines(required(args, "--group-order"));
    const uint64_t top_number_u64 = positive_u64(required(args, "--top-number"), "--top-number");
    if (top_number_u64 > std::numeric_limits<std::size_t>::max()) fail("--top-number超出系统可表示范围");
    const std::size_t top_number = static_cast<std::size_t>(top_number_u64);
    const fs::path output_path = required(args, "--output");
    if (samples.empty()) fail("单倍型分析没有样本");
    if (group_names.empty()) fail("单倍型分析没有Group");
    for (const auto& gene : genes) {
        if (gene.chrom != chromosome) fail("单倍型基因区间包含非目标染色体: " + gene.chrom);
    }

    std::unordered_map<std::string, std::size_t> group_index;
    for (std::size_t index = 0; index < group_names.size(); ++index) {
        if (!group_index.emplace(group_names[index], index).second) {
            fail("Group顺序列表包含重复值: " + group_names[index]);
        }
    }
    std::ifstream sample_group_input(required(args, "--sample-group"));
    if (!sample_group_input) fail("无法读取单倍型样本分组表");
    std::vector<std::size_t> sample_groups;
    std::string line;
    while (std::getline(sample_group_input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() != 2 || sample_groups.size() >= samples.size() || fields[0] != samples[sample_groups.size()]) {
            fail("单倍型样本分组表必须与--samples的样本名称及顺序完全一致: " + line);
        }
        const auto iterator = group_index.find(fields[1]);
        if (iterator == group_index.end()) fail("单倍型样本所属Group不在Group顺序列表中: " + fields[1]);
        sample_groups.push_back(iterator->second);
    }
    if (sample_groups.size() != samples.size()) fail("单倍型样本分组表未完整覆盖--samples");

    struct ActiveHaplotypeGene {
        std::size_t index = 0;
        std::vector<std::string> codes;
    };
    std::vector<ActiveHaplotypeGene> active_genes;
    std::vector<std::string> result_rows(genes.size());
    std::vector<uint8_t> finalized(genes.size(), 0);
    std::vector<GenotypeState> states(samples.size());
    std::size_t gene_cursor = 0;
    uint64_t record_count = 0;
    uint64_t assignment_count = 0;
    uint64_t last_position = 0;
    bool have_last_position = false;

    auto finalize_gene = [&](std::size_t gene_index, const std::vector<std::string>& codes) {
        if (finalized[gene_index]) fail("单倍型基因被重复完成: " + genes[gene_index].id);
        std::unordered_map<std::string, uint64_t> all_counts;
        std::vector<uint64_t> group_valid(group_names.size(), 0);
        std::vector<std::unordered_map<std::string, uint64_t>> group_counts(group_names.size());
        uint64_t all_valid = 0;
        for (std::size_t sample = 0; sample < samples.size(); ++sample) {
            const std::string code = codes.empty() ? std::string{} : codes[sample];
            ++all_valid;
            ++all_counts[code];
            ++group_valid[sample_groups[sample]];
            ++group_counts[sample_groups[sample]][code];
        }
        std::vector<std::pair<std::string, uint64_t>> global_ranked(all_counts.begin(), all_counts.end());
        std::sort(global_ranked.begin(), global_ranked.end(), [](const auto& first, const auto& second) {
            if (first.second != second.second) return first.second > second.second;
            return first.first < second.first;
        });
        if (global_ranked.size() > top_number) global_ranked.resize(top_number);

        std::ostringstream row;
        const auto& gene = genes[gene_index];
        row << gene.id << '\t' << gene.chrom << '\t' << gene.start << '\t' << gene.end << '\t' << gene.length
            << '\t' << all_valid;
        for (std::size_t rank = 0; rank < top_number; ++rank) {
            row << '\t';
            if (rank >= global_ranked.size()) row << "NA";
            else row << haplotype_frequency_cell(global_ranked[rank].second, all_valid);
        }
        for (std::size_t group = 0; group < group_names.size(); ++group) {
            row << '\t' << group_valid[group];
            for (std::size_t rank = 0; rank < top_number; ++rank) {
                row << '\t';
                if (rank >= global_ranked.size() || group_valid[group] == 0) {
                    row << "NA";
                } else {
                    const auto iterator = group_counts[group].find(global_ranked[rank].first);
                    const uint64_t count = iterator == group_counts[group].end() ? 0 : iterator->second;
                    row << haplotype_frequency_cell(count, group_valid[group]);
                }
            }
        }
        result_rows[gene_index] = row.str();
        finalized[gene_index] = 1;
    };

    while (std::getline(std::cin, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() != samples.size() + 2) {
            fail("bcftools单倍型GT数据列数异常：预期" + std::to_string(samples.size() + 2) +
                 "列，实际" + std::to_string(fields.size()) + "列");
        }
        if (fields[0] != chromosome) fail("单倍型VCF记录染色体不一致: " + fields[0]);
        uint64_t position = 0;
        try {
            std::size_t parsed = 0;
            position = std::stoull(fields[1], &parsed);
            if (parsed != fields[1].size() || position == 0) fail("POS无效");
        } catch (...) {
            fail("单倍型VCF POS无效: " + fields[1]);
        }
        if (have_last_position && position < last_position) {
            fail("单倍型VCF位置未按坐标排序: " + std::to_string(position) +
                 " < " + std::to_string(last_position));
        }
        last_position = position;
        have_last_position = true;
        ++record_count;
        auto iterator = active_genes.begin();
        while (iterator != active_genes.end()) {
            if (genes[iterator->index].end < position) {
                finalize_gene(iterator->index, iterator->codes);
                iterator = active_genes.erase(iterator);
            } else {
                ++iterator;
            }
        }
        while (gene_cursor < genes.size() && genes[gene_cursor].start <= position) {
            if (genes[gene_cursor].end >= position) {
                ActiveHaplotypeGene active;
                active.index = gene_cursor;
                active.codes.resize(samples.size());
                active_genes.push_back(std::move(active));
            }
            ++gene_cursor;
        }
        bool used = false;
        for (const auto& active : active_genes) {
            if (gene_contains_position(genes[active.index], position)) { used = true; break; }
        }
        if (!used) continue;
        for (std::size_t sample = 0; sample < samples.size(); ++sample) states[sample] = classify_gt(fields[sample + 2]);
        for (auto& active : active_genes) {
            if (!gene_contains_position(genes[active.index], position)) continue;
            for (std::size_t sample = 0; sample < samples.size(); ++sample) {
                switch (states[sample]) {
                    case GenotypeState::HomRef: active.codes[sample].push_back('0'); break;
                    case GenotypeState::Heterozygous: active.codes[sample].push_back('N'); break;
                    case GenotypeState::HomAlt: active.codes[sample].push_back('2'); break;
                    case GenotypeState::Missing:
                        active.codes[sample].push_back('N');
                        break;
                }
            }
            ++assignment_count;
        }
    }
    for (auto& active : active_genes) {
        finalize_gene(active.index, active.codes);
    }
    const std::vector<std::string> empty_codes;
    for (std::size_t index = 0; index < genes.size(); ++index) {
        if (!finalized[index]) finalize_gene(index, empty_codes);
    }

    std::ofstream output(output_path, std::ios::trunc);
    if (!output) fail("无法创建单倍型统计表: " + output_path.string());
    output << "GeneID\tChr\tStart\tEnd\tLength\tAllSampleNumber";
    for (std::size_t rank = 1; rank <= top_number; ++rank) output << "\tAllSampleNo" << rank;
    for (const auto& group : group_names) {
        output << '\t' << group << "SampleNumber";
        for (std::size_t rank = 1; rank <= top_number; ++rank) output << '\t' << group << "SampleNo" << rank;
    }
    output << '\n';
    for (const auto& row : result_rows) output << row << '\n';
    output.close();
    if (!output) fail("写入单倍型统计表失败");
    log_info("完成单倍型统计 " + chromosome + ": genes=" + std::to_string(genes.size()) +
             ", samples=" + std::to_string(samples.size()) + ", VCF_records=" + std::to_string(record_count) +
             ", gene_assignments=" + std::to_string(assignment_count) + ", top=" + std::to_string(top_number));
    return 0;
}

static std::pair<std::string, std::string> parse_pair(const std::string& pair) {
    const std::size_t comma = pair.find(',');
    if (comma == std::string::npos || comma == 0 || comma + 1 >= pair.size() ||
        pair.find(',', comma + 1) != std::string::npos) {
        fail("样本对格式无效: " + pair);
    }
    return {pair.substr(0, comma), pair.substr(comma + 1)};
}

struct BinaryGeneCounts {
    std::ifstream ndiff;
    std::ifstream nmiss;
    std::vector<uint16_t> ndiff_row;
    std::vector<uint16_t> nmiss_row;
};

static BinaryGeneCounts open_binary_counts(const std::string& ndiff_path,
                                           const std::string& nmiss_path,
                                           std::size_t gene_count,
                                           std::size_t pair_count) {
    const uint64_t expected_bytes = static_cast<uint64_t>(gene_count) * pair_count * sizeof(uint16_t);
    if (!fs::is_regular_file(ndiff_path) || !fs::is_regular_file(nmiss_path) ||
        fs::file_size(ndiff_path) != expected_bytes || fs::file_size(nmiss_path) != expected_bytes) {
        fail("基因计数二进制文件大小与基因数/样本对数不匹配");
    }
    BinaryGeneCounts counts;
    counts.ndiff.open(ndiff_path, std::ios::binary);
    counts.nmiss.open(nmiss_path, std::ios::binary);
    if (!counts.ndiff || !counts.nmiss) fail("无法读取基因计数二进制文件");
    counts.ndiff_row.resize(pair_count);
    counts.nmiss_row.resize(pair_count);
    return counts;
}

static void read_binary_row(BinaryGeneCounts& counts) {
    const std::streamsize bytes = static_cast<std::streamsize>(counts.ndiff_row.size() * sizeof(uint16_t));
    counts.ndiff.read(reinterpret_cast<char*>(counts.ndiff_row.data()), bytes);
    counts.nmiss.read(reinterpret_cast<char*>(counts.nmiss_row.data()), bytes);
    if (!counts.ndiff || !counts.nmiss) fail("读取基因计数二进制行失败");
}

struct GroupIndexEntry {
    std::string first_sample;
    std::size_t pair_count = 0;
};

static std::vector<GroupIndexEntry> read_group_index(const std::string& path) {
    std::ifstream input(path);
    if (!input) fail("无法读取样本组索引: " + path);
    std::vector<GroupIndexEntry> groups;
    std::unordered_set<std::string> seen;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() < 2 || fields[0].empty() || !seen.insert(fields[0]).second) {
            fail("样本组索引格式无效或第一样本ID重复: " + line);
        }
        const uint64_t count = positive_u64(fields[1], "样本组样本对数");
        if (count > std::numeric_limits<std::size_t>::max()) fail("样本组样本对数超出平台限制");
        groups.push_back({fields[0], static_cast<std::size_t>(count)});
    }
    if (groups.empty()) fail("样本组索引为空: " + path);
    return groups;
}

struct RunningMoments {
    uint64_t count = 0;
    double mean = 0.0;
    double m2 = 0.0;

    void add(double value) {
        ++count;
        const double delta = value - mean;
        mean += delta / static_cast<double>(count);
        m2 += delta * (value - mean);
    }
};

static int compare_groups_mode(const Arguments& args) {
    const auto cache_genes = read_gene_bed(required(args, "--gene-bed"));
    const auto genes = read_gene_bed(required(args, "--selected-gene-bed"));
    const auto cache_groups = read_group_index(required(args, "--group-index"));
    const fs::path pair_group_dir = required(args, "--pair-group-dir");
    const auto group_names = read_nonempty_lines(required(args, "--group-order"));
    const std::string sample_group_path = required(args, "--sample-group");
    const fs::path output_path = required(args, "--output");
    if (group_names.size() < 2) fail("组间比较至少需要两个Group");

    std::unordered_map<std::string, std::size_t> group_index;
    group_index.reserve(group_names.size() * 2);
    for (std::size_t index = 0; index < group_names.size(); ++index) {
        if (!group_index.emplace(group_names[index], index).second) {
            fail("Group顺序列表包含重复值: " + group_names[index]);
        }
    }

    std::ifstream sample_group_input(sample_group_path);
    if (!sample_group_input) fail("无法读取标准化样本分组表: " + sample_group_path);
    std::unordered_map<std::string, std::size_t> sample_to_group;
    std::vector<uint64_t> group_sample_counts(group_names.size(), 0);
    std::string line;
    while (std::getline(sample_group_input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() != 2 || fields[0].empty() || fields[1].empty()) {
            fail("标准化样本分组表必须恰好两列: " + line);
        }
        const auto group_iterator = group_index.find(fields[1]);
        if (group_iterator == group_index.end()) fail("样本所属Group不在Group顺序列表中: " + fields[1]);
        if (!sample_to_group.emplace(fields[0], group_iterator->second).second) {
            fail("标准化样本分组表包含重复SampleID: " + fields[0]);
        }
        ++group_sample_counts[group_iterator->second];
    }
    for (std::size_t group = 0; group < group_names.size(); ++group) {
        if (group_sample_counts[group] < 2) fail("Group少于2个样本: " + group_names[group]);
    }

    std::unordered_map<std::string, std::size_t> cache_gene_rows;
    cache_gene_rows.reserve(cache_genes.size() * 2);
    for (std::size_t row = 0; row < cache_genes.size(); ++row) cache_gene_rows.emplace(cache_genes[row].id, row);
    std::vector<std::size_t> cache_to_selected(cache_genes.size(), std::numeric_limits<std::size_t>::max());
    for (std::size_t selected_row = 0; selected_row < genes.size(); ++selected_row) {
        const auto iterator = cache_gene_rows.find(genes[selected_row].id);
        if (iterator == cache_gene_rows.end()) fail("分析基因未在GeneCounter缓存中找到: " + genes[selected_row].id);
        const std::size_t cache_row = iterator->second;
        const auto& cached = cache_genes[cache_row];
        if (!same_gene_definition(cached, genes[selected_row])) {
            fail("分析基因与GeneCounter缓存坐标不一致: " + genes[selected_row].id);
        }
        cache_to_selected[cache_row] = selected_row;
    }

    if (genes.size() > std::numeric_limits<std::size_t>::max() / group_names.size()) {
        fail("逐基因分组统计数组大小溢出");
    }
    std::vector<RunningMoments> moments(genes.size() * group_names.size());
    std::vector<uint64_t> found_pair_counts(group_names.size(), 0);
    uint64_t valid_cells = 0;
    uint64_t invalid_cells = 0;

    for (const auto& cache_group : cache_groups) {
        const fs::path group_dir = pair_group_dir / cache_group.first_sample;
        const auto group_pairs = read_nonempty_lines((group_dir / (cache_group.first_sample + ".SamplePair.list")).string());
        if (group_pairs.size() != cache_group.pair_count) fail("样本组样本对数与索引不一致: " + cache_group.first_sample);
        std::vector<std::size_t> selected_columns;
        std::vector<std::size_t> selected_column_targets;
        selected_columns.reserve(group_pairs.size());
        selected_column_targets.reserve(group_pairs.size());
        for (std::size_t column = 0; column < group_pairs.size(); ++column) {
            const auto pair = parse_pair(group_pairs[column]);
            const auto first = sample_to_group.find(pair.first);
            const auto second = sample_to_group.find(pair.second);
            if (first == sample_to_group.end() || second == sample_to_group.end()) continue;
            if (first->second == second->second) {
                selected_columns.push_back(column);
                selected_column_targets.push_back(first->second);
                ++found_pair_counts[first->second];
            }
        }
        if (selected_columns.empty()) continue;

        auto counts = open_binary_counts(
            (group_dir / (cache_group.first_sample + ".Ndiff.u16le.bin")).string(),
            (group_dir / (cache_group.first_sample + ".Nmiss.u16le.bin")).string(),
            cache_genes.size(), group_pairs.size());
        for (std::size_t cache_row = 0; cache_row < cache_genes.size(); ++cache_row) {
            read_binary_row(counts);
            const std::size_t selected_row = cache_to_selected[cache_row];
            if (selected_row == std::numeric_limits<std::size_t>::max()) continue;
            const uint64_t length = genes[selected_row].length;
            for (std::size_t selected_column = 0; selected_column < selected_columns.size(); ++selected_column) {
                const std::size_t column = selected_columns[selected_column];
                const uint64_t nmiss = counts.nmiss_row[column];
                if (nmiss >= length) {
                    ++invalid_cells;
                    continue;
                }
                const double value = static_cast<double>(counts.ndiff_row[column]) * 1000.0 /
                                     static_cast<double>(length - nmiss);
                if (!std::isfinite(value) || value < 0.0) fail("逐基因DSR_per_kb计算得到无效值: " + genes[selected_row].id);
                const std::size_t target = selected_column_targets[selected_column];
                moments[selected_row * group_names.size() + target].add(value);
                ++valid_cells;
            }
        }
    }
    for (std::size_t group = 0; group < group_names.size(); ++group) {
        const uint64_t expected_pairs = group_sample_counts[group] * (group_sample_counts[group] - 1) / 2;
        if (found_pair_counts[group] != expected_pairs) {
            fail("Group=" + group_names[group] + "的组内样本对覆盖不完整: expected=" +
                 std::to_string(expected_pairs) + ", found=" + std::to_string(found_pair_counts[group]));
        }
    }
    std::ofstream output(output_path, std::ios::trunc);
    if (!output) fail("无法创建组间逐基因统计表: " + output_path.string());
    output << "GeneID\tChr\tStart\tEnd\tRegionLength";
    for (const auto& group_name : group_names) {
        output << '\t' << group_name << ".ValidPairCount"
               << '\t' << group_name << ".Mean_DSRperKb"
               << '\t' << group_name << ".Variance_DSRperKb";
    }
    output << '\n' << std::setprecision(15);
    for (std::size_t row = 0; row < genes.size(); ++row) {
        const uint64_t length = genes[row].length;
        output << genes[row].id << '\t' << genes[row].chrom << '\t' << genes[row].start << '\t'
               << genes[row].end << '\t' << length;
        for (std::size_t group = 0; group < group_names.size(); ++group) {
            const auto& value = moments[row * group_names.size() + group];
            output << '\t' << value.count << '\t';
            if (value.count == 0) output << "NA"; else output << value.mean;
            output << '\t';
            if (value.count < 2) output << "NA"; else output << value.m2 / static_cast<double>(value.count - 1);
        }
        output << '\n';
    }
    output.close();
    if (!output) fail("写入组间逐基因统计表失败");
    log_info("完成Group内逐基因流式统计: groups=" + std::to_string(group_names.size()) +
             ", genes=" + std::to_string(genes.size()) + ", valid_cells=" + std::to_string(valid_cells) +
             ", invalid_cells=" + std::to_string(invalid_cells));
    return 0;
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) fail("用法: ggComp.plus.counter <extract-regions|count|haplotype|compare-groups> [参数]");
        const std::string mode = argv[1];
        const Arguments args = parse_arguments(argc, argv, 2);
        if (mode == "extract-regions") return extract_regions_mode(args);
        if (mode == "count") return count_mode(args);
        if (mode == "haplotype") return haplotype_mode(args);
        if (mode == "compare-groups") return compare_groups_mode(args);
        fail("未知C++辅助程序模式: " + mode);
    } catch (const std::exception& error) {
        const std::string line = "[ERROR " + timestamp() + "] " + error.what() + "\n";
        (void)::write(STDERR_FILENO, line.data(), line.size());
        return 1;
    }
}
