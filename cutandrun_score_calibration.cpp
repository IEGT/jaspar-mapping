#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr std::int64_t NEGATIVE_INFINITY_BIN =
    std::numeric_limits<std::int64_t>::min();
constexpr std::int64_t POSITIVE_INFINITY_BIN =
    std::numeric_limits<std::int64_t>::max();

struct CoverageSpec {
    std::string sampleId;
    std::filesystem::path path;
};

struct Arguments {
    std::filesystem::path plusParquet;
    std::filesystem::path minusParquet;
    std::filesystem::path contextPlusParquet;
    std::filesystem::path contextMinusParquet;
    std::filesystem::path featureParquet;
    std::filesystem::path outputDir;
    std::vector<CoverageSpec> coverage;
    std::string motifId = "MA0861.2";
    std::string scoreMode = "log2_relative_risk";
    std::string chrom = "1";
    double pseudocount = 1.0;
    double binWidth = 0.2;
    std::int64_t motifLength = 16;
    std::string contextMotifId = "MA1961.2";
    double contextPseudocount = 1.0;
    std::int64_t contextMotifLength = 11;
    std::int64_t contextFlank = 150;
    double minimumAnchorScore = -10.0;
    std::uint64_t progressEvery = 10'000'000;

    bool hasContextMotif() const {
        return !contextPlusParquet.empty() || !contextMinusParquet.empty();
    }
};

void printHelp(const char* program) {
    std::cout
        << "Usage: " << program << " --plus-parquet FILE --minus-parquet FILE\n"
        << "       --coverage SAMPLE_ID=FILE [--coverage SAMPLE_ID=FILE ...]\n"
        << "       --output-dir DIR [options]\n\n"
        << "Stream orientation-paired dense Parquet scores against one or more sorted,\n"
        << "non-overlapping positive-depth bedGraph tracks. Only compact score-bin and\n"
        << "threshold summaries are written; no per-position table is materialized.\n\n"
        << "Options:\n"
        << "  --motif-id ID          Motif identifier (default: MA0861.2)\n"
        << "  --motif-length BP      Scored alignment span (default: 16)\n"
        << "  --chrom CHROM          Chromosome, with or without chr prefix (default: 1)\n"
        << "  --score-mode MODE      Score provenance label (default: log2_relative_risk)\n"
        << "  --pseudocount N        Score provenance value (default: 1)\n"
        << "  --bin-width N          Fixed score-bin width (default: 0.2)\n"
        << "  --context-plus-parquet FILE   Optional plus-orientation dense context scores\n"
        << "  --context-minus-parquet FILE  Optional minus-orientation dense context scores\n"
        << "  --context-motif-id ID         Context motif identifier (default: MA1961.2)\n"
        << "  --context-motif-length BP     Context alignment span (default: 11)\n"
        << "  --context-pseudocount N       Context score provenance (default: 1)\n"
        << "  --context-flank BP            Maximum center distance on either side (default: 150)\n"
        << "  --minimum-anchor-score N      Lowest anchor score in joint histogram (default: -10)\n"
        << "  --feature-parquet FILE        Optional exact per-anchor context/coverage table\n"
        << "  --progress-every N     Report every N alignment starts; 0 disables (default: 10000000)\n"
        << "  -h, --help             Show this help\n\n"
        << "Coverage is positive only where bedGraph depth is finite and greater than zero.\n"
        << "Overlapping bedGraph rows are rejected because depth would be ambiguous. Adjacent\n"
        << "rows are joined for immersion. A scored span [start,start+L) is supported only\n"
        << "when component_start < start and component_end > start+L. Its effective depth\n"
        << "is then the maximum bedGraph depth across that span; otherwise it is zero.\n";
}

std::string requireValue(int& index, int argc, char** argv,
                         const std::string& option) {
    if (++index >= argc) {
        throw std::runtime_error(option + " requires a value");
    }
    return argv[index];
}

std::int64_t parsePositiveInteger(const std::string& value,
                                  const std::string& option) {
    std::size_t parsed = 0;
    const long long result = std::stoll(value, &parsed);
    if (parsed != value.size() || result <= 0) {
        throw std::runtime_error(option + " expects a positive integer, got '" +
                                 value + "'");
    }
    return result;
}

std::uint64_t parseNonnegativeInteger(const std::string& value,
                                      const std::string& option) {
    if (value.empty() || value.front() == '-') {
        throw std::runtime_error(option + " expects a non-negative integer, got '" +
                                 value + "'");
    }
    std::size_t parsed = 0;
    const unsigned long long result = std::stoull(value, &parsed);
    if (parsed != value.size()) {
        throw std::runtime_error(option + " expects a non-negative integer, got '" +
                                 value + "'");
    }
    return result;
}

double parseFiniteDouble(const std::string& value, const std::string& option) {
    std::size_t parsed = 0;
    const double result = std::stod(value, &parsed);
    if (parsed != value.size() || !std::isfinite(result)) {
        throw std::runtime_error(option + " expects a finite number, got '" + value +
                                 "'");
    }
    return result;
}

Arguments parseArguments(int argc, char** argv) {
    Arguments arguments;
    for (int index = 1; index < argc; ++index) {
        const std::string option = argv[index];
        if (option == "-h" || option == "--help") {
            printHelp(argv[0]);
            std::exit(0);
        } else if (option == "--plus-parquet") {
            arguments.plusParquet = requireValue(index, argc, argv, option);
        } else if (option == "--minus-parquet") {
            arguments.minusParquet = requireValue(index, argc, argv, option);
        } else if (option == "--context-plus-parquet") {
            arguments.contextPlusParquet = requireValue(index, argc, argv, option);
        } else if (option == "--context-minus-parquet") {
            arguments.contextMinusParquet = requireValue(index, argc, argv, option);
        } else if (option == "--feature-parquet") {
            arguments.featureParquet = requireValue(index, argc, argv, option);
        } else if (option == "--output-dir") {
            arguments.outputDir = requireValue(index, argc, argv, option);
        } else if (option == "--coverage") {
            const std::string specification = requireValue(index, argc, argv, option);
            const std::size_t separator = specification.find('=');
            if (separator == std::string::npos || separator == 0 ||
                separator + 1 == specification.size()) {
                throw std::runtime_error(
                    "--coverage expects SAMPLE_ID=FILE, got '" + specification + "'");
            }
            arguments.coverage.push_back(
                {specification.substr(0, separator), specification.substr(separator + 1)});
        } else if (option == "--motif-id") {
            arguments.motifId = requireValue(index, argc, argv, option);
        } else if (option == "--motif-length") {
            arguments.motifLength = parsePositiveInteger(
                requireValue(index, argc, argv, option), option);
        } else if (option == "--context-motif-id") {
            arguments.contextMotifId = requireValue(index, argc, argv, option);
        } else if (option == "--context-motif-length") {
            arguments.contextMotifLength = parsePositiveInteger(
                requireValue(index, argc, argv, option), option);
        } else if (option == "--context-flank") {
            arguments.contextFlank = parsePositiveInteger(
                requireValue(index, argc, argv, option), option);
        } else if (option == "--chrom") {
            arguments.chrom = requireValue(index, argc, argv, option);
        } else if (option == "--score-mode") {
            arguments.scoreMode = requireValue(index, argc, argv, option);
        } else if (option == "--pseudocount") {
            arguments.pseudocount = parseFiniteDouble(
                requireValue(index, argc, argv, option), option);
        } else if (option == "--context-pseudocount") {
            arguments.contextPseudocount = parseFiniteDouble(
                requireValue(index, argc, argv, option), option);
        } else if (option == "--minimum-anchor-score") {
            arguments.minimumAnchorScore = parseFiniteDouble(
                requireValue(index, argc, argv, option), option);
        } else if (option == "--bin-width") {
            arguments.binWidth = parseFiniteDouble(
                requireValue(index, argc, argv, option), option);
            if (arguments.binWidth <= 0.0) {
                throw std::runtime_error("--bin-width must be greater than zero");
            }
        } else if (option == "--progress-every") {
            arguments.progressEvery = parseNonnegativeInteger(
                requireValue(index, argc, argv, option), option);
        } else {
            throw std::runtime_error("unknown argument: " + option);
        }
    }

    if (arguments.plusParquet.empty() || arguments.minusParquet.empty() ||
        arguments.outputDir.empty() || arguments.coverage.empty()) {
        throw std::runtime_error(
            "--plus-parquet, --minus-parquet, --coverage, and --output-dir are required");
    }
    if (arguments.contextPlusParquet.empty() !=
        arguments.contextMinusParquet.empty()) {
        throw std::runtime_error(
            "--context-plus-parquet and --context-minus-parquet must be supplied together");
    }
    if (!arguments.featureParquet.empty() && !arguments.hasContextMotif()) {
        throw std::runtime_error(
            "--feature-parquet requires both context Parquet inputs");
    }
    std::set<std::string> sampleIds;
    for (const CoverageSpec& coverage : arguments.coverage) {
        if (!sampleIds.insert(coverage.sampleId).second) {
            throw std::runtime_error("duplicate coverage sample ID: " + coverage.sampleId);
        }
        if (coverage.sampleId.find_first_of("\t\r\n") != std::string::npos) {
            throw std::runtime_error("coverage sample IDs may not contain tabs or newlines");
        }
    }
    return arguments;
}

std::string normalizedChrom(std::string chrom) {
    if (chrom.size() >= 3 &&
        (chrom[0] == 'c' || chrom[0] == 'C') &&
        (chrom[1] == 'h' || chrom[1] == 'H') &&
        (chrom[2] == 'r' || chrom[2] == 'R')) {
        chrom.erase(0, 3);
    }
    return chrom;
}

struct CoverageRun {
    std::int64_t start;
    std::int64_t end;
    double depth;
};

struct CoverageComponent {
    std::int64_t start;
    std::int64_t end;
};

struct Observation {
    bool supported = false;
    double depth = 0.0;
    std::size_t componentIndex = 0;
};

struct TrackBinStats {
    std::uint64_t supported = 0;
    long double depthSum = 0.0;
    long double log1pDepthSum = 0.0;
    double maximumDepth = 0.0;
};

class CoverageTrack {
public:
    explicit CoverageTrack(CoverageSpec specification)
        : specification_(std::move(specification)) {}

    void load(const std::string& requestedChrom) {
        std::ifstream input(specification_.path);
        if (!input) {
            throw std::runtime_error("cannot open coverage bedGraph: " +
                                     specification_.path.string());
        }

        const std::string targetChrom = normalizedChrom(requestedChrom);
        std::string line;
        std::uint64_t lineNumber = 0;
        while (std::getline(input, line)) {
            ++lineNumber;
            if (line.empty() || line[0] == '#' || line.rfind("track", 0) == 0 ||
                line.rfind("browser", 0) == 0) {
                continue;
            }
            std::istringstream fields(line);
            std::string chrom;
            std::int64_t start = 0;
            std::int64_t end = 0;
            double depth = 0.0;
            if (!(fields >> chrom >> start >> end >> depth)) {
                throw std::runtime_error("invalid bedGraph row at " +
                                         specification_.path.string() + ":" +
                                         std::to_string(lineNumber));
            }
            if (normalizedChrom(chrom) != targetChrom) {
                continue;
            }
            if (start < 0 || end <= start || !std::isfinite(depth) || depth <= 0.0) {
                throw std::runtime_error("invalid positive-depth bedGraph interval at " +
                                         specification_.path.string() + ":" +
                                         std::to_string(lineNumber));
            }
            if (!runs_.empty()) {
                const CoverageRun& previous = runs_.back();
                if (start < previous.start) {
                    throw std::runtime_error("bedGraph is not coordinate-sorted at " +
                                             specification_.path.string() + ":" +
                                             std::to_string(lineNumber));
                }
                if (start < previous.end) {
                    throw std::runtime_error(
                        "overlapping bedGraph rows are not supported at " +
                        specification_.path.string() + ":" +
                        std::to_string(lineNumber));
                }
            }
            runs_.push_back({start, end, depth});
        }
        if (runs_.empty()) {
            throw std::runtime_error("coverage bedGraph has no positive intervals on " +
                                     requestedChrom + ": " +
                                     specification_.path.string());
        }

        for (const CoverageRun& run : runs_) {
            if (components_.empty() || run.start > components_.back().end) {
                components_.push_back({run.start, run.end});
            } else {
                components_.back().end = std::max(components_.back().end, run.end);
            }
        }
        componentBestBin_.resize(components_.size(), NEGATIVE_INFINITY_BIN);
        componentHasScore_.resize(components_.size(), false);
    }

    Observation observe(const std::int64_t motifStart,
                        const std::int64_t motifLength) {
        const std::int64_t motifEnd = motifStart + motifLength;

        while (nextRunToAdd_ < runs_.size() &&
               runs_[nextRunToAdd_].start < motifEnd) {
            activeDepths_.insert(runs_[nextRunToAdd_].depth);
            ++nextRunToAdd_;
        }
        while (nextRunToRemove_ < nextRunToAdd_ &&
               runs_[nextRunToRemove_].end <= motifStart) {
            const auto depth = activeDepths_.find(runs_[nextRunToRemove_].depth);
            if (depth == activeDepths_.end()) {
                throw std::runtime_error("internal coverage-depth state became inconsistent");
            }
            activeDepths_.erase(depth);
            ++nextRunToRemove_;
        }
        while (componentCursor_ < components_.size() &&
               components_[componentCursor_].end <= motifStart) {
            ++componentCursor_;
        }

        if (componentCursor_ >= components_.size()) {
            return {};
        }
        const CoverageComponent& component = components_[componentCursor_];
        if (component.start < motifStart && component.end > motifEnd) {
            if (activeDepths_.empty()) {
                throw std::runtime_error(
                    "strictly immersed motif unexpectedly has no active depth");
            }
            return {true, *activeDepths_.rbegin(), componentCursor_};
        }
        return {};
    }

    void recordComponentBest(const std::size_t componentIndex,
                             const std::int64_t binIndex) {
        if (!componentHasScore_[componentIndex] ||
            binIndex > componentBestBin_[componentIndex]) {
            componentBestBin_[componentIndex] = binIndex;
            componentHasScore_[componentIndex] = true;
        }
    }

    const CoverageSpec& specification() const { return specification_; }
    std::size_t intervalCount() const { return runs_.size(); }
    std::size_t componentCount() const { return components_.size(); }

    std::uint64_t coveredBases() const {
        std::uint64_t result = 0;
        for (const CoverageComponent& component : components_) {
            result += static_cast<std::uint64_t>(component.end - component.start);
        }
        return result;
    }

    std::uint64_t possibleImmersedStarts(const std::int64_t motifLength) const {
        std::uint64_t result = 0;
        for (const CoverageComponent& component : components_) {
            const std::int64_t count = component.end - component.start - motifLength - 1;
            if (count > 0) {
                result += static_cast<std::uint64_t>(count);
            }
        }
        return result;
    }

    std::map<std::int64_t, std::uint64_t> componentBestHistogram() const {
        std::map<std::int64_t, std::uint64_t> histogram;
        for (std::size_t index = 0; index < componentBestBin_.size(); ++index) {
            if (componentHasScore_[index]) {
                ++histogram[componentBestBin_[index]];
            }
        }
        return histogram;
    }

private:
    CoverageSpec specification_;
    std::vector<CoverageRun> runs_;
    std::vector<CoverageComponent> components_;
    std::size_t nextRunToAdd_ = 0;
    std::size_t nextRunToRemove_ = 0;
    std::size_t componentCursor_ = 0;
    std::multiset<double> activeDepths_;
    std::vector<std::int64_t> componentBestBin_;
    std::vector<bool> componentHasScore_;
};

class DenseParquetReader {
public:
    explicit DenseParquetReader(const std::filesystem::path& path) {
        parquet::arrow::FileReaderBuilder builder;
        const arrow::Status openStatus = builder.OpenFile(path.string(), false);
        if (!openStatus.ok()) {
            throw std::runtime_error("opening Parquet input " + path.string() + ": " +
                                     openStatus.ToString());
        }
        auto readerResult = builder.Build();
        if (!readerResult.ok()) {
            throw std::runtime_error("building Parquet reader for " + path.string() +
                                     ": " + readerResult.status().ToString());
        }
        fileReader_ = std::move(readerResult).ValueOrDie();
        auto batchResult = fileReader_->GetRecordBatchReader();
        if (!batchResult.ok()) {
            throw std::runtime_error("creating Parquet batch reader for " +
                                     path.string() + ": " +
                                     batchResult.status().ToString());
        }
        batchReader_ = std::move(batchResult).ValueOrDie();
    }

    std::shared_ptr<arrow::RecordBatch> next() {
        std::shared_ptr<arrow::RecordBatch> batch;
        const arrow::Status status = batchReader_->ReadNext(&batch);
        if (!status.ok()) {
            throw std::runtime_error("reading dense Parquet batch: " + status.ToString());
        }
        return batch;
    }

private:
    std::unique_ptr<parquet::arrow::FileReader> fileReader_;
    std::unique_ptr<arrow::RecordBatchReader> batchReader_;
};

struct Histogram {
    explicit Histogram(const std::size_t trackCount) : byTrack(trackCount) {}

    std::size_t slotFor(const std::int64_t binIndex) {
        const auto existing = slotByBin.find(binIndex);
        if (existing != slotByBin.end()) {
            return existing->second;
        }
        const std::size_t slot = binIndices.size();
        slotByBin.emplace(binIndex, slot);
        binIndices.push_back(binIndex);
        total.push_back(0);
        for (auto& trackBins : byTrack) {
            trackBins.emplace_back();
        }
        return slot;
    }

    std::unordered_map<std::int64_t, std::size_t> slotByBin;
    std::vector<std::int64_t> binIndices;
    std::vector<std::uint64_t> total;
    std::vector<std::vector<TrackBinStats>> byTrack;
};

struct JointKey {
    std::int64_t anchorBin;
    std::int64_t contextBin;

    bool operator==(const JointKey& other) const {
        return anchorBin == other.anchorBin && contextBin == other.contextBin;
    }
};

struct JointKeyHash {
    std::size_t operator()(const JointKey& key) const {
        const std::size_t first = std::hash<std::int64_t>{}(key.anchorBin);
        const std::size_t second = std::hash<std::int64_t>{}(key.contextBin);
        return first ^ (second + 0x9e3779b9U + (first << 6U) + (first >> 2U));
    }
};

struct JointGeometryStats {
    std::uint64_t sameOrientation = 0;
    std::uint64_t oppositeOrientation = 0;
    long double signedCenterDistanceTwiceSum = 0.0;
    long double absoluteCenterDistanceTwiceSum = 0.0;
};

struct JointTrackBinStats {
    TrackBinStats all;
    std::uint64_t supportedSameOrientation = 0;
    std::uint64_t supportedOppositeOrientation = 0;
};

struct JointHistogram {
    explicit JointHistogram(const std::size_t trackCount) : byTrack(trackCount) {}

    std::size_t slotFor(const JointKey key) {
        const auto existing = slotByKey.find(key);
        if (existing != slotByKey.end()) {
            return existing->second;
        }
        const std::size_t slot = keys.size();
        slotByKey.emplace(key, slot);
        keys.push_back(key);
        total.push_back(0);
        geometry.emplace_back();
        for (auto& trackBins : byTrack) {
            trackBins.emplace_back();
        }
        return slot;
    }

    std::unordered_map<JointKey, std::size_t, JointKeyHash> slotByKey;
    std::vector<JointKey> keys;
    std::vector<std::uint64_t> total;
    std::vector<JointGeometryStats> geometry;
    std::vector<std::vector<JointTrackBinStats>> byTrack;
};

std::int64_t scoreBin(const double score, const double binWidth) {
    if (std::isinf(score)) {
        return score < 0.0 ? NEGATIVE_INFINITY_BIN : POSITIVE_INFINITY_BIN;
    }
    if (!std::isfinite(score)) {
        throw std::runtime_error("dense score contains NaN");
    }
    const double bin = std::floor(score / binWidth);
    if (bin <= static_cast<double>(NEGATIVE_INFINITY_BIN + 1) ||
        bin >= static_cast<double>(POSITIVE_INFINITY_BIN - 1)) {
        throw std::runtime_error("score bin exceeds signed 64-bit range");
    }
    return static_cast<std::int64_t>(bin);
}

std::string thresholdText(const std::int64_t binIndex, const double binWidth) {
    if (binIndex == NEGATIVE_INFINITY_BIN) {
        return "-Inf";
    }
    if (binIndex == POSITIVE_INFINITY_BIN) {
        return "Inf";
    }
    std::ostringstream output;
    output << std::setprecision(12) << binIndex * binWidth;
    return output.str();
}

struct DenseColumns {
    std::shared_ptr<arrow::Int64Array> starts;
    std::shared_ptr<arrow::ListArray> lists;
    std::shared_ptr<arrow::FloatArray> values;
};

DenseColumns denseColumns(const std::shared_ptr<arrow::RecordBatch>& batch) {
    const int startIndex = batch->schema()->GetFieldIndex("block_start");
    const int scoreIndex = batch->schema()->GetFieldIndex("scores");
    if (startIndex < 0 || scoreIndex < 0) {
        throw std::runtime_error("dense Parquet must contain block_start and scores");
    }
    if (batch->column(startIndex)->type_id() != arrow::Type::INT64 ||
        batch->column(scoreIndex)->type_id() != arrow::Type::LIST) {
        throw std::runtime_error("dense Parquet has unexpected column types");
    }
    auto starts = std::static_pointer_cast<arrow::Int64Array>(batch->column(startIndex));
    auto lists = std::static_pointer_cast<arrow::ListArray>(batch->column(scoreIndex));
    if (lists->values()->type_id() != arrow::Type::FLOAT) {
        throw std::runtime_error("dense Parquet scores must be float32 lists");
    }
    return {starts, lists,
            std::static_pointer_cast<arrow::FloatArray>(lists->values())};
}

struct DenseScore {
    std::int64_t start = 0;
    double score = 0.0;
    char strand = '.';
    bool valid = false;
};

class PairedDenseScoreStream {
public:
    PairedDenseScoreStream(const std::filesystem::path& plusPath,
                           const std::filesystem::path& minusPath)
        : plusReader_(plusPath), minusReader_(minusPath) {}

    bool next(DenseScore& result) {
        while (true) {
            if (!plusBatch_ || row_ >= plusBatch_->num_rows()) {
                if (!loadBatch()) {
                    return false;
                }
            }
            if (!rowReady_) {
                prepareRow();
            }
            if (offset_ >= rowLength_) {
                ++row_;
                rowReady_ = false;
                continue;
            }

            result.start = blockStart_ + offset_;
            if (result.start <= lastStart_) {
                throw std::runtime_error(
                    "dense Parquet scores are not strictly coordinate-ordered");
            }
            const bool plusValid = !plusColumns_.values->IsNull(plusOffset_ + offset_);
            const bool minusValid = !minusColumns_.values->IsNull(minusOffset_ + offset_);
            result.valid = plusValid || minusValid;
            if (!result.valid) {
                result.score = 0.0;
                result.strand = '.';
            } else if (plusValid &&
                       (!minusValid || plusColumns_.values->Value(plusOffset_ + offset_) >=
                                           minusColumns_.values->Value(minusOffset_ + offset_))) {
                result.score = plusColumns_.values->Value(plusOffset_ + offset_);
                result.strand = '+';
            } else {
                result.score = minusColumns_.values->Value(minusOffset_ + offset_);
                result.strand = '-';
            }
            lastStart_ = result.start;
            ++offset_;
            return true;
        }
    }

private:
    bool loadBatch() {
        while (true) {
            plusBatch_ = plusReader_.next();
            minusBatch_ = minusReader_.next();
            if (!plusBatch_ || !minusBatch_) {
                if (plusBatch_ || minusBatch_) {
                    throw std::runtime_error(
                        "plus and minus Parquet streams have different batch counts");
                }
                return false;
            }
            if (plusBatch_->num_rows() != minusBatch_->num_rows()) {
                throw std::runtime_error(
                    "plus and minus Parquet batches have different row counts");
            }
            plusColumns_ = denseColumns(plusBatch_);
            minusColumns_ = denseColumns(minusBatch_);
            row_ = 0;
            rowReady_ = false;
            if (plusBatch_->num_rows() > 0) {
                return true;
            }
        }
    }

    void prepareRow() {
        if (plusColumns_.starts->IsNull(row_) || minusColumns_.starts->IsNull(row_) ||
            plusColumns_.lists->IsNull(row_) || minusColumns_.lists->IsNull(row_)) {
            throw std::runtime_error("dense Parquet contains a null block");
        }
        blockStart_ = plusColumns_.starts->Value(row_);
        if (blockStart_ != minusColumns_.starts->Value(row_)) {
            throw std::runtime_error("plus and minus Parquet block starts do not match");
        }
        rowLength_ = plusColumns_.lists->value_length(row_);
        if (rowLength_ != minusColumns_.lists->value_length(row_)) {
            throw std::runtime_error(
                "plus and minus Parquet score-list lengths do not match");
        }
        if (rowLength_ <= 0 || blockStart_ <= lastStart_) {
            throw std::runtime_error(
                "dense Parquet contains an empty or out-of-order score block");
        }
        plusOffset_ = plusColumns_.lists->value_offset(row_);
        minusOffset_ = minusColumns_.lists->value_offset(row_);
        offset_ = 0;
        rowReady_ = true;
    }

    DenseParquetReader plusReader_;
    DenseParquetReader minusReader_;
    std::shared_ptr<arrow::RecordBatch> plusBatch_;
    std::shared_ptr<arrow::RecordBatch> minusBatch_;
    DenseColumns plusColumns_;
    DenseColumns minusColumns_;
    std::int64_t row_ = 0;
    std::int64_t blockStart_ = 0;
    std::int64_t rowLength_ = 0;
    std::int64_t plusOffset_ = 0;
    std::int64_t minusOffset_ = 0;
    std::int64_t offset_ = 0;
    std::int64_t lastStart_ = -1;
    bool rowReady_ = false;
};

struct ContextCandidate {
    DenseScore score;
    std::int64_t centerTwice;
};

void requireArrowStatus(const arrow::Status& status, const std::string& operation) {
    if (!status.ok()) {
        throw std::runtime_error(operation + ": " + status.ToString());
    }
}

std::string parquetColumnSuffix(const std::string& sampleId) {
    std::string result;
    result.reserve(sampleId.size() + 7);
    for (const unsigned char character : sampleId) {
        result.push_back(std::isalnum(character) || character == '_'
                             ? static_cast<char>(character)
                             : '_');
    }
    if (result.empty() || std::isdigit(static_cast<unsigned char>(result.front()))) {
        result.insert(0, "sample_");
    }
    return result;
}

class AnchorFeatureParquetWriter {
public:
    AnchorFeatureParquetWriter(const std::filesystem::path& path,
                               const std::vector<CoverageTrack>& tracks,
                               const std::size_t rowsPerBatch = 65'536)
        : rowsPerBatch_(std::max<std::size_t>(1, rowsPerBatch)) {
        std::vector<std::shared_ptr<arrow::Field>> fields = {
            arrow::field("chrom", arrow::utf8(), false),
            arrow::field("anchor_start", arrow::int64(), false),
            arrow::field("anchor_end", arrow::int64(), false),
            arrow::field("anchor_score", arrow::float32(), false),
            arrow::field("anchor_strand", arrow::int8(), false),
            arrow::field("context_start", arrow::int64(), false),
            arrow::field("context_end", arrow::int64(), false),
            arrow::field("context_score", arrow::float32(), false),
            arrow::field("context_strand", arrow::int8(), false),
            arrow::field("context_center_distance_bp", arrow::float32(), false),
            arrow::field("same_orientation", arrow::boolean(), false)
        };

        std::set<std::string> suffixes;
        for (const CoverageTrack& track : tracks) {
            const std::string suffix =
                parquetColumnSuffix(track.specification().sampleId);
            if (!suffixes.insert(suffix).second) {
                throw std::runtime_error(
                    "coverage sample IDs collide as Parquet column names: " + suffix);
            }
            sampleSuffixes_.push_back(suffix);
            fields.push_back(arrow::field("supported_" + suffix,
                                          arrow::boolean(), false));
            fields.push_back(arrow::field("depth_" + suffix,
                                          arrow::float32(), false));
        }
        schema_ = arrow::schema(std::move(fields));

        auto outputResult = arrow::io::FileOutputStream::Open(path.string());
        if (!outputResult.ok()) {
            throw std::runtime_error("opening feature Parquet output " +
                                     path.string() + ": " +
                                     outputResult.status().ToString());
        }
        outputStream_ = std::move(outputResult).ValueOrDie();

        parquet::WriterProperties::Builder writerPropertiesBuilder;
        writerPropertiesBuilder.compression(parquet::Compression::ZSTD);
        const std::shared_ptr<parquet::WriterProperties> writerProperties =
            writerPropertiesBuilder.build();
        parquet::ArrowWriterProperties::Builder arrowPropertiesBuilder;
        arrowPropertiesBuilder.store_schema();
        const std::shared_ptr<parquet::ArrowWriterProperties> arrowProperties =
            arrowPropertiesBuilder.build();
        auto writerResult = parquet::arrow::FileWriter::Open(
            *schema_, arrow::default_memory_pool(), outputStream_,
            writerProperties, arrowProperties);
        if (!writerResult.ok()) {
            throw std::runtime_error("opening feature Parquet writer: " +
                                     writerResult.status().ToString());
        }
        writer_ = std::move(writerResult).ValueOrDie();
        resetBuilders();
    }

    void append(const Arguments& arguments,
                const DenseScore& anchor,
                const ContextCandidate& context,
                const std::vector<Observation>& observations) {
        if (observations.size() != supportBuilders_.size()) {
            throw std::runtime_error(
                "internal feature observation count does not match coverage tracks");
        }
        const std::int64_t anchorCenterTwice =
            2 * anchor.start + arguments.motifLength;
        const std::int64_t centerDistanceTwice =
            context.centerTwice - anchorCenterTwice;

        requireArrowStatus(chromBuilder_->Append(arguments.chrom),
                           "appending feature chromosome");
        requireArrowStatus(anchorStartBuilder_->Append(anchor.start),
                           "appending feature anchor start");
        requireArrowStatus(anchorEndBuilder_->Append(
                               anchor.start + arguments.motifLength),
                           "appending feature anchor end");
        requireArrowStatus(anchorScoreBuilder_->Append(
                               static_cast<float>(anchor.score)),
                           "appending feature anchor score");
        requireArrowStatus(anchorStrandBuilder_->Append(strandCode(anchor.strand)),
                           "appending feature anchor strand");
        requireArrowStatus(contextStartBuilder_->Append(context.score.start),
                           "appending feature context start");
        requireArrowStatus(contextEndBuilder_->Append(
                               context.score.start + arguments.contextMotifLength),
                           "appending feature context end");
        requireArrowStatus(contextScoreBuilder_->Append(
                               static_cast<float>(context.score.score)),
                           "appending feature context score");
        requireArrowStatus(contextStrandBuilder_->Append(
                               strandCode(context.score.strand)),
                           "appending feature context strand");
        requireArrowStatus(contextDistanceBuilder_->Append(
                               static_cast<float>(centerDistanceTwice) / 2.0F),
                           "appending feature context distance");
        requireArrowStatus(sameOrientationBuilder_->Append(
                               anchor.strand == context.score.strand),
                           "appending feature orientation relation");
        for (std::size_t index = 0; index < observations.size(); ++index) {
            requireArrowStatus(supportBuilders_[index]->Append(
                                   observations[index].supported),
                               "appending feature coverage support");
            requireArrowStatus(depthBuilders_[index]->Append(
                                   static_cast<float>(observations[index].depth)),
                               "appending feature coverage depth");
        }

        ++pendingRows_;
        ++rowsWritten_;
        if (pendingRows_ >= rowsPerBatch_) {
            flush();
        }
    }

    void close() {
        if (writer_) {
            flush();
            requireArrowStatus(writer_->Close(), "closing feature Parquet writer");
            writer_.reset();
        }
        if (outputStream_) {
            requireArrowStatus(outputStream_->Close(),
                               "closing feature Parquet output");
            outputStream_.reset();
        }
    }

    std::uint64_t rowsWritten() const { return rowsWritten_; }

private:
    static std::int8_t strandCode(const char strand) {
        if (strand == '+') {
            return 1;
        }
        if (strand == '-') {
            return -1;
        }
        return 0;
    }

    void resetBuilders() {
        arrow::MemoryPool* const pool = arrow::default_memory_pool();
        chromBuilder_ = std::make_unique<arrow::StringBuilder>(pool);
        anchorStartBuilder_ = std::make_unique<arrow::Int64Builder>(pool);
        anchorEndBuilder_ = std::make_unique<arrow::Int64Builder>(pool);
        anchorScoreBuilder_ = std::make_unique<arrow::FloatBuilder>(pool);
        anchorStrandBuilder_ = std::make_unique<arrow::Int8Builder>(pool);
        contextStartBuilder_ = std::make_unique<arrow::Int64Builder>(pool);
        contextEndBuilder_ = std::make_unique<arrow::Int64Builder>(pool);
        contextScoreBuilder_ = std::make_unique<arrow::FloatBuilder>(pool);
        contextStrandBuilder_ = std::make_unique<arrow::Int8Builder>(pool);
        contextDistanceBuilder_ = std::make_unique<arrow::FloatBuilder>(pool);
        sameOrientationBuilder_ = std::make_unique<arrow::BooleanBuilder>(pool);
        supportBuilders_.clear();
        depthBuilders_.clear();
        for (std::size_t index = 0; index < sampleSuffixes_.size(); ++index) {
            supportBuilders_.push_back(
                std::make_unique<arrow::BooleanBuilder>(pool));
            depthBuilders_.push_back(std::make_unique<arrow::FloatBuilder>(pool));
        }
        pendingRows_ = 0;
    }

    void flush() {
        if (pendingRows_ == 0) {
            return;
        }
        std::vector<std::shared_ptr<arrow::Array>> columns;
        columns.reserve(11 + 2 * sampleSuffixes_.size());
        finish(*chromBuilder_, columns, "chromosome");
        finish(*anchorStartBuilder_, columns, "anchor start");
        finish(*anchorEndBuilder_, columns, "anchor end");
        finish(*anchorScoreBuilder_, columns, "anchor score");
        finish(*anchorStrandBuilder_, columns, "anchor strand");
        finish(*contextStartBuilder_, columns, "context start");
        finish(*contextEndBuilder_, columns, "context end");
        finish(*contextScoreBuilder_, columns, "context score");
        finish(*contextStrandBuilder_, columns, "context strand");
        finish(*contextDistanceBuilder_, columns, "context distance");
        finish(*sameOrientationBuilder_, columns, "orientation relation");
        for (std::size_t index = 0; index < sampleSuffixes_.size(); ++index) {
            finish(*supportBuilders_[index], columns, "coverage support");
            finish(*depthBuilders_[index], columns, "coverage depth");
        }
        const std::shared_ptr<arrow::RecordBatch> batch =
            arrow::RecordBatch::Make(schema_, static_cast<std::int64_t>(pendingRows_),
                                     std::move(columns));
        requireArrowStatus(writer_->WriteRecordBatch(*batch),
                           "writing feature Parquet record batch");
        resetBuilders();
    }

    template <typename Builder>
    static void finish(Builder& builder,
                       std::vector<std::shared_ptr<arrow::Array>>& columns,
                       const std::string& label) {
        std::shared_ptr<arrow::Array> array;
        requireArrowStatus(builder.Finish(&array),
                           "finishing feature " + label + " array");
        columns.push_back(std::move(array));
    }

    const std::size_t rowsPerBatch_;
    std::size_t pendingRows_ = 0;
    std::uint64_t rowsWritten_ = 0;
    std::vector<std::string> sampleSuffixes_;
    std::shared_ptr<arrow::Schema> schema_;
    std::shared_ptr<arrow::io::FileOutputStream> outputStream_;
    std::unique_ptr<parquet::arrow::FileWriter> writer_;
    std::unique_ptr<arrow::StringBuilder> chromBuilder_;
    std::unique_ptr<arrow::Int64Builder> anchorStartBuilder_;
    std::unique_ptr<arrow::Int64Builder> anchorEndBuilder_;
    std::unique_ptr<arrow::FloatBuilder> anchorScoreBuilder_;
    std::unique_ptr<arrow::Int8Builder> anchorStrandBuilder_;
    std::unique_ptr<arrow::Int64Builder> contextStartBuilder_;
    std::unique_ptr<arrow::Int64Builder> contextEndBuilder_;
    std::unique_ptr<arrow::FloatBuilder> contextScoreBuilder_;
    std::unique_ptr<arrow::Int8Builder> contextStrandBuilder_;
    std::unique_ptr<arrow::FloatBuilder> contextDistanceBuilder_;
    std::unique_ptr<arrow::BooleanBuilder> sameOrientationBuilder_;
    std::vector<std::unique_ptr<arrow::BooleanBuilder>> supportBuilders_;
    std::vector<std::unique_ptr<arrow::FloatBuilder>> depthBuilders_;
};

void processDenseScores(const Arguments& arguments,
                        std::vector<CoverageTrack>& tracks,
                        Histogram& histogram,
                        JointHistogram* jointHistogram,
                        AnchorFeatureParquetWriter* featureWriter,
                        std::uint64_t& alignmentStarts,
                        std::uint64_t& validScores,
                        std::uint64_t& jointValidScores) {
    PairedDenseScoreStream anchorStream(arguments.plusParquet, arguments.minusParquet);
    std::unique_ptr<PairedDenseScoreStream> contextStream;
    std::deque<ContextCandidate> contextMaximum;
    DenseScore nextContext;
    bool contextAvailable = false;
    if (jointHistogram != nullptr) {
        contextStream = std::make_unique<PairedDenseScoreStream>(
            arguments.contextPlusParquet, arguments.contextMinusParquet);
        contextAvailable = contextStream->next(nextContext);
    }
    std::uint64_t nextProgress = arguments.progressEvery;
    std::vector<Observation> observations(tracks.size());

    DenseScore anchor;
    while (anchorStream.next(anchor)) {
        ++alignmentStarts;

        const ContextCandidate* bestContext = nullptr;
        if (jointHistogram != nullptr) {
            const std::int64_t anchorCenterTwice =
                2 * anchor.start + arguments.motifLength;
            const std::int64_t maximumCenterTwice =
                anchorCenterTwice + 2 * arguments.contextFlank;
            while (contextAvailable) {
                const std::int64_t contextCenterTwice =
                    2 * nextContext.start + arguments.contextMotifLength;
                if (contextCenterTwice > maximumCenterTwice) {
                    break;
                }
                if (nextContext.valid) {
                    while (!contextMaximum.empty() &&
                           contextMaximum.back().score.score < nextContext.score) {
                        contextMaximum.pop_back();
                    }
                    contextMaximum.push_back({nextContext, contextCenterTwice});
                }
                contextAvailable = contextStream->next(nextContext);
            }
            const std::int64_t minimumCenterTwice =
                anchorCenterTwice - 2 * arguments.contextFlank;
            while (!contextMaximum.empty() &&
                   contextMaximum.front().centerTwice < minimumCenterTwice) {
                contextMaximum.pop_front();
            }
            if (!contextMaximum.empty()) {
                bestContext = &contextMaximum.front();
            }
        }

        if (!anchor.valid) {
            continue;
        }
        const std::int64_t binIndex = scoreBin(anchor.score, arguments.binWidth);
        const std::size_t slot = histogram.slotFor(binIndex);
        ++histogram.total[slot];
        ++validScores;

        std::size_t jointSlot = 0;
        bool useJointSlot = false;
        bool sameOrientation = false;
        if (jointHistogram != nullptr && bestContext != nullptr &&
            anchor.score >= arguments.minimumAnchorScore) {
            const std::int64_t contextBin =
                scoreBin(bestContext->score.score, arguments.binWidth);
            jointSlot = jointHistogram->slotFor({binIndex, contextBin});
            ++jointHistogram->total[jointSlot];
            ++jointValidScores;
            sameOrientation = anchor.strand == bestContext->score.strand;
            JointGeometryStats& geometry = jointHistogram->geometry[jointSlot];
            if (sameOrientation) {
                ++geometry.sameOrientation;
            } else {
                ++geometry.oppositeOrientation;
            }
            const std::int64_t centerDistanceTwice =
                bestContext->centerTwice - (2 * anchor.start + arguments.motifLength);
            geometry.signedCenterDistanceTwiceSum += centerDistanceTwice;
            geometry.absoluteCenterDistanceTwiceSum +=
                std::abs(centerDistanceTwice);
            useJointSlot = true;
        }

        for (std::size_t trackIndex = 0; trackIndex < tracks.size(); ++trackIndex) {
            const Observation observation =
                tracks[trackIndex].observe(anchor.start, arguments.motifLength);
            observations[trackIndex] = observation;
            if (!observation.supported) {
                continue;
            }
            TrackBinStats& stats = histogram.byTrack[trackIndex][slot];
            ++stats.supported;
            stats.depthSum += observation.depth;
            stats.log1pDepthSum += std::log1p(observation.depth);
            stats.maximumDepth = std::max(stats.maximumDepth, observation.depth);
            tracks[trackIndex].recordComponentBest(observation.componentIndex, binIndex);

            if (useJointSlot) {
                JointTrackBinStats& jointStats =
                    jointHistogram->byTrack[trackIndex][jointSlot];
                ++jointStats.all.supported;
                jointStats.all.depthSum += observation.depth;
                jointStats.all.log1pDepthSum += std::log1p(observation.depth);
                jointStats.all.maximumDepth =
                    std::max(jointStats.all.maximumDepth, observation.depth);
                if (sameOrientation) {
                    ++jointStats.supportedSameOrientation;
                } else {
                    ++jointStats.supportedOppositeOrientation;
                }
            }
        }
        if (featureWriter != nullptr && useJointSlot) {
            featureWriter->append(arguments, anchor, *bestContext, observations);
        }

        if (arguments.progressEvery > 0 && alignmentStarts >= nextProgress) {
            std::cerr << "I: processed " << alignmentStarts
                      << " alignment starts; valid scores=" << validScores;
            if (jointHistogram != nullptr) {
                std::cerr << "; joint candidates=" << jointValidScores;
            }
            std::cerr << '\n';
            while (nextProgress <= alignmentStarts) {
                nextProgress += arguments.progressEvery;
            }
        }
    }
}

double ratio(const long double numerator, const long double denominator) {
    return denominator == 0.0L
        ? std::numeric_limits<double>::quiet_NaN()
        : static_cast<double>(numerator / denominator);
}

struct SummaryRow {
    double rocAuc = std::numeric_limits<double>::quiet_NaN();
    double averagePrecision = std::numeric_limits<double>::quiet_NaN();
    std::string youdenThreshold;
    double youden = -std::numeric_limits<double>::infinity();
    std::string f1Threshold;
    double f1 = -std::numeric_limits<double>::infinity();
};

void writeOutputs(const Arguments& arguments,
                  const std::vector<CoverageTrack>& tracks,
                  const Histogram& histogram,
                  const std::uint64_t alignmentStarts,
                  const std::uint64_t validScores) {
    const std::filesystem::path histogramPath =
        arguments.outputDir / "score_histogram.tsv";
    const std::filesystem::path curvePath =
        arguments.outputDir / "threshold_curve.tsv";
    const std::filesystem::path summaryPath = arguments.outputDir / "summary.tsv";
    const std::filesystem::path configPath = arguments.outputDir / "run_config.tsv";

    std::ofstream histogramOutput(histogramPath);
    std::ofstream curveOutput(curvePath);
    std::ofstream summaryOutput(summaryPath);
    std::ofstream configOutput(configPath);
    if (!histogramOutput || !curveOutput || !summaryOutput || !configOutput) {
        throw std::runtime_error("cannot create calibration output files");
    }

    std::vector<std::size_t> order(histogram.binIndices.size());
    for (std::size_t index = 0; index < order.size(); ++index) {
        order[index] = index;
    }
    std::sort(order.begin(), order.end(), [&](const std::size_t left,
                                               const std::size_t right) {
        return histogram.binIndices[left] > histogram.binIndices[right];
    });

    histogramOutput
        << "sample_id\tmotif_id\tscore_mode\tpseudocount\tbin_width\tbin_index"
        << "\tthreshold\tn_total\tn_supported\teffective_depth_sum"
        << "\tlog1p_effective_depth_sum\tmaximum_effective_depth\n";
    curveOutput
        << "sample_id\tmotif_id\tscore_mode\tpseudocount\tbin_width\tbin_index"
        << "\tthreshold\tselected_motifs\tsupported_selected_motifs"
        << "\tunsupported_selected_motifs\tmotif_precision\tmotif_recall"
        << "\tmotif_false_positive_rate\tcoverage_component_recall"
        << "\tmean_effective_depth\tmean_supported_depth"
        << "\tmean_log1p_effective_depth\tsupport_prevalence"
        << "\tprecision_enrichment\tf1\tyouden_j\n";
    summaryOutput
        << "sample_id\tmotif_id\tscore_mode\tpseudocount\tbin_width"
        << "\tn_alignment_starts\tn_valid_motifs\tn_supported_motifs"
        << "\tsupport_prevalence\tn_coverage_intervals\tn_coverage_components"
        << "\tcovered_bases\tpossible_immersed_starts\tn_score_bins"
        << "\troc_auc_binned\taverage_precision\tyouden_threshold\tyouden_j"
        << "\tf1_threshold\tf1\n";
    configOutput << "sample_id\tcoverage_bedgraph\tchrom\tmotif_id\tmotif_length"
                 << "\tscore_mode\tpseudocount\tbin_width\tplus_parquet"
                 << "\tminus_parquet\tcoverage_merge_rule\timmersion_rule"
                 << "\tdepth_rule\torientation_rule\n";

    const std::uint64_t totalMotifs = validScores;
    for (std::size_t trackIndex = 0; trackIndex < tracks.size(); ++trackIndex) {
        const CoverageTrack& track = tracks[trackIndex];
        std::uint64_t totalSupported = 0;
        std::uint64_t totalUnsupported = 0;
        long double aucNumerator = 0.0;
        std::uint64_t lowerUnsupported = 0;

        for (auto iterator = order.rbegin(); iterator != order.rend(); ++iterator) {
            const std::size_t slot = *iterator;
            const TrackBinStats& stats = histogram.byTrack[trackIndex][slot];
            const std::uint64_t unsupported = histogram.total[slot] - stats.supported;
            aucNumerator += static_cast<long double>(stats.supported) *
                (lowerUnsupported + 0.5L * unsupported);
            lowerUnsupported += unsupported;
            totalSupported += stats.supported;
            totalUnsupported += unsupported;
        }
        const double prevalence = ratio(totalSupported, totalMotifs);
        SummaryRow summary;
        if (totalSupported > 0 && totalUnsupported > 0) {
            summary.rocAuc = ratio(
                aucNumerator,
                static_cast<long double>(totalSupported) * totalUnsupported);
        }

        const std::map<std::int64_t, std::uint64_t> componentBest =
            track.componentBestHistogram();
        std::uint64_t selected = 0;
        std::uint64_t supportedSelected = 0;
        std::uint64_t selectedComponents = 0;
        long double selectedDepth = 0.0;
        long double selectedLogDepth = 0.0;
        double averagePrecision = 0.0;

        for (const std::size_t slot : order) {
            const std::int64_t binIndex = histogram.binIndices[slot];
            const std::string threshold = thresholdText(binIndex, arguments.binWidth);
            const TrackBinStats& stats = histogram.byTrack[trackIndex][slot];
            selected += histogram.total[slot];
            supportedSelected += stats.supported;
            selectedDepth += stats.depthSum;
            selectedLogDepth += stats.log1pDepthSum;
            const auto componentCount = componentBest.find(binIndex);
            if (componentCount != componentBest.end()) {
                selectedComponents += componentCount->second;
            }

            const std::uint64_t unsupportedSelected = selected - supportedSelected;
            const double precision = ratio(supportedSelected, selected);
            const double recall = ratio(supportedSelected, totalSupported);
            const double falsePositiveRate = ratio(unsupportedSelected,
                                                   totalUnsupported);
            const double componentRecall = ratio(selectedComponents,
                                                 track.componentCount());
            const double meanDepth = ratio(selectedDepth, selected);
            const double meanSupportedDepth = ratio(selectedDepth,
                                                    supportedSelected);
            const double meanLogDepth = ratio(selectedLogDepth, selected);
            const double enrichment = prevalence == 0.0
                ? std::numeric_limits<double>::quiet_NaN()
                : precision / prevalence;
            const double f1 = (precision + recall) == 0.0
                ? std::numeric_limits<double>::quiet_NaN()
                : 2.0 * precision * recall / (precision + recall);
            const double youden = recall - falsePositiveRate;

            if (totalSupported > 0) {
                averagePrecision +=
                    static_cast<double>(stats.supported) / totalSupported * precision;
            }
            if (std::isfinite(youden) && youden > summary.youden) {
                summary.youden = youden;
                summary.youdenThreshold = threshold;
            }
            if (std::isfinite(f1) && f1 > summary.f1) {
                summary.f1 = f1;
                summary.f1Threshold = threshold;
            }

            histogramOutput << track.specification().sampleId << '\t'
                << arguments.motifId << '\t' << arguments.scoreMode << '\t'
                << std::setprecision(12) << arguments.pseudocount << '\t'
                << arguments.binWidth << '\t' << binIndex << '\t' << threshold
                << '\t' << histogram.total[slot] << '\t' << stats.supported << '\t'
                << static_cast<double>(stats.depthSum) << '\t'
                << static_cast<double>(stats.log1pDepthSum) << '\t'
                << stats.maximumDepth << '\n';

            curveOutput << track.specification().sampleId << '\t'
                << arguments.motifId << '\t' << arguments.scoreMode << '\t'
                << std::setprecision(12) << arguments.pseudocount << '\t'
                << arguments.binWidth << '\t' << binIndex << '\t' << threshold
                << '\t' << selected << '\t' << supportedSelected << '\t'
                << unsupportedSelected << '\t' << precision << '\t' << recall
                << '\t' << falsePositiveRate << '\t' << componentRecall << '\t'
                << meanDepth << '\t' << meanSupportedDepth << '\t' << meanLogDepth
                << '\t' << prevalence << '\t' << enrichment << '\t' << f1
                << '\t' << youden << '\n';
        }
        summary.averagePrecision = averagePrecision;

        summaryOutput << track.specification().sampleId << '\t' << arguments.motifId
            << '\t' << arguments.scoreMode << '\t' << std::setprecision(12)
            << arguments.pseudocount << '\t' << arguments.binWidth << '\t'
            << alignmentStarts << '\t' << totalMotifs << '\t' << totalSupported
            << '\t' << prevalence << '\t' << track.intervalCount() << '\t'
            << track.componentCount() << '\t' << track.coveredBases() << '\t'
            << track.possibleImmersedStarts(arguments.motifLength) << '\t'
            << histogram.binIndices.size() << '\t' << summary.rocAuc << '\t'
            << summary.averagePrecision << '\t' << summary.youdenThreshold << '\t'
            << summary.youden << '\t' << summary.f1Threshold << '\t' << summary.f1
            << '\n';

        configOutput << track.specification().sampleId << '\t'
            << std::filesystem::absolute(track.specification().path).string() << '\t'
            << arguments.chrom << '\t' << arguments.motifId << '\t'
            << arguments.motifLength << '\t' << arguments.scoreMode << '\t'
            << arguments.pseudocount << '\t' << arguments.binWidth << '\t'
            << std::filesystem::absolute(arguments.plusParquet).string() << '\t'
            << std::filesystem::absolute(arguments.minusParquet).string() << '\t'
            << "overlap_or_adjacent_union\tcomponent_start < motif_start AND "
               "component_end > motif_end\tmaximum_bedgraph_depth_within_span_or_zero"
            << "\tmaximum_orientation_score\n";
    }
}

void writeJointOutputs(const Arguments& arguments,
                       const std::vector<CoverageTrack>& tracks,
                       const JointHistogram& histogram,
                       const std::uint64_t jointValidScores) {
    const std::filesystem::path histogramPath =
        arguments.outputDir / "joint_score_histogram.tsv";
    const std::filesystem::path configPath =
        arguments.outputDir / "joint_run_config.tsv";
    std::ofstream output(histogramPath);
    std::ofstream config(configPath);
    if (!output || !config) {
        throw std::runtime_error("cannot create joint calibration output files");
    }

    output
        << "sample_id\tanchor_motif_id\tcontext_motif_id\tscore_mode"
        << "\tanchor_pseudocount\tcontext_pseudocount\tbin_width"
        << "\tanchor_bin_index\tanchor_threshold\tcontext_bin_index"
        << "\tcontext_threshold\tn_total\tn_supported\teffective_depth_sum"
        << "\tlog1p_effective_depth_sum\tmaximum_effective_depth"
        << "\tn_same_orientation\tn_opposite_orientation"
        << "\tn_supported_same_orientation\tn_supported_opposite_orientation"
        << "\tmean_signed_center_distance_bp\tmean_absolute_center_distance_bp\n";
    config
        << "sample_id\tcoverage_bedgraph\tchrom\tanchor_motif_id"
        << "\tanchor_motif_length\tcontext_motif_id\tcontext_motif_length"
        << "\tscore_mode\tanchor_pseudocount\tcontext_pseudocount\tbin_width"
        << "\tcontext_flank_bp\tminimum_anchor_score\tn_joint_candidates"
        << "\tanchor_plus_parquet\tanchor_minus_parquet\tcontext_plus_parquet"
        << "\tcontext_minus_parquet\tfeature_parquet\tcontext_rule"
        << "\tcontext_tie_rule\n";

    std::vector<std::size_t> order(histogram.keys.size());
    for (std::size_t index = 0; index < order.size(); ++index) {
        order[index] = index;
    }
    std::sort(order.begin(), order.end(), [&](const std::size_t left,
                                               const std::size_t right) {
        const JointKey& leftKey = histogram.keys[left];
        const JointKey& rightKey = histogram.keys[right];
        if (leftKey.anchorBin != rightKey.anchorBin) {
            return leftKey.anchorBin > rightKey.anchorBin;
        }
        return leftKey.contextBin > rightKey.contextBin;
    });

    for (std::size_t trackIndex = 0; trackIndex < tracks.size(); ++trackIndex) {
        const CoverageTrack& track = tracks[trackIndex];
        for (const std::size_t slot : order) {
            const JointKey key = histogram.keys[slot];
            const JointTrackBinStats& stats = histogram.byTrack[trackIndex][slot];
            const JointGeometryStats& geometry = histogram.geometry[slot];
            const double meanSignedDistance =
                ratio(geometry.signedCenterDistanceTwiceSum,
                      2.0L * histogram.total[slot]);
            const double meanAbsoluteDistance =
                ratio(geometry.absoluteCenterDistanceTwiceSum,
                      2.0L * histogram.total[slot]);
            output << track.specification().sampleId << '\t'
                   << arguments.motifId << '\t' << arguments.contextMotifId << '\t'
                   << arguments.scoreMode << '\t' << std::setprecision(12)
                   << arguments.pseudocount << '\t' << arguments.contextPseudocount
                   << '\t' << arguments.binWidth << '\t' << key.anchorBin << '\t'
                   << thresholdText(key.anchorBin, arguments.binWidth) << '\t'
                   << key.contextBin << '\t'
                   << thresholdText(key.contextBin, arguments.binWidth) << '\t'
                   << histogram.total[slot] << '\t' << stats.all.supported << '\t'
                   << static_cast<double>(stats.all.depthSum) << '\t'
                   << static_cast<double>(stats.all.log1pDepthSum) << '\t'
                   << stats.all.maximumDepth << '\t' << geometry.sameOrientation << '\t'
                   << geometry.oppositeOrientation << '\t'
                   << stats.supportedSameOrientation << '\t'
                   << stats.supportedOppositeOrientation << '\t'
                   << meanSignedDistance << '\t' << meanAbsoluteDistance << '\n';
        }

        config << track.specification().sampleId << '\t'
               << std::filesystem::absolute(track.specification().path).string() << '\t'
               << arguments.chrom << '\t' << arguments.motifId << '\t'
               << arguments.motifLength << '\t' << arguments.contextMotifId << '\t'
               << arguments.contextMotifLength << '\t' << arguments.scoreMode << '\t'
               << arguments.pseudocount << '\t' << arguments.contextPseudocount << '\t'
               << arguments.binWidth << '\t' << arguments.contextFlank << '\t'
               << arguments.minimumAnchorScore << '\t' << jointValidScores << '\t'
               << std::filesystem::absolute(arguments.plusParquet).string() << '\t'
               << std::filesystem::absolute(arguments.minusParquet).string() << '\t'
               << std::filesystem::absolute(arguments.contextPlusParquet).string() << '\t'
               << std::filesystem::absolute(arguments.contextMinusParquet).string() << '\t'
               << (arguments.featureParquet.empty()
                       ? std::string()
                       : std::filesystem::absolute(arguments.featureParquet).string())
               << '\t'
               << "maximum orientation-collapsed context score within center distance"
               << "\tearliest genomic start then plus orientation\n";
    }
}

void prepareOutputDirectory(const Arguments& arguments) {
    std::vector<std::string> expected = {
        "score_histogram.tsv", "threshold_curve.tsv", "summary.tsv",
        "run_config.tsv"};
    if (arguments.hasContextMotif()) {
        expected.push_back("joint_score_histogram.tsv");
        expected.push_back("joint_run_config.tsv");
    }
    const std::filesystem::path featurePath = arguments.featureParquet.empty()
        ? std::filesystem::path()
        : std::filesystem::absolute(arguments.featureParquet).lexically_normal();
    for (const std::string& name : expected) {
        const std::filesystem::path path = arguments.outputDir / name;
        if (!featurePath.empty() &&
            featurePath == std::filesystem::absolute(path).lexically_normal()) {
            throw std::runtime_error(
                "feature Parquet path conflicts with standard output: " +
                path.string());
        }
        if (std::filesystem::exists(path)) {
            throw std::runtime_error("refusing to replace existing output: " +
                                     path.string());
        }
    }
    if (!arguments.featureParquet.empty()) {
        if (std::filesystem::exists(arguments.featureParquet)) {
            throw std::runtime_error("refusing to replace existing output: " +
                                     arguments.featureParquet.string());
        }
        const std::filesystem::path parent = arguments.featureParquet.parent_path();
        if (!parent.empty()) {
            std::filesystem::create_directories(parent);
        }
    }
    std::filesystem::create_directories(arguments.outputDir);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Arguments arguments = parseArguments(argc, argv);
        if (!std::filesystem::is_regular_file(arguments.plusParquet)) {
            throw std::runtime_error("plus Parquet file not found: " +
                                     arguments.plusParquet.string());
        }
        if (!std::filesystem::is_regular_file(arguments.minusParquet)) {
            throw std::runtime_error("minus Parquet file not found: " +
                                     arguments.minusParquet.string());
        }
        if (arguments.hasContextMotif()) {
            if (!std::filesystem::is_regular_file(arguments.contextPlusParquet)) {
                throw std::runtime_error("context plus Parquet file not found: " +
                                         arguments.contextPlusParquet.string());
            }
            if (!std::filesystem::is_regular_file(arguments.contextMinusParquet)) {
                throw std::runtime_error("context minus Parquet file not found: " +
                                         arguments.contextMinusParquet.string());
            }
        }
        prepareOutputDirectory(arguments);

        std::vector<CoverageTrack> tracks;
        tracks.reserve(arguments.coverage.size());
        for (const CoverageSpec& specification : arguments.coverage) {
            if (!std::filesystem::is_regular_file(specification.path)) {
                throw std::runtime_error("coverage bedGraph file not found: " +
                                         specification.path.string());
            }
            tracks.emplace_back(specification);
            tracks.back().load(arguments.chrom);
            std::cerr << "I: " << specification.sampleId << ": loaded "
                      << tracks.back().intervalCount() << " intervals in "
                      << tracks.back().componentCount() << " components\n";
        }

        Histogram histogram(tracks.size());
        std::unique_ptr<JointHistogram> jointHistogram;
        std::unique_ptr<AnchorFeatureParquetWriter> featureWriter;
        if (arguments.hasContextMotif()) {
            jointHistogram = std::make_unique<JointHistogram>(tracks.size());
        }
        if (!arguments.featureParquet.empty()) {
            featureWriter = std::make_unique<AnchorFeatureParquetWriter>(
                arguments.featureParquet, tracks);
        }
        std::uint64_t alignmentStarts = 0;
        std::uint64_t validScores = 0;
        std::uint64_t jointValidScores = 0;
        processDenseScores(arguments, tracks, histogram, jointHistogram.get(),
                           featureWriter.get(),
                           alignmentStarts, validScores, jointValidScores);
        if (featureWriter != nullptr) {
            featureWriter->close();
            if (featureWriter->rowsWritten() != jointValidScores) {
                throw std::runtime_error(
                    "feature Parquet row count does not match joint candidates");
            }
        }
        writeOutputs(arguments, tracks, histogram, alignmentStarts, validScores);
        if (jointHistogram != nullptr) {
            writeJointOutputs(arguments, tracks, *jointHistogram, jointValidScores);
        }
        std::cerr << "I: completed " << alignmentStarts << " alignment starts; "
                  << validScores << " had a score in " << histogram.binIndices.size()
                  << " bins";
        if (jointHistogram != nullptr) {
            std::cerr << "; " << jointValidScores << " joint candidates in "
                      << jointHistogram->keys.size() << " score-bin pairs";
        }
        if (featureWriter != nullptr) {
            std::cerr << "; " << featureWriter->rowsWritten()
                      << " exact feature rows written";
        }
        std::cerr << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "E: " << error.what() << '\n';
        return 1;
    }
}
