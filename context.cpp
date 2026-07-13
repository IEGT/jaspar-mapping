#include "context_core.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

constexpr std::size_t DEFAULT_BATCH_SIZE = 32;
constexpr std::size_t MAX_TEMPORARY_BATCH_FILES = 64;
// The spool is deliberately fixed-width so its size can be checked before a
// potentially very wide legacy context table is generated.
constexpr std::uintmax_t SERIALIZED_ANNOTATION_BYTES =
    3 * sizeof(double) + sizeof(std::uint64_t) + 2 * sizeof(std::uint8_t);

struct Options {
    std::int64_t flankBp = 150;
    std::size_t batchSize = DEFAULT_BATCH_SIZE;
    std::filesystem::path temporaryDirectory = std::filesystem::temp_directory_path();
    bool excludeSameNamedLocus = true;
    bool verbose = false;
    std::string mainBedFile;
    std::vector<std::string> otherBedFiles;
};

struct BatchFile {
    std::filesystem::path path;
    std::size_t columnCount = 0;
};

class TemporaryFiles {
public:
    // Remove private spool files on success and on exception; final output is
    // written by the Makefile through a separate atomic temporary file.
    ~TemporaryFiles() {
        for (const BatchFile& file : files_) {
            std::error_code error;
            std::filesystem::remove(file.path, error);
        }
    }

    void add(BatchFile file) { files_.push_back(std::move(file)); }
    const std::vector<BatchFile>& files() const { return files_; }

private:
    std::vector<BatchFile> files_;
};

void printHelp(const char* programName) {
    std::cout
        << "Usage: " << programName
        << " [options] <main.bed> <other1.bed> [<other2.bed> ...]\n\n"
        << "Annotate each anchor in main.bed with the strongest and the number of\n"
        << "neighboring motif occurrences from each other BED file. Distances are\n"
        << "between BED interval centers. The output is a legacy wide table on stdout.\n\n"
        << "Options:\n"
        << "  --flank BP          Center-to-center radius on each side (default: 150)\n"
        << "  --batch-size N      Secondary files held per batch (default: 32)\n"
        << "  --temp-dir DIR      Temporary wide-table spool directory\n"
        << "  --include-self      Include an exact same-named anchor locus (default: exclude)\n"
        << "  -v, --verbose       Report batch progress\n"
        << "  -h, --help          Show this help and exit\n\n"
        << "Inputs must be BED6+ sorted by chromosome and start; main.bed should\n"
        << "contain one fixed-width anchor motif. _Shift is the\n"
        << "motif-center displacement after orienting the anchor to '+'.\n"
        << "_GenomicShift retains the unflipped genomic displacement.\n";
}

std::uint64_t parseUnsigned(const std::string& text, const std::string& option) {
    std::size_t consumed = 0;
    std::uint64_t value = 0;
    try {
        value = std::stoull(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument(option + " expects a non-negative integer");
    }
    if (consumed != text.size()) {
        throw std::invalid_argument(option + " expects a non-negative integer");
    }
    return value;
}

Options parseArguments(const int argc, const char* const argv[]) {
    Options options;
    std::vector<std::string> positional;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        auto requireValue = [&](const std::string& option) -> std::string {
            if (++index >= argc) throw std::invalid_argument(option + " requires a value");
            return argv[index];
        };

        if (argument == "-h" || argument == "--help") {
            printHelp(argv[0]);
            std::exit(0);
        }
        if (argument == "-v" || argument == "--verbose") {
            options.verbose = true;
        } else if (argument == "--include-self") {
            options.excludeSameNamedLocus = false;
        } else if (argument == "--flank") {
            const std::uint64_t value = parseUnsigned(requireValue(argument), argument);
            if (value > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max() / 2)) {
                throw std::invalid_argument("--flank is too large");
            }
            options.flankBp = static_cast<std::int64_t>(value);
        } else if (argument == "--batch-size") {
            const std::uint64_t value = parseUnsigned(requireValue(argument), argument);
            if (value == 0 || value > std::numeric_limits<std::size_t>::max()) {
                throw std::invalid_argument("--batch-size must be greater than zero");
            }
            options.batchSize = static_cast<std::size_t>(value);
        } else if (argument == "--temp-dir") {
            options.temporaryDirectory = requireValue(argument);
        } else if (!argument.empty() && argument.front() == '-') {
            throw std::invalid_argument("unknown option: " + argument);
        } else {
            positional.push_back(argument);
        }
    }

    if (positional.size() < 2) {
        throw std::invalid_argument("a main BED and at least one neighboring BED are required");
    }
    options.mainBedFile = positional.front();
    options.otherBedFiles.assign(positional.begin() + 1, positional.end());
    return options;
}

std::string humanBytes(const std::uintmax_t bytes) {
    static const char* units[] = {"B", "KiB", "MiB", "GiB", "TiB"};
    double value = static_cast<double>(bytes);
    std::size_t unit = 0;
    while (value >= 1024.0 && unit + 1 < std::size(units)) {
        value /= 1024.0;
        ++unit;
    }
    std::ostringstream output;
    output << std::fixed << std::setprecision(unit == 0 ? 0 : 1) << value << ' ' << units[unit];
    return output.str();
}

std::uintmax_t estimatedTemporaryBytes(const std::size_t rows, const std::size_t columns) {
    if (rows != 0 && columns > std::numeric_limits<std::uintmax_t>::max() / rows) {
        throw std::overflow_error("temporary-space estimate overflowed");
    }
    const std::uintmax_t cells = static_cast<std::uintmax_t>(rows) * columns;
    if (cells > std::numeric_limits<std::uintmax_t>::max() / SERIALIZED_ANNOTATION_BYTES) {
        throw std::overflow_error("temporary-space estimate overflowed");
    }
    return cells * SERIALIZED_ANNOTATION_BYTES;
}

void verifyTemporarySpace(const Options& options, const std::size_t rowCount) {
    std::filesystem::create_directories(options.temporaryDirectory);
    const std::uintmax_t estimate = estimatedTemporaryBytes(rowCount, options.otherBedFiles.size());
    const std::filesystem::space_info available = std::filesystem::space(options.temporaryDirectory);
    constexpr std::uintmax_t safetyMargin = 64ULL * 1024ULL * 1024ULL;
    std::cerr << "I: Legacy wide context needs approximately " << humanBytes(estimate)
              << " of temporary space in '" << options.temporaryDirectory.string()
              << "' (available: " << humanBytes(available.available) << ").\n";
    if (estimate > available.available
        || available.available - estimate < safetyMargin) {
        throw std::runtime_error("insufficient temporary disk space for legacy wide context output");
    }
}

std::filesystem::path uniqueBatchPath(const std::filesystem::path& directory,
                                      const std::size_t batchIndex) {
    const auto clockValue = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::random_device random;
    for (int attempt = 0; attempt < 100; ++attempt) {
        const std::string name = "jaspar-context-"
                                 + std::to_string(static_cast<unsigned long long>(clockValue))
                                 + "-" + std::to_string(batchIndex)
                                 + "-" + std::to_string(random()) + ".bin";
        const std::filesystem::path path = directory / name;
        if (!std::filesystem::exists(path)) return path;
    }
    throw std::runtime_error("could not allocate a unique temporary context filename");
}

template <typename Value>
void writeBinary(std::ostream& output, const Value& value) {
    output.write(reinterpret_cast<const char*>(&value), sizeof(value));
    if (!output) throw std::runtime_error("failed writing temporary context data");
}

template <typename Value>
void readBinary(std::istream& input, Value& value) {
    input.read(reinterpret_cast<char*>(&value), sizeof(value));
    if (!input) throw std::runtime_error("failed reading temporary context data");
}

void writeAnnotation(std::ostream& output, const motif_context::ContextAnnotation& annotation) {
    const std::uint8_t hasMatch = annotation.hasMatch ? 1 : 0;
    const std::uint8_t strandEqual = annotation.strandEqual ? 1 : 0;
    writeBinary(output, annotation.genomicCenterShift);
    writeBinary(output, annotation.anchorOrientedCenterShift);
    writeBinary(output, annotation.score);
    writeBinary(output, annotation.count);
    writeBinary(output, hasMatch);
    writeBinary(output, strandEqual);
}

motif_context::ContextAnnotation readAnnotation(std::istream& input) {
    motif_context::ContextAnnotation annotation;
    std::uint8_t hasMatch = 0;
    std::uint8_t strandEqual = 0;
    readBinary(input, annotation.genomicCenterShift);
    readBinary(input, annotation.anchorOrientedCenterShift);
    readBinary(input, annotation.score);
    readBinary(input, annotation.count);
    readBinary(input, hasMatch);
    readBinary(input, strandEqual);
    annotation.hasMatch = hasMatch != 0;
    annotation.strandEqual = strandEqual != 0;
    return annotation;
}

std::vector<std::string> contextColumnPrefixes(const std::vector<std::string>& files) {
    std::vector<std::string> prefixes;
    std::unordered_set<std::string> observed;
    prefixes.reserve(files.size());
    for (const std::string& file : files) {
        const std::string prefix = motif_context::basenameWithoutContextSuffix(file);
        if (!observed.insert(prefix).second) {
            throw std::runtime_error("context output column prefix collision: " + prefix);
        }
        prefixes.push_back(prefix);
    }
    return prefixes;
}

std::size_t effectiveBatchSize(const Options& options) {
    const std::size_t minimumForFileLimit =
        (options.otherBedFiles.size() + MAX_TEMPORARY_BATCH_FILES - 1)
        / MAX_TEMPORARY_BATCH_FILES;
    return std::max(options.batchSize, std::max<std::size_t>(1, minimumForFileLimit));
}

void buildBatchFiles(const Options& options,
                     const std::vector<motif_context::BedEntry>& anchors,
                     TemporaryFiles& temporaryFiles) {
    const std::size_t batchSize = effectiveBatchSize(options);
    std::size_t batchIndex = 0;
    for (std::size_t first = 0; first < options.otherBedFiles.size(); first += batchSize) {
        const std::size_t last = std::min(first + batchSize, options.otherBedFiles.size());
        std::vector<std::vector<motif_context::ContextAnnotation>> columns;
        columns.reserve(last - first);
        for (std::size_t index = first; index < last; ++index) {
            if (options.verbose) {
                std::cerr << "I: Context " << (index + 1) << "/"
                          << options.otherBedFiles.size() << ": "
                          << options.otherBedFiles[index] << '\n';
            }
            columns.push_back(motif_context::annotateFromBedFile(
                anchors, options.otherBedFiles[index], options.flankBp,
                options.excludeSameNamedLocus, false));
        }

        const std::filesystem::path path =
            uniqueBatchPath(options.temporaryDirectory, batchIndex++);
        std::ofstream output(path, std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("failed to create temporary file: " + path.string());
        for (std::size_t row = 0; row < anchors.size(); ++row) {
            for (const auto& column : columns) writeAnnotation(output, column[row]);
        }
        output.close();
        if (!output) throw std::runtime_error("failed to close temporary file: " + path.string());
        temporaryFiles.add({path, last - first});
    }
}

void printDistance(const double value) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) < 1e-12) {
        std::cout << static_cast<std::int64_t>(rounded);
    } else {
        std::cout << std::fixed << std::setprecision(1) << value << std::defaultfloat;
    }
}

void writeWideOutput(const Options& options,
                     const motif_context::BedTable& anchors,
                     const std::vector<std::string>& prefixes,
                     const TemporaryFiles& temporaryFiles) {
    std::vector<std::ifstream> batchInputs;
    batchInputs.reserve(temporaryFiles.files().size());
    for (const BatchFile& file : temporaryFiles.files()) {
        batchInputs.emplace_back(file.path, std::ios::binary);
        if (!batchInputs.back()) {
            throw std::runtime_error("failed opening temporary file: " + file.path.string());
        }
    }

    if (!anchors.header.empty()) {
        std::cout << anchors.header
                  << "\tContextFlankBp\tContextDistanceMetric\tContextSelfMatchPolicy";
        for (const std::string& prefix : prefixes) {
            std::cout << '\t' << prefix << "_Shift"
                      << '\t' << prefix << "_GenomicShift"
                      << '\t' << prefix << "_Score"
                      << '\t' << prefix << "_StrandEqual"
                      << '\t' << prefix << "_NumInWindow";
        }
        std::cout << '\n';
    }

    for (const motif_context::BedEntry& anchor : anchors.entries) {
        std::cout << anchor.line << '\t' << options.flankBp
                  << "\tmotif_center\t"
                  << (options.excludeSameNamedLocus ? "exclude_same_named_locus"
                                                    : "include_self");
        for (std::size_t batch = 0; batch < batchInputs.size(); ++batch) {
            for (std::size_t column = 0;
                 column < temporaryFiles.files()[batch].columnCount; ++column) {
                const motif_context::ContextAnnotation annotation =
                    readAnnotation(batchInputs[batch]);
                if (!annotation.hasMatch) {
                    std::cout << "\tNA\tNA\tNA\tNA\t0";
                    continue;
                }
                std::cout << '\t';
                printDistance(annotation.anchorOrientedCenterShift);
                std::cout << '\t';
                printDistance(annotation.genomicCenterShift);
                std::cout << '\t' << std::setprecision(15) << annotation.score
                          << '\t' << (annotation.strandEqual ? 1 : 0)
                          << '\t' << annotation.count;
            }
        }
        std::cout << '\n';
    }
}

} // namespace

int main(const int argc, const char* const argv[]) {
    try {
        const Options options = parseArguments(argc, argv);
        const std::vector<std::string> prefixes =
            contextColumnPrefixes(options.otherBedFiles);
        const motif_context::BedTable anchors =
            motif_context::readBedTable(options.mainBedFile, options.verbose);
        if (anchors.entries.empty()) {
            throw std::runtime_error("main BED contains no anchor rows");
        }

        verifyTemporarySpace(options, anchors.entries.size());
        TemporaryFiles temporaryFiles;
        buildBatchFiles(options, anchors.entries, temporaryFiles);
        writeWideOutput(options, anchors, prefixes, temporaryFiles);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "E: " << error.what() << '\n';
        std::cerr << "Try --help for usage.\n";
        return 1;
    }
}
