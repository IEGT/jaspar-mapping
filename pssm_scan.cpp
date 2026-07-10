#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <filesystem> // test and creating for directories, defines filesystem::path
#include <algorithm> // For removing brackets, replace
#include <array>
#include <iomanip>   // For formatted output / setting precision
#include <getopt.h> // For argument handling, defines getopt_long and optarg
#include <cctype>
#include <cmath>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <csignal>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <utility>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

#ifdef PSSM_SCAN_WITH_PARQUET
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>
#endif

#include "progress.h"
#include "pssm.h"
#include "pssm_scan_core.h"
#include "compressed_file_reader.h"

namespace {

typedef std::set<std::string> chromosome_set_type;

int beVerbose = 0;
int showDebug = 0;
#ifdef _WIN32
volatile LONG progressStatusRequested = 0;
#else
volatile std::sig_atomic_t progressStatusRequested = 0;
#endif
static constexpr size_t DEFAULT_SCORE_BLOCK_SIZE = 65536;

#ifdef _WIN32
BOOL WINAPI handleProgressControlEvent(const DWORD controlType) {
    if (controlType != CTRL_BREAK_EVENT) {
        return FALSE;
    }

    InterlockedExchange(&progressStatusRequested, 1);
    return TRUE;
}
#else
void handleProgressSignal(int) {
    progressStatusRequested = 1;
}

void installProgressSignalHandler(const int signalNumber, const char* signalName) {
    struct sigaction action;
    std::memset(&action, 0, sizeof(action));
    action.sa_handler = handleProgressSignal;
    sigemptyset(&action.sa_mask);
#ifdef SA_RESTART
    action.sa_flags = SA_RESTART;
#endif

    if (sigaction(signalNumber, &action, nullptr) != 0) {
        std::cerr << "W: Failed to install " << signalName << " progress handler: "
                  << std::strerror(errno) << std::endl;
    }
}
#endif

void installProgressRequestHandlers() {
#ifdef _WIN32
    if (SetConsoleCtrlHandler(handleProgressControlEvent, TRUE) == 0) {
        std::cerr << "W: Failed to install CTRL+Break progress handler: Windows error "
                  << GetLastError() << std::endl;
    }
#else
    installProgressSignalHandler(SIGUSR1, "SIGUSR1");
#ifdef SIGINFO
    installProgressSignalHandler(SIGINFO, "SIGINFO");
#endif
#endif
}

bool consumeProgressStatusRequest() {
#ifdef _WIN32
    return InterlockedExchange(&progressStatusRequested, 0) != 0;
#else
    if (progressStatusRequested == 0) {
        return false;
    }
    progressStatusRequested = 0;
    return true;
#endif
}

// global variable to control verbosity
std::unordered_map<char, const double> backgroundFrequencies = {
    {'A', 0.25},  // Assuming equal background probabilities; adjust as needed
    {'C', 0.25},
    {'G', 0.25},
    {'T', 0.25},
    {'a', 0.25},
    {'c', 0.25},
    {'g', 0.25},
    {'t', 0.25},
    {'N', 1},
    {'n', 1}
};

// Region constraining the scan for motifs
class Region {
    public:
        std::string chromosome;
        long from;
        long to;
        std::string name;

        static std::vector<Region> parseRegionsFile(const std::string& regionsFile);
};

enum class StrandMode {
    Plus,
    Minus,
    Both
};

struct EncodedChromosome {
    std::string sequence;
    std::vector<std::uint8_t> plusCodes;
    std::vector<std::uint8_t> minusCodes;
};

struct FastaIndexEntry {
    std::string name;
    std::uint64_t length = 0;
    std::uint64_t offset = 0;
    std::uint64_t basesPerLine = 0;
    std::uint64_t bytesPerLine = 0;
};

struct GzipIndexEntry {
    std::uint64_t compressedOffset = 0;
    std::uint64_t uncompressedOffset = 0;
};

std::vector<Region> Region::parseRegionsFile(const std::string& regionsFile) {
    std::ifstream inFile(regionsFile);
    std::vector<Region> regions;

    if (!inFile.is_open()) {
        std::cerr << "E: Error opening regions file: " << regionsFile << std::endl;
        return regions;
    }

    std::string line;
    bool headerSkipped = false;
    while (std::getline(inFile, line)) {
        line = PSSM::trim(line);
        if (line.empty()) continue;

        // Optionally skip a header
        if (!headerSkipped && (line.find("Chromosome") != std::string::npos || line.find("chromosome") != std::string::npos)) {
            headerSkipped = true;
            continue;
        }

        std::stringstream ss(line);
        Region region;

        if (!(ss >> region.chromosome >> region.from >> region.to >> region.name)) {
            std::cerr << "Error parsing line in regions file: " << line << std::endl;
            continue;
        }

        regions.push_back(region);
    }

    inFile.close();
    return regions;
}

std::string strandModeName(StrandMode strandMode) {
    switch (strandMode) {
        case StrandMode::Plus: return "+";
        case StrandMode::Minus: return "-";
        case StrandMode::Both: return "both";
    }
    return "both";
}

std::string strandModeFileLabel(StrandMode strandMode) {
    switch (strandMode) {
        case StrandMode::Plus: return "positive";
        case StrandMode::Minus: return "negative";
        case StrandMode::Both: return "both";
    }
    return "both";
}

std::string coordinateModeName(CoordinateMode coordinateMode) {
    switch (coordinateMode) {
        case CoordinateMode::Legacy: return "legacy";
        case CoordinateMode::Bed: return "bed";
    }
    return "legacy";
}

bool parseStrandMode(const std::string& value, StrandMode& strandMode) {
    if (value == "+" || value == "plus" || value == "positive") {
        strandMode = StrandMode::Plus;
        return true;
    }
    if (value == "-" || value == "minus" || value == "negative") {
        strandMode = StrandMode::Minus;
        return true;
    }
    if (value == "both" || value == "+-" || value == "-+") {
        strandMode = StrandMode::Both;
        return true;
    }
    return false;
}

bool parseCoordinateMode(const std::string& value, CoordinateMode& coordinateMode) {
    std::string normalized = value;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    if (normalized == "legacy") {
        coordinateMode = CoordinateMode::Legacy;
        return true;
    }
    if (normalized == "bed") {
        coordinateMode = CoordinateMode::Bed;
        return true;
    }
    return false;
}

bool parseLongStrict(const char* text, long& value) {
    errno = 0;
    char* end = nullptr;
    const long parsed = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || (end != nullptr && *end != '\0')) {
        return false;
    }
    value = parsed;
    return true;
}

EncodedChromosome encodeChromosome(std::string sequence, bool includeMinusStrand) {
    EncodedChromosome encoded;
    encoded.sequence = std::move(sequence);
    encoded.plusCodes.reserve(encoded.sequence.size());
    for (char base : encoded.sequence) {
        encoded.plusCodes.push_back(codeForBase(base));
    }

    if (includeMinusStrand) {
        encoded.minusCodes.resize(encoded.plusCodes.size());
        for (size_t i = 0; i < encoded.plusCodes.size(); ++i) {
            encoded.minusCodes[encoded.plusCodes.size() - 1 - i] = complementCode(encoded.plusCodes[i]);
        }
    }

    return encoded;
}

std::string sequenceWindowForOutput(const std::string& sequence, size_t start, size_t motifLength, bool reverseComplementWindow) {
    if (!reverseComplementWindow) {
        return sequence.substr(start, motifLength);
    }
    std::string window;
    window.reserve(motifLength);
    for (size_t i = 0; i < motifLength; ++i) {
        window.push_back(complementBase(sequence[start + motifLength - 1 - i]));
    }
    return window;
}

std::string formatOptionalDouble(double value) {
    if (!std::isfinite(value)) {
        return "NA";
    }
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6) << value;
    return ss.str();
}

std::string formatScoreBoundForHelp(double value) {
    if (!std::isfinite(value)) {
        return "unbounded";
    }
    return formatOptionalDouble(value);
}

std::string formatPseudocountForHelp(double value) {
    if (!std::isfinite(value)) {
        return "auto";
    }
    return formatOptionalDouble(value);
}

double defaultPseudocountForScoreMode(const std::string& scoreMode) {
    return scoreMode == "log_odds" ? 1.0 : 0.0;
}

std::string strandPartitionLabel(const std::string& strand) {
    return strand == "-" ? "minus" : "plus";
}

std::string denseScoreOutputExtension() {
#ifdef PSSM_SCAN_WITH_PARQUET
    return ".parquet";
#else
    return ".tsv";
#endif
}

std::string denseScoreFormatName() {
#ifdef PSSM_SCAN_WITH_PARQUET
    return "Parquet";
#else
    return "TSV";
#endif
}

std::filesystem::path denseScoreOutputPath(const std::string& outdir, const PSSM& pssm,
                                           const std::string& pssmFile, const std::string& scoreMode,
                                           double pseudocount, const std::string& chromosome,
                                           const std::string& strand, long from, long to, bool skipN) {
    std::filesystem::path outputPath = outdir;
    outputPath /= "tables";
    outputPath /= motifDatasetLabelFromPSSMFile(pssmFile);
    outputPath /= "motif_score_dense";
    outputPath /= "motif_id=" + pssm.motifID;
    outputPath /= "score_mode=" + scoreMode;
    outputPath /= "pseudocount=" + formatDoubleForFileLabel(pseudocount);
    outputPath /= "chrom=" + chromosome;
    outputPath /= "strand=" + strandPartitionLabel(strand);
    outputPath /= denseScorePartFilename(from, to, skipN, denseScoreOutputExtension());
    return outputPath;
}

#ifdef PSSM_SCAN_WITH_PARQUET
class DenseScoreParquetWriter {
public:
    explicit DenseScoreParquetWriter(const size_t blocksPerRecordBatch = 64)
        : blocksPerRecordBatch_(std::max<size_t>(1, blocksPerRecordBatch)) {}

    bool open(const std::filesystem::path& outputFilePath, std::string& error) {
        schema_ = arrow::schema({
            arrow::field("block_start", arrow::int64()),
            arrow::field("scores", arrow::list(arrow::float32()))
        });

        auto outputResult = arrow::io::FileOutputStream::Open(outputFilePath.string());
        if (!outputResult.ok()) {
            error = "opening Parquet output: " + outputResult.status().ToString();
            return false;
        }
        outputStream_ = *outputResult;

        parquet::WriterProperties::Builder writerPropertiesBuilder;
        writerPropertiesBuilder.compression(parquet::Compression::ZSTD);
        auto writerProperties = writerPropertiesBuilder.build();

        parquet::ArrowWriterProperties::Builder arrowPropertiesBuilder;
        arrowPropertiesBuilder.store_schema();
        auto arrowProperties = arrowPropertiesBuilder.build();

        auto writerResult = parquet::arrow::FileWriter::Open(
            *schema_, arrow::default_memory_pool(), outputStream_,
            writerProperties, arrowProperties);
        if (!writerResult.ok()) {
            error = "opening Parquet writer: " + writerResult.status().ToString();
            return false;
        }
        writer_ = std::move(writerResult).ValueOrDie();

        resetBuilders();
        return true;
    }

    bool writeBlock(const ScoreBlock& block, std::string& error) {
        if (!writer_) {
            error = "Parquet writer is not open.";
            return false;
        }

        auto status = blockStartBuilder_->Append(static_cast<std::int64_t>(block.blockStart));
        if (!status.ok()) {
            error = "appending dense block start: " + status.ToString();
            return false;
        }

        status = scoresBuilder_->Append();
        if (!status.ok()) {
            error = "appending dense score list: " + status.ToString();
            return false;
        }

        auto* valueBuilder = static_cast<arrow::FloatBuilder*>(scoresBuilder_->value_builder());
        for (double score : block.scores) {
            if (isSkippedScore(score)) {
                status = valueBuilder->AppendNull();
            } else {
                status = valueBuilder->Append(static_cast<float>(score));
            }
            if (!status.ok()) {
                error = "appending dense score: " + status.ToString();
                return false;
            }
        }

        pendingBlocks_++;
        if (pendingBlocks_ >= blocksPerRecordBatch_) {
            return flush(error);
        }
        return true;
    }

    bool close(std::string& error) {
        if (writer_) {
            if (!flush(error)) {
                return false;
            }
            auto status = writer_->Close();
            if (!status.ok()) {
                error = "closing Parquet writer: " + status.ToString();
                return false;
            }
            writer_.reset();
        }
        if (outputStream_) {
            auto status = outputStream_->Close();
            if (!status.ok()) {
                error = "closing Parquet output: " + status.ToString();
                return false;
            }
            outputStream_.reset();
        }
        return true;
    }

private:
    void resetBuilders() {
        auto* pool = arrow::default_memory_pool();
        blockStartBuilder_ = std::make_unique<arrow::Int64Builder>(pool);
        scoresBuilder_ = std::make_unique<arrow::ListBuilder>(
            pool, std::make_shared<arrow::FloatBuilder>(pool));
        pendingBlocks_ = 0;
    }

    bool flush(std::string& error) {
        if (pendingBlocks_ == 0) {
            return true;
        }

        std::shared_ptr<arrow::Array> blockStartArray;
        auto status = blockStartBuilder_->Finish(&blockStartArray);
        if (!status.ok()) {
            error = "finishing dense block_start array: " + status.ToString();
            return false;
        }

        std::shared_ptr<arrow::Array> scoresArray;
        status = scoresBuilder_->Finish(&scoresArray);
        if (!status.ok()) {
            error = "finishing dense scores array: " + status.ToString();
            return false;
        }

        auto batch = arrow::RecordBatch::Make(
            schema_, static_cast<std::int64_t>(pendingBlocks_),
            {blockStartArray, scoresArray});
        status = writer_->WriteRecordBatch(*batch);
        if (!status.ok()) {
            error = "writing dense Parquet record batch: " + status.ToString();
            return false;
        }

        resetBuilders();
        return true;
    }

    size_t blocksPerRecordBatch_;
    size_t pendingBlocks_ = 0;
    std::shared_ptr<arrow::Schema> schema_;
    std::shared_ptr<arrow::io::FileOutputStream> outputStream_;
    std::unique_ptr<parquet::arrow::FileWriter> writer_;
    std::unique_ptr<arrow::Int64Builder> blockStartBuilder_;
    std::unique_ptr<arrow::ListBuilder> scoresBuilder_;
};
#endif

void writeDenseScoreBlock(std::ostream& outFile, const ScoreBlock& block) {
    outFile << block.blockStart << "\t[";
    for (size_t i = 0; i < block.scores.size(); ++i) {
        if (i > 0) {
            outFile << ",";
        }
        if (isSkippedScore(block.scores[i])) {
            outFile << "NULL";
        } else {
            outFile << std::setprecision(9) << block.scores[i];
        }
    }
    outFile << "]\n";
}

std::string motifBaseID(const std::string& motifID) {
    const size_t versionSeparator = motifID.find('.');
    return versionSeparator == std::string::npos ? motifID : motifID.substr(0, versionSeparator);
}

std::string joinStrings(const std::vector<std::string>& values, const std::string& separator) {
    std::ostringstream joined;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            joined << separator;
        }
        joined << values[i];
    }
    return joined.str();
}

std::vector<std::string> findMotifVersionsInPSSMFile(const std::string& pssmFile, const std::string& targetMotifID) {
    std::vector<std::string> matches;
    const std::string targetBaseID = motifBaseID(targetMotifID);
    if (targetBaseID.empty()) {
        return matches;
    }

    std::ifstream inFile(pssmFile);
    if (!inFile.is_open()) {
        return matches;
    }

    std::string line;
    while (std::getline(inFile, line)) {
        line = PSSM::trim(line);
        if (line.empty() || line[0] != '>') {
            continue;
        }

        std::stringstream header(line.substr(1));
        std::string motifID;
        std::string motifName;
        header >> motifID >> motifName;
        if (motifID.empty() || motifID == targetMotifID || motifBaseID(motifID) != targetBaseID) {
            continue;
        }

        std::ostringstream label;
        label << motifID;
        if (!motifName.empty()) {
            label << " (" << motifName << ")";
        }
        matches.push_back(label.str());
    }

    return matches;
}

void maybePrintRequestedProgress(const char* operation, const PSSM& pssm,
                                 const std::string& chromosome, const std::string& strand,
                                 size_t posStart, size_t posEnd,
                                 size_t motifLength, size_t currentPos) {
    if (!consumeProgressStatusRequest()) {
        return;
    }

    size_t totalWindows = 0;
    if (posEnd >= motifLength && posStart <= posEnd - motifLength) {
        totalWindows = posEnd - motifLength - posStart + 1;
    }

    size_t completedWindows = 0;
    if (totalWindows > 0 && currentPos >= posStart) {
        completedWindows = std::min(totalWindows, currentPos - posStart + 1);
    }

    const double progress = totalWindows == 0 ? 1.0 :
        static_cast<double>(completedWindows) / static_cast<double>(totalWindows);

    std::ostringstream message;
    message << "\nI: progress request"
            << " operation=" << operation
            << " motif_id=" << pssm.motifID
            << " motif_name=" << pssm.motifName
            << " chr=" << chromosome
            << " strand=" << strand
            << " pos=" << currentPos
            << " windows=" << completedWindows << "/" << totalWindows
            << " progress=" << std::fixed << std::setprecision(3) << progress * 100.0 << "%";
    std::cerr << message.str() << std::endl;
}

void writeSparseHitHeader(std::ostream& outFile, bool showSequence) {
    outFile << "Chromosome\tFrom\tTo\tName\tScore\tStrand";
    if (showSequence) {
        outFile << "\tMatch";
    }
    outFile << "\tScoreMode\tPseudocount\tPWMRelativeScore\tCoordinateMode\n";
}

/** \brief Function to slide the PSSM across the DNA sequence and calculate scores
 * @param chromosome - the chromosome ID
 * @param sequence - the DNA sequence to scan
 * @param strand - the strand of the sequence
 * @param pssm - the PSSM to use for scoring
 * @param outFile - the output file stream to write results to
 * @param skipN - whether to skip windows containing 'N'
 * @param threshold - the minimal score to achieve to print
 * @param from - the minimal position on the chromosome to consider
 * @param to - the maximal position on the chromosome to consider
 * @param name - the name of the PSSM to be used in the output
 * @param showHeader - whether to show the header in the output file
 * @param scoreMode - score mode used for the Score column
 * @return 0 if successful, else -1
 */
int scanSequence(const std::string& chromosome, const EncodedChromosome& encodedChromosome, const std::string& strand,
                 const PSSM& pssm, const FlatPSSM& flatPssm,
                 std::ofstream& outFile, double threshold, double minPwmRelativeScore,
                 double maxPwmRelativeScore, long from, long to, bool showHeader,
                 bool showSequence, const std::string& scoreMode, double pseudocount, const ScoreRange& scoreRange,
                 CoordinateMode coordinateMode) {
    //size_t motifLength = pssm.begin()->second.size();
    const size_t motifLength = pssm.motifLength;
    const std::string& sequence = encodedChromosome.sequence;
    const size_t sequenceLength = sequence.size();
    const size_t reportInterval = std::max<size_t>(1, sequenceLength / 100 / 10);  // Update progress every 0.1%

    if (showDebug) std::cerr << "D: Sequence length=" << sequenceLength << ". report interval=" << reportInterval << std::endl;


    if (showHeader) {
        writeSparseHitHeader(outFile, showSequence);
    }

    size_t posStart=0L;
    if (from>0L) {
        posStart = (size_t) from;
        if (beVerbose) std::cerr << "I: Scanning chromosome " << chromosome << " from position " << from << " for motif " << pssm.motifName << std::endl;
    }

    size_t posEnd=sequenceLength;
    if (to > 0L && static_cast<size_t>(to) < sequenceLength) {
        posEnd = (size_t) to;
        if (beVerbose) std::cerr << "I: Scanning chromosome " << chromosome << " up to position " << to << " for motif " << pssm.motifName << std::endl;
    }

    // Slide the window over the sequence
    if (beVerbose) {
        std::cerr << "I: Scanning chromosome " << chromosome << " for motif " << pssm.motifName << " from " << posStart << " to " << posEnd << " - " << motifLength << std::endl;
    }

    if (beVerbose) displayProgressBar(0.0);
    if (posEnd < motifLength || posStart > posEnd - motifLength) {
        std::cerr << "W: Requested scan interval is shorter than motif " << pssm.motifName << "." << std::endl;
        return 0;
    }

    const bool reverseComplementWindow = strand == "-";
    const std::vector<std::uint8_t>& codes = reverseComplementWindow ? encodedChromosome.minusCodes : encodedChromosome.plusCodes;
    if (reverseComplementWindow && codes.empty()) {
        std::cerr << "E: Minus-strand scan requested, but chromosome " << chromosome << " lacks reverse-complement codes." << std::endl;
        return -1;
    }
    const size_t statusCheckMask = 0x3fff;
    for (size_t i = posStart; i <= posEnd - motifLength; ++i) {
        const double score = calculateScoreAtGenomicStart(codes, sequenceLength, reverseComplementWindow, i, flatPssm);
        const double relativeScore = pwmRelativeScore(score, scoreRange);
        // Skip output if the window contained 'N' or invalid nucleotides
        // std::cerr << "D: Score: " << score << std::endl;

        // Progress indicator
        if (beVerbose && i % reportInterval == 0) {
            const size_t denominator = posEnd - motifLength - posStart;
            const double progress = denominator == 0 ? 1.0 : static_cast<double>(i - posStart) / denominator;
            //auto now = std::chrono::high_resolution_clock::now();
            //auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
            //std::cout << "Progress: " << std::fixed << std::setprecision(2) << progress * 100 << "% - Elapsed time: " << elapsed << " seconds\n";
            displayProgressBar(progress);
        }
        if (((i - posStart) & statusCheckMask) == 0) {
            maybePrintRequestedProgress("scan", pssm, chromosome, strand, posStart, posEnd, motifLength, i);
        }

        if (score < threshold) {
            continue;
        }
        if (score < SENTINEL_SCORE / 10.0) {
            continue;
        }
        if (std::isfinite(minPwmRelativeScore) && (!std::isfinite(relativeScore) || relativeScore < minPwmRelativeScore)) {
            continue;
        }
        if (std::isfinite(maxPwmRelativeScore) && (!std::isfinite(relativeScore) || relativeScore > maxPwmRelativeScore)) {
            continue;
        }

        const size_t outputStart = outputStartForCoordinateMode(coordinateMode, i);
        const size_t outputEnd = outputEndForCoordinateMode(coordinateMode, i, motifLength);
        outFile << chromosome << "\t" << outputStart << "\t" << outputEnd << "\t" << pssm.motifName
                              << "\t" << std::fixed << std::setprecision(3) << score << "\t" << strand;
        if (showSequence) {
            outFile << "\t" << sequenceWindowForOutput(sequence, i, motifLength, reverseComplementWindow);
        }
        outFile << "\t" << scoreMode
                << "\t" << std::fixed << std::setprecision(6) << pseudocount
                << "\t" << std::fixed << std::setprecision(6) << relativeScore
                << "\t" << coordinateModeName(coordinateMode);
        outFile << '\n';

    }
    maybePrintRequestedProgress("scan", pssm, chromosome, strand, posStart, posEnd, motifLength, posEnd - motifLength);
    return 0;
}

/** \brief Slide the PSSM across one sequence and collect score distribution bins.
 */
int scanScoreDistribution(const std::string& chromosome, const EncodedChromosome& encodedChromosome, const std::string& strand,
                          const PSSM& pssm, const FlatPSSM& flatPssm, long from, long to,
                          ScoreDistribution& distribution) {
    const size_t motifLength = pssm.motifLength;
    const std::string& sequence = encodedChromosome.sequence;
    const size_t sequenceLength = sequence.size();
    const size_t reportInterval = std::max<size_t>(1, sequenceLength / 100 / 10);

    if (showDebug) std::cerr << "D: Distribution sequence length=" << sequenceLength << ". report interval=" << reportInterval << std::endl;

    size_t posStart = 0L;
    if (from > 0L) {
        posStart = static_cast<size_t>(from);
        if (beVerbose) std::cerr << "I: Collecting score distribution on chromosome " << chromosome << " from position " << from << " for motif " << pssm.motifName << std::endl;
    }

    size_t posEnd = sequenceLength;
    if (to > 0L && to < static_cast<long>(sequenceLength)) {
        posEnd = static_cast<size_t>(to);
        if (beVerbose) std::cerr << "I: Collecting score distribution on chromosome " << chromosome << " up to position " << to << " for motif " << pssm.motifName << std::endl;
    }

    if (posEnd < motifLength || posStart > posEnd - motifLength) {
        std::cerr << "W: Requested distribution interval is shorter than motif " << pssm.motifName << "." << std::endl;
        return 0;
    }

    if (beVerbose) displayProgressBar(0.0);
    const bool reverseComplementWindow = strand == "-";
    const std::vector<std::uint8_t>& codes = reverseComplementWindow ? encodedChromosome.minusCodes : encodedChromosome.plusCodes;
    if (reverseComplementWindow && codes.empty()) {
        std::cerr << "E: Minus-strand distribution requested, but chromosome " << chromosome << " lacks reverse-complement codes." << std::endl;
        return -1;
    }
    const size_t statusCheckMask = 0x3fff;
    const size_t lastWindowStart = posEnd - motifLength;
    for (size_t blockStart = posStart; blockStart <= lastWindowStart; blockStart += DEFAULT_SCORE_BLOCK_SIZE) {
        const size_t windowCount = std::min(DEFAULT_SCORE_BLOCK_SIZE, lastWindowStart - blockStart + 1);
        const ScoreBlock block = calculateScoreBlock(codes, sequenceLength, reverseComplementWindow, blockStart, windowCount, flatPssm);
        for (double score : block.scores) {
            distribution.add(score);
        }

        if (beVerbose && blockStart % reportInterval == 0) {
            const size_t denominator = posEnd - motifLength - posStart;
            const double progress = denominator == 0 ? 1.0 : static_cast<double>(blockStart - posStart) / denominator;
            displayProgressBar(progress);
        }
        if (((blockStart - posStart) & statusCheckMask) == 0) {
            maybePrintRequestedProgress("score_distribution", pssm, chromosome, strand, posStart, posEnd, motifLength, blockStart);
        }
    }

    maybePrintRequestedProgress("score_distribution", pssm, chromosome, strand, posStart, posEnd, motifLength, posEnd - motifLength);
    return 0;
}

/** \brief Slide the PSSM across one sequence and write dense score blocks.
 */
int scanDenseScores(const std::string& chromosome, const EncodedChromosome& encodedChromosome, const std::string& strand,
                    const PSSM& pssm, const FlatPSSM& flatPssm, long from, long to,
                    size_t denseBlockSize, const std::filesystem::path& outputFilePath) {
    if (denseBlockSize == 0) {
        std::cerr << "E: --dense-block-size must be greater than 0." << std::endl;
        return -1;
    }

    const size_t motifLength = pssm.motifLength;
    const std::string& sequence = encodedChromosome.sequence;
    const size_t sequenceLength = sequence.size();
    const size_t reportInterval = std::max<size_t>(1, sequenceLength / 100 / 10);

    size_t posStart = 0L;
    if (from > 0L) {
        posStart = static_cast<size_t>(from);
        if (beVerbose) std::cerr << "I: Writing dense scores on chromosome " << chromosome << " from position " << from << " for motif " << pssm.motifName << std::endl;
    }

    size_t posEnd = sequenceLength;
    if (to > 0L && to < static_cast<long>(sequenceLength)) {
        posEnd = static_cast<size_t>(to);
        if (beVerbose) std::cerr << "I: Writing dense scores on chromosome " << chromosome << " up to position " << to << " for motif " << pssm.motifName << std::endl;
    }

    if (posEnd < motifLength || posStart > posEnd - motifLength) {
        std::cerr << "W: Requested dense-score interval is shorter than motif " << pssm.motifName << "." << std::endl;
        return 0;
    }

    const bool reverseComplementWindow = strand == "-";
    const std::vector<std::uint8_t>& codes = reverseComplementWindow ? encodedChromosome.minusCodes : encodedChromosome.plusCodes;
    if (reverseComplementWindow && codes.empty()) {
        std::cerr << "E: Minus-strand dense-score scan requested, but chromosome " << chromosome << " lacks reverse-complement codes." << std::endl;
        return -1;
    }

    std::filesystem::create_directories(outputFilePath.parent_path());
#ifdef PSSM_SCAN_WITH_PARQUET
    DenseScoreParquetWriter denseWriter;
    std::string denseWriterError;
    if (!denseWriter.open(outputFilePath, denseWriterError)) {
        std::cerr << "E: Error opening dense score output file '" << outputFilePath.string()
                  << "': " << denseWriterError << std::endl;
        return -1;
    }
#else
    std::ofstream outFile(outputFilePath);
    if (!outFile.is_open()) {
        std::cerr << "E: Error opening dense score output file '" << outputFilePath.string() << "'." << std::endl;
        return -1;
    }
    outFile << "block_start\tscores\n";
#endif

    if (beVerbose) displayProgressBar(0.0);
    const size_t statusCheckMask = 0x3fff;
    const size_t lastWindowStart = posEnd - motifLength;
    std::uint64_t validWindows = 0;
    std::uint64_t skippedWindows = 0;
    for (size_t blockStart = posStart; blockStart <= lastWindowStart; blockStart += denseBlockSize) {
        const size_t windowCount = std::min(denseBlockSize, lastWindowStart - blockStart + 1);
        const ScoreBlock block = calculateScoreBlock(codes, sequenceLength, reverseComplementWindow, blockStart, windowCount, flatPssm);
#ifdef PSSM_SCAN_WITH_PARQUET
        if (!denseWriter.writeBlock(block, denseWriterError)) {
            std::cerr << "E: Error writing dense score output file '" << outputFilePath.string()
                      << "': " << denseWriterError << std::endl;
            return -1;
        }
#else
        writeDenseScoreBlock(outFile, block);
#endif
        validWindows += block.validWindows;
        skippedWindows += block.skippedWindows;

        if (beVerbose && blockStart % reportInterval == 0) {
            const size_t denominator = posEnd - motifLength - posStart;
            const double progress = denominator == 0 ? 1.0 : static_cast<double>(blockStart - posStart) / denominator;
            displayProgressBar(progress);
        }
        if (((blockStart - posStart) & statusCheckMask) == 0) {
            maybePrintRequestedProgress("dense_scores", pssm, chromosome, strand, posStart, posEnd, motifLength, blockStart);
        }
    }

#ifdef PSSM_SCAN_WITH_PARQUET
    if (!denseWriter.close(denseWriterError)) {
        std::cerr << "E: Error closing dense score output file '" << outputFilePath.string()
                  << "': " << denseWriterError << std::endl;
        return -1;
    }
#endif
    maybePrintRequestedProgress("dense_scores", pssm, chromosome, strand, posStart, posEnd, motifLength, posEnd - motifLength);
    if (beVerbose) {
        std::cerr << "I: Dense " << denseScoreFormatName() << " score blocks saved to " << outputFilePath.string()
                  << " (valid_windows=" << validWindows
                  << ", skipped_windows=" << skippedWindows << ")." << std::endl;
    }
    return 0;
}

void writeScoreDistribution(std::ofstream& outFile, const ScoreDistribution& distribution, const PSSM& pssm,
                            const std::string& chromosome, const std::string& strand, const std::string& scoreMode,
                            double pseudocount) {
    outFile << "MotifID\tMotifName\tChromosome\tStrand\tScoreMode\tPseudocount\tBinScheme\tBinWidth\tValidWindows\tSkippedWindows\tMinScore\tMaxScore\tMeanScore\tScoreBinStart\tScoreBinEnd\tBinCount\n";
    if (distribution.skippedWindows > 0) {
        outFile << pssm.motifID << "\t" << pssm.motifName << "\t" << chromosome << "\t" << strand << "\t" << scoreMode
                << "\t" << std::fixed << std::setprecision(6) << pseudocount
                << "\tsentinel"
                << "\tNA"
                << "\t" << distribution.validWindows
                << "\t" << distribution.skippedWindows
                << "\t" << formatOptionalDouble(distribution.minScore)
                << "\t" << formatOptionalDouble(distribution.maxScore)
                << "\t" << formatOptionalDouble(distribution.meanScore())
                << "\t-Inf"
                << "\t-10000"
                << "\t" << distribution.skippedWindows
                << '\n';
    }
    for (const auto& [bin, count] : distribution.bins) {
        outFile << pssm.motifID << "\t" << pssm.motifName << "\t" << chromosome << "\t" << strand << "\t" << scoreMode
                << "\t" << std::fixed << std::setprecision(6) << pseudocount
                << "\t" << distribution.binScheme
                << "\t" << std::fixed << std::setprecision(6) << bin.width
                << "\t" << distribution.validWindows
                << "\t" << distribution.skippedWindows
                << "\t" << formatOptionalDouble(distribution.minScore)
                << "\t" << formatOptionalDouble(distribution.maxScore)
                << "\t" << formatOptionalDouble(distribution.meanScore())
                << "\t" << bin.start
                << "\t" << bin.end
                << "\t" << count
                << '\n';
    }
}

typedef std::vector<FastaIndexEntry> fasta_index_type;
typedef std::vector<GzipIndexEntry> gzip_index_type;

std::uint16_t readLittleEndianUint16(const unsigned char* bytes) {
    return static_cast<std::uint16_t>(bytes[0]) |
           (static_cast<std::uint16_t>(bytes[1]) << 8);
}

std::uint64_t readLittleEndianUint64(const unsigned char* bytes) {
    std::uint64_t value = 0;
    for (int i = 7; i >= 0; --i) {
        value = (value << 8) | bytes[i];
    }
    return value;
}

int readFastaIndexFile(const std::string& indexFile, fasta_index_type& fastaIndex) {
    std::ifstream inFile(indexFile);
    if (!inFile.is_open()) {
        return -1;
    }

    std::set<std::string> chromosomeNames;
    std::string line;
    while (std::getline(inFile, line)) {
        line = PSSM::trim(line);
        if (line.empty()) continue;

        std::stringstream ss(line);
        FastaIndexEntry entry;
        if (!(ss >> entry.name >> entry.length >> entry.offset >> entry.basesPerLine >> entry.bytesPerLine)) {
            std::cerr << "E: Could not parse FASTA index line: " << line << std::endl;
            return -1;
        }
        if (entry.basesPerLine == 0 || entry.bytesPerLine < entry.basesPerLine) {
            std::cerr << "E: Invalid FASTA index line geometry: " << line << std::endl;
            return -1;
        }
        if (!chromosomeNames.insert(entry.name).second) {
            std::cerr << "E: Duplicate chromosome " << entry.name
                      << " in FASTA index " << indexFile << std::endl;
            return -1;
        }
        fastaIndex.push_back(std::move(entry));
    }

    return fastaIndex.empty() ? -1 : 0;
}

int selectFastaIndexEntries(const fasta_index_type& fastaIndex,
                            const chromosome_set_type& targetChromosomes,
                            std::vector<FastaIndexEntry>& selectedEntries) {
    chromosome_set_type missingChromosomes = targetChromosomes;
    selectedEntries.clear();
    selectedEntries.reserve(targetChromosomes.empty() ? fastaIndex.size() : targetChromosomes.size());

    for (const FastaIndexEntry& entry : fastaIndex) {
        if (targetChromosomes.empty() || targetChromosomes.count(entry.name) != 0) {
            selectedEntries.push_back(entry);
            missingChromosomes.erase(entry.name);
        }
    }

    if (!missingChromosomes.empty()) {
        std::vector<std::string> names(missingChromosomes.begin(), missingChromosomes.end());
        std::cerr << "E: Chromosome(s) not found in FASTA index: "
                  << joinStrings(names, ", ") << std::endl;
        return -1;
    }
    return selectedEntries.empty() ? -1 : 0;
}

int readGzipIndexFile(const std::string& indexFile, gzip_index_type& gzipIndex) {
    std::ifstream inFile(indexFile, std::ios::binary);
    if (!inFile.is_open()) {
        return -1;
    }

    std::array<unsigned char, 8> countBytes{};
    inFile.read(reinterpret_cast<char*>(countBytes.data()), countBytes.size());
    if (inFile.gcount() != static_cast<std::streamsize>(countBytes.size())) {
        std::cerr << "E: Could not read BGZF .gzi entry count from " << indexFile << std::endl;
        return -1;
    }

    const std::uint64_t entryCount = readLittleEndianUint64(countBytes.data());
    gzipIndex.reserve(static_cast<size_t>(entryCount) + 1);
    gzipIndex.push_back(GzipIndexEntry{0, 0});

    for (std::uint64_t i = 0; i < entryCount; ++i) {
        std::array<unsigned char, 16> pairBytes{};
        inFile.read(reinterpret_cast<char*>(pairBytes.data()), pairBytes.size());
        if (inFile.gcount() != static_cast<std::streamsize>(pairBytes.size())) {
            std::cerr << "E: Could not read BGZF .gzi entry " << i << " from " << indexFile << std::endl;
            return -1;
        }
        gzipIndex.push_back(GzipIndexEntry{
            readLittleEndianUint64(pairBytes.data()),
            readLittleEndianUint64(pairBytes.data() + 8)
        });
    }

    std::sort(gzipIndex.begin(), gzipIndex.end(), [](const GzipIndexEntry& a, const GzipIndexEntry& b) {
        return a.uncompressedOffset < b.uncompressedOffset;
    });

    if (beVerbose) std::cerr << "I: Read " << gzipIndex.size() << " BGZF index anchors from " << indexFile << std::endl;
    return 0;
}

GzipIndexEntry gzipAnchorForOffset(const gzip_index_type& gzipIndex, std::uint64_t uncompressedOffset) {
    GzipIndexEntry selected{0, 0};
    for (const auto& entry : gzipIndex) {
        if (entry.uncompressedOffset > uncompressedOffset) {
            break;
        }
        selected = entry;
    }
    return selected;
}

class BgzfBlockReader {
public:
    BgzfBlockReader(const std::string& filename, const gzip_index_type& gzipIndex, std::uint64_t startOffset) {
        file_.open(filename, std::ios::binary);
        if (!file_.is_open()) {
            throw std::runtime_error("Failed to open BGZF file: " + filename);
        }

        const GzipIndexEntry anchor = gzipAnchorForOffset(gzipIndex, startOffset);
        file_.seekg(static_cast<std::streamoff>(anchor.compressedOffset), std::ios::beg);
        if (!file_) {
            throw std::runtime_error("Failed to seek BGZF file: " + filename);
        }
        currentUncompressedOffset_ = anchor.uncompressedOffset;
        if (startOffset > currentUncompressedOffset_ && !skip(startOffset - currentUncompressedOffset_)) {
            throw std::runtime_error("Failed to seek within BGZF stream: " + filename);
        }
    }

    bool read(char* destination, size_t bytesToRead) {
        size_t copied = 0;
        while (copied < bytesToRead) {
            if (blockPosition_ >= block_.size()) {
                if (!loadNextBlock()) {
                    return false;
                }
                if (block_.empty()) {
                    continue;
                }
            }
            const size_t available = block_.size() - blockPosition_;
            const size_t take = std::min(available, bytesToRead - copied);
            std::memcpy(destination + copied, block_.data() + blockPosition_, take);
            blockPosition_ += take;
            copied += take;
            currentUncompressedOffset_ += take;
        }
        return true;
    }

    bool skip(std::uint64_t bytesToSkip) {
        std::array<char, 8192> buffer{};
        while (bytesToSkip > 0) {
            const size_t take = static_cast<size_t>(std::min<std::uint64_t>(bytesToSkip, buffer.size()));
            if (!read(buffer.data(), take)) {
                return false;
            }
            bytesToSkip -= take;
        }
        return true;
    }

private:
    bool loadNextBlock() {
        block_.clear();
        blockPosition_ = 0;

        std::array<unsigned char, 12> header{};
        file_.read(reinterpret_cast<char*>(header.data()), header.size());
        if (file_.gcount() == 0) {
            return false;
        }
        if (file_.gcount() != static_cast<std::streamsize>(header.size())) {
            throw std::runtime_error("Truncated BGZF header");
        }
        if (header[0] != 31 || header[1] != 139 || header[2] != 8 || !(header[3] & 4)) {
            throw std::runtime_error("Compressed FASTA is not BGZF gzip with an extra field");
        }

        const std::uint16_t extraLength = readLittleEndianUint16(header.data() + 10);
        std::vector<unsigned char> extra(extraLength);
        file_.read(reinterpret_cast<char*>(extra.data()), extra.size());
        if (file_.gcount() != static_cast<std::streamsize>(extra.size())) {
            throw std::runtime_error("Truncated BGZF extra field");
        }

        bool foundBlockSize = false;
        std::uint16_t blockSizeMinusOne = 0;
        size_t extraPosition = 0;
        while (extraPosition + 4 <= extra.size()) {
            const unsigned char subfield1 = extra[extraPosition];
            const unsigned char subfield2 = extra[extraPosition + 1];
            const std::uint16_t subfieldLength = readLittleEndianUint16(extra.data() + extraPosition + 2);
            extraPosition += 4;
            if (extraPosition + subfieldLength > extra.size()) {
                throw std::runtime_error("Invalid BGZF extra subfield length");
            }
            if (subfield1 == 'B' && subfield2 == 'C' && subfieldLength == 2) {
                blockSizeMinusOne = readLittleEndianUint16(extra.data() + extraPosition);
                foundBlockSize = true;
                break;
            }
            extraPosition += subfieldLength;
        }
        if (!foundBlockSize) {
            throw std::runtime_error("BGZF block size subfield not found");
        }

        const size_t headerSize = header.size() + extra.size();
        const size_t totalBlockSize = static_cast<size_t>(blockSizeMinusOne) + 1;
        if (totalBlockSize < headerSize) {
            throw std::runtime_error("Invalid BGZF block size");
        }

        std::vector<unsigned char> compressedBlock(totalBlockSize);
        std::memcpy(compressedBlock.data(), header.data(), header.size());
        std::memcpy(compressedBlock.data() + header.size(), extra.data(), extra.size());
        const size_t remainingCompressedBytes = totalBlockSize - headerSize;
        file_.read(reinterpret_cast<char*>(compressedBlock.data() + headerSize), static_cast<std::streamsize>(remainingCompressedBytes));
        if (file_.gcount() != static_cast<std::streamsize>(remainingCompressedBytes)) {
            throw std::runtime_error("Truncated BGZF block");
        }

        z_stream stream{};
        if (inflateInit2(&stream, 31) != Z_OK) {
            throw std::runtime_error("inflateInit2 failed for BGZF block");
        }

        std::array<unsigned char, 65536> output{};
        stream.next_in = compressedBlock.data();
        stream.avail_in = static_cast<uInt>(compressedBlock.size());
        stream.next_out = output.data();
        stream.avail_out = static_cast<uInt>(output.size());

        const int inflateResult = inflate(&stream, Z_FINISH);
        if (inflateResult != Z_STREAM_END) {
            inflateEnd(&stream);
            throw std::runtime_error("Failed to inflate BGZF block");
        }
        const size_t outputBytes = stream.total_out;
        inflateEnd(&stream);

        block_.assign(reinterpret_cast<char*>(output.data()), reinterpret_cast<char*>(output.data()) + outputBytes);
        return true;
    }

    std::ifstream file_;
    std::vector<char> block_;
    size_t blockPosition_ = 0;
    std::uint64_t currentUncompressedOffset_ = 0;
};

int readIndexedPlainFastaChromosome(const std::string& fastaFile, const FastaIndexEntry& entry,
                                    std::string& sequence) {
    std::ifstream inFile(fastaFile, std::ios::binary);
    if (!inFile.is_open()) {
        std::cerr << "E: Error opening FASTA file: " << fastaFile << std::endl;
        return -1;
    }

    inFile.seekg(static_cast<std::streamoff>(entry.offset), std::ios::beg);
    if (!inFile) {
        std::cerr << "E: Could not seek to chromosome " << entry.name << " in FASTA file " << fastaFile << std::endl;
        return -1;
    }

    sequence.clear();
    sequence.reserve(static_cast<size_t>(entry.length));
    std::uint64_t remainingBases = entry.length;
    const std::uint64_t lineBreakBytes = entry.bytesPerLine - entry.basesPerLine;

    while (remainingBases > 0) {
        const size_t basesToRead = static_cast<size_t>(std::min<std::uint64_t>(entry.basesPerLine, remainingBases));
        std::string buffer(basesToRead, '\0');
        inFile.read(buffer.data(), static_cast<std::streamsize>(basesToRead));
        if (inFile.gcount() != static_cast<std::streamsize>(basesToRead)) {
            std::cerr << "E: Truncated FASTA sequence while reading chromosome " << entry.name << std::endl;
            return -1;
        }
        sequence += buffer;
        remainingBases -= basesToRead;
        if (remainingBases > 0 && lineBreakBytes > 0) {
            inFile.seekg(static_cast<std::streamoff>(lineBreakBytes), std::ios::cur);
            if (!inFile) {
                std::cerr << "E: Could not skip FASTA line ending while reading chromosome " << entry.name << std::endl;
                return -1;
            }
        }
    }

    return 0;
}

int readIndexedBgzfFastaChromosome(const std::string& fastaFile, const FastaIndexEntry& entry,
                                   const gzip_index_type& gzipIndex, std::string& sequence) {
    try {
        BgzfBlockReader reader(fastaFile, gzipIndex, entry.offset);
        sequence.clear();
        sequence.reserve(static_cast<size_t>(entry.length));
        std::uint64_t remainingBases = entry.length;
        const std::uint64_t lineBreakBytes = entry.bytesPerLine - entry.basesPerLine;

        while (remainingBases > 0) {
            const size_t basesToRead = static_cast<size_t>(std::min<std::uint64_t>(entry.basesPerLine, remainingBases));
            std::string buffer(basesToRead, '\0');
            if (!reader.read(buffer.data(), basesToRead)) {
                std::cerr << "E: Truncated BGZF FASTA sequence while reading chromosome " << entry.name << std::endl;
                return -1;
            }
            sequence += buffer;
            remainingBases -= basesToRead;
            if (remainingBases > 0 && lineBreakBytes > 0 && !reader.skip(lineBreakBytes)) {
                std::cerr << "E: Could not skip BGZF FASTA line ending while reading chromosome " << entry.name << std::endl;
                return -1;
            }
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "W: Could not use BGZF .gzi random access for " << fastaFile << ": " << e.what() << std::endl;
        return -1;
    }
}

int readIndexedFastaChromosome(const std::string& fastaFile,
                               CompressedFileReader::FileType fileType,
                               const FastaIndexEntry& entry,
                               const gzip_index_type& gzipIndex,
                               std::string& sequence) {
    if (fileType == CompressedFileReader::PLAIN) {
        return readIndexedPlainFastaChromosome(fastaFile, entry, sequence);
    }
    return readIndexedBgzfFastaChromosome(fastaFile, entry, gzipIndex, sequence);
}

std::string sanitizeMotifNameForFile(std::string motifName) {
    for (const char character : std::string(":/()")) {
        std::replace(motifName.begin(), motifName.end(), character, '-');
    }
    return motifName;
}

struct PreparedMotif {
    std::string motifID;
    std::unique_ptr<PSSM> pssm;
    FlatPSSM flatPssm;
    ScoreRange scoreRange;
    std::string motifNameForFile;
    std::unique_ptr<ScoreDistribution> distribution;
};

void printHelp(const std::string& programName, const std::string& genomeFile, const std::string& fastaIndexFile, const std::string& gzipIndexFile, const std::string& pssmFile, const std::string& targetMotifID, double threshold, double minPwmRelativeScore, double maxPwmRelativeScore, double pseudocount, const std::string& chromosome, long from, long to, const std::string& regions, const std::string& outdir, bool showSequence, const std::string& scoreMode, bool scoreDistribution, const std::string& distributionBinWidth, bool denseScores, size_t denseBlockSize, StrandMode strandMode, CoordinateMode coordinateMode) {
    std::cout << "Usage: " << programName << " [-v] [-c chromosome] [-t toBp] [-f fromBp] [-g genome_file] [-p pssm_file] [-m motif_id] [--score-mode log2_relative_risk|log_odds] [--pseudocount value] [--strand +|-|both] [--coordinate-mode legacy|bed] [--score-distribution] [--distribution-bin-width adaptive|width] [--dense-scores] [--dense-block-size windows] [--skip-N | --neutral-N] [--skip-normalization]" << std::endl;
    std::cout << " -v, --verbose        Allow verbose output (set to " << beVerbose << ")" << std::endl;
    std::cout << " -d, --debug          Allow debug output (set to " << showDebug << ")" << std::endl;
    std::cout << " -g, --genome         Path to genome FASTA file (set to '" << genomeFile  << "')" << std::endl;
    std::cout << " --fasta-index        Required FASTA .fai index (default '" << (fastaIndexFile.empty() ? genomeFile + ".fai" : fastaIndexFile) << "')" << std::endl;
    std::cout << " --gzip-index         Required BGZF .gzi index for FASTA.gz (default '" << (gzipIndexFile.empty() ? genomeFile + ".gzi" : gzipIndexFile) << "')" << std::endl;
    std::cout << " -p, --pssm           Path to JASPAR PSSM file (set to '" << pssmFile << "')" << std::endl;
    std::cout << " -m, --motif          Target motif ID from JASPAR file (set to '" << targetMotifID << "')" << std::endl;
    std::cout << " -l, --threshold      Minimal score to achieve to print (set to " << threshold << ")" << std::endl;
    std::cout << " --min-pwm-relative-score Minimal PWM-relative score to print (set to " << formatScoreBoundForHelp(minPwmRelativeScore) << ")" << std::endl;
    std::cout << " --max-pwm-relative-score Maximal PWM-relative score to print (set to " << formatScoreBoundForHelp(maxPwmRelativeScore) << ")" << std::endl;
    std::cout << " -c, --chr            Single chromosome to consider (set to '" << chromosome << "')" << std::endl;
    std::cout << " -f, --from           Minimal position on chromosome to consider (set to " << from << ")" << std::endl;
    std::cout << " -t, --to             Maximal position on chromosome to consider (set to " << to << ")" << std::endl;
    std::cout << " -r, --regions        Regions constraint file (set to '" << regions << "')" << std::endl;
    std::cout << " -o, --outdir         Directory to create output files in (set to '" << outdir << "')" << std::endl;
    std::cout << " -s, --show-sequence  Adds the sequence matched by the motif as another field in the output (set to " << showSequence << ")" << std::endl;
    std::cout << " --score-mode         PSSM score mode: log2_relative_risk or log_odds (set to '" << scoreMode << "')" << std::endl;
    std::cout << " --pseudocount        Count added to each A/C/G/T motif entry before normalization (set to " << formatPseudocountForHelp(pseudocount) << "; auto is 0 for log2_relative_risk and 1 for log_odds)" << std::endl;
    std::cout << " --strand             Strand to scan: +, -, or both (set to '" << strandModeName(strandMode) << "')" << std::endl;
    std::cout << " --coordinate-mode    Output coordinate convention: legacy or bed (set to '" << coordinateModeName(coordinateMode) << "')" << std::endl;
    std::cout << " --score-distribution Write a score histogram instead of BED hits; default strand is + unless --strand is set (set to " << scoreDistribution << ")" << std::endl;
    std::cout << " --distribution-bin-width Score histogram bin width or adaptive ladder (set to '" << distributionBinWidth << "')" << std::endl;
    std::cout << " --dense-scores       Write dense " << denseScoreFormatName() << " score blocks for one motif and one chromosome (set to " << denseScores << ")" << std::endl;
    std::cout << " --dense-block-size   Dense alignment scores per output block (set to " << denseBlockSize << ")" << std::endl;
    std::cout << " --skip-N             Skip windows containing 'N'" << std::endl;
    std::cout << " --neutral-N          Treat 'N' as neutral (contribute 0 to the score)" << std::endl;
    std::cout << " -N, --skip-normalization Skip log-normalisation, will affect scoring." << std::endl;
#ifdef _WIN32
    std::cout << " Ctrl+Break            Request one progress line on stderr in a Windows console" << std::endl;
#else
    std::cout << " SIGUSR1              Request one progress line on stderr from a running process, e.g. kill -USR1 <pid>" << std::endl;
#ifdef SIGINFO
    std::cout << " SIGINFO              Also requests one progress line where available, e.g. Ctrl-T on BSD/macOS" << std::endl;
#endif
#endif
    std::cout << " -h, --help           Display this help message" << std::endl;
}

} // namespace

int main(int argc, char* argv[]) {
    installProgressRequestHandlers();

    // Default file names
    std::string genomeFile = "Homo_sapiens.GRCh38.dna.primary_assembly.fasta";
    std::string fastaIndexFile;
    std::string gzipIndexFile;
    std::string pssmFile = "JASPAR2026_CORE_non-redundant_pfms_jaspar.txt";
    std::string regionsFile = "";
    std::string targetMotifID;
    double threshold = SENTINEL_SCORE; // Very small number, i.e. effectively no threshold
    bool thresholdSet = false;
    double minPwmRelativeScore = -std::numeric_limits<double>::infinity();
    double maxPwmRelativeScore = std::numeric_limits<double>::infinity();
    double pseudocount = std::numeric_limits<double>::quiet_NaN();
    bool showSequence = false;
    bool skipN = true;  // Default to skipping N
    bool skipNormalization = false;  // Default to not skip normalization
    bool showHelp = false; // Do not show help after all variables have been assigned
    long targetFrom = -1L, targetTo = -1L ;
    std::string targetChromosome ;
    std::vector<Region> regions;
    std::string outdir = ".";
    std::string scoreMode = "log2_relative_risk";
    bool scoreDistribution = false;
    std::string distributionBinWidth = "adaptive";
    bool denseScores = false;
    size_t denseBlockSize = DEFAULT_SCORE_BLOCK_SIZE;
    StrandMode strandMode = StrandMode::Both;
    bool strandModeSet = false;
    CoordinateMode coordinateMode = CoordinateMode::Legacy;

    // Option flags and variables for getopt
    int option;
    static struct option long_options[] = {
        {"genome", required_argument, 0, 'g'},
        {"fasta-index", required_argument, 0, 0},
        {"gzip-index", required_argument, 0, 0},
        {"pssm", required_argument, 0, 'p'},
        {"motif", required_argument, 0, 'm'},
        {"threshold", required_argument, 0, 'l'},
        {"min-pwm-relative-score", required_argument, 0, 0},
        {"max-pwm-relative-score", required_argument, 0, 0},
        {"pwm-relative-threshold", required_argument, 0, 0},
        {"chr", required_argument, 0, 'c'},
        {"from", required_argument, 0, 'f'},
        {"to", required_argument, 0, 't'},
        {"regions", required_argument, 0, 'r'},  // regions file
        {"outdir", required_argument, 0, 'o'},  // output directory
        {"show-sequence", no_argument, 0, 's'},  // output directory
        {"score-mode", required_argument, 0, 0},
        {"pseudocount", required_argument, 0, 0},
        {"strand", required_argument, 0, 0},
        {"coordinate-mode", required_argument, 0, 0},
        {"score-distribution", no_argument, 0, 0},
        {"distribution-bin-width", required_argument, 0, 0},
        {"dense-scores", no_argument, 0, 0},
        {"dense-block-size", required_argument, 0, 0},
        {"skip-N", no_argument, 0, 0},
        {"neutral-N", no_argument, 0, 0},
        {"skip-normalization", no_argument, 0, 'N'},
        {"verbose", no_argument, 0, 'v'},
        {"debug", no_argument, 0, 'd'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    // Parse command line arguments using getopt
    int option_index = 0;
    while ((option = getopt_long(argc, argv, "g:p:m:l:c:f:t:o:Nsvh", long_options, &option_index)) != -1) {
        switch (option) {
            case 'g':
                genomeFile = optarg;
                break;
            case 'p':
                pssmFile = optarg;
                break;
            case 'm':
                targetMotifID = optarg;
                if (beVerbose) std::cerr << "I: Setting motif ID to " << targetMotifID << std::endl;
                break;
            case 'r':
                regionsFile = optarg;
                break;
            case 'c':
                targetChromosome = optarg;
                if (beVerbose) std::cerr << "I: Only showing matches on chromosome " << targetChromosome << "." << std::endl;
                break;
            case 'f':
                if (!parseLongStrict(optarg, targetFrom)) {
                    std::cerr << "E: --from expects an integer, got '" << optarg << "'." << std::endl;
                    showHelp = 1;
                    break;
                }
                if (beVerbose) std::cerr << "I: Only showing matches with position downstream of " << targetFrom << "." << std::endl;
                break;
            case 't':
                if (!parseLongStrict(optarg, targetTo)) {
                    std::cerr << "E: --to expects an integer, got '" << optarg << "'." << std::endl;
                    showHelp = 1;
                    break;
                }
                if (beVerbose) std::cerr << "I: Only showing matches with position up to " << targetTo << "." << std::endl;
                break;
            case 'l':
                if (!parseDoubleStrict(optarg, threshold)) {
                    std::cerr << "E: --threshold expects a finite numeric score, got '" << optarg << "'." << std::endl;
                    showHelp = 1;
                    break;
                }
                thresholdSet = true;
                if (beVerbose) std::cerr << "I: Only showing matches with score > " << threshold << "." << std::endl;
                break;
            case 'o':
                outdir = optarg;
                if (!std::filesystem::is_directory(outdir)) {
                    if (std::filesystem::create_directory(outdir)) {
                        std::cerr << "I: Created output directory '" << outdir << "'." << std::endl;
                    } else {
                        std::cerr << "E: Output directory '" << outdir << "' not existing and could not be created." << std::endl;
                    }
                }
                break;
            case 0:
                if (std::string(long_options[option_index].name) == "skip-N") {
                    skipN = true;
                } else if (std::string(long_options[option_index].name) == "neutral-N") {
                    skipN = false;
                } else if (std::string(long_options[option_index].name) == "fasta-index") {
                    fastaIndexFile = optarg;
                } else if (std::string(long_options[option_index].name) == "gzip-index") {
                    gzipIndexFile = optarg;
                } else if (std::string(long_options[option_index].name) == "score-mode") {
                    const std::string canonicalScoreMode = PSSM::canonicalScoreModeName(optarg);
                    if (canonicalScoreMode.empty()) {
                        std::cerr << "E: Unsupported score mode '" << optarg << "'. Use log2_relative_risk or log_odds." << std::endl;
                        showHelp = 1;
                    } else {
                        scoreMode = canonicalScoreMode;
                    }
                } else if (std::string(long_options[option_index].name) == "pseudocount") {
                    if (!parseDoubleStrict(optarg, pseudocount) || pseudocount < 0.0) {
                        std::cerr << "E: --pseudocount expects a finite non-negative numeric value, got '" << optarg << "'." << std::endl;
                        showHelp = 1;
                    }
                } else if (std::string(long_options[option_index].name) == "strand") {
                    if (!parseStrandMode(optarg, strandMode)) {
                        std::cerr << "E: Unsupported strand mode '" << optarg << "'. Use +, -, or both." << std::endl;
                        showHelp = 1;
                    } else {
                        strandModeSet = true;
                    }
                } else if (std::string(long_options[option_index].name) == "coordinate-mode") {
                    if (!parseCoordinateMode(optarg, coordinateMode)) {
                        std::cerr << "E: Unsupported coordinate mode '" << optarg << "'. Use legacy or bed." << std::endl;
                        showHelp = 1;
                    }
                } else if (std::string(long_options[option_index].name) == "min-pwm-relative-score" ||
                           std::string(long_options[option_index].name) == "pwm-relative-threshold") {
                    if (!parseDoubleStrict(optarg, minPwmRelativeScore)) {
                        std::cerr << "E: --min-pwm-relative-score expects a finite numeric score, got '" << optarg << "'." << std::endl;
                        showHelp = 1;
                        break;
                    }
                    if (minPwmRelativeScore < 0.0 || minPwmRelativeScore > 1.0) {
                        std::cerr << "E: --min-pwm-relative-score must be between 0 and 1." << std::endl;
                        showHelp = 1;
                    }
                } else if (std::string(long_options[option_index].name) == "max-pwm-relative-score") {
                    if (!parseDoubleStrict(optarg, maxPwmRelativeScore)) {
                        std::cerr << "E: --max-pwm-relative-score expects a finite numeric score, got '" << optarg << "'." << std::endl;
                        showHelp = 1;
                        break;
                    }
                    if (maxPwmRelativeScore < 0.0 || maxPwmRelativeScore > 1.0) {
                        std::cerr << "E: --max-pwm-relative-score must be between 0 and 1." << std::endl;
                        showHelp = 1;
                    }
                } else if (std::string(long_options[option_index].name) == "score-distribution") {
                    scoreDistribution = true;
                } else if (std::string(long_options[option_index].name) == "distribution-bin-width") {
                    distributionBinWidth = optarg;
                    std::string normalizedBinWidth = distributionBinWidth;
                    std::transform(normalizedBinWidth.begin(), normalizedBinWidth.end(), normalizedBinWidth.begin(), [](unsigned char c) {
                        return static_cast<char>(std::tolower(c));
                    });
                    double parsedBinWidth = 0.0;
                    if (normalizedBinWidth != "adaptive" && (!parseDoubleStrict(distributionBinWidth.c_str(), parsedBinWidth) || parsedBinWidth <= 0)) {
                        std::cerr << "E: --distribution-bin-width must be 'adaptive' or greater than 0." << std::endl;
                        showHelp = 1;
                    }
                } else if (std::string(long_options[option_index].name) == "dense-scores") {
                    denseScores = true;
                } else if (std::string(long_options[option_index].name) == "dense-block-size") {
                    long parsedDenseBlockSize = 0;
                    if (!parseLongStrict(optarg, parsedDenseBlockSize) || parsedDenseBlockSize <= 0L) {
                        std::cerr << "E: --dense-block-size expects a positive integer, got '" << optarg << "'." << std::endl;
                        showHelp = 1;
                    } else {
                        denseBlockSize = static_cast<size_t>(parsedDenseBlockSize);
                    }
                }
                break;
            case 's':
                showSequence=true;
                break;
            case 'N':
                skipNormalization=true;
                break;
            case 'v':
                beVerbose = 1;
                break;
            case 'd':
                beVerbose = 2;
		showDebug = 1;
                break;
            case 'h':
                printHelp(argv[0],genomeFile,fastaIndexFile,gzipIndexFile,pssmFile,targetMotifID,threshold,minPwmRelativeScore,maxPwmRelativeScore,pseudocount,targetChromosome,targetFrom,targetTo,regionsFile,outdir,showSequence,scoreMode,scoreDistribution,distributionBinWidth,denseScores,denseBlockSize,strandMode,coordinateMode);
                return 0;
            default:
                showHelp = 1;
                break;
        }

        if (showHelp) {
            printHelp(argv[0],genomeFile,fastaIndexFile,gzipIndexFile,pssmFile,targetMotifID,threshold,minPwmRelativeScore,maxPwmRelativeScore,pseudocount,targetChromosome,targetFrom,targetTo,regionsFile,outdir,showSequence,scoreMode,scoreDistribution,distributionBinWidth,denseScores,denseBlockSize,strandMode,coordinateMode);
            return 1;
        }
    }

    if (minPwmRelativeScore > maxPwmRelativeScore) {
        std::cerr << "E: --min-pwm-relative-score must not be greater than --max-pwm-relative-score." << std::endl;
        printHelp(argv[0],genomeFile,fastaIndexFile,gzipIndexFile,pssmFile,targetMotifID,threshold,minPwmRelativeScore,maxPwmRelativeScore,pseudocount,targetChromosome,targetFrom,targetTo,regionsFile,outdir,showSequence,scoreMode,scoreDistribution,distributionBinWidth,denseScores,denseBlockSize,strandMode,coordinateMode);
        return 1;
    }
    if (scoreDistribution && denseScores) {
        std::cerr << "E: --score-distribution and --dense-scores are separate output modes; choose one." << std::endl;
        return 1;
    }

/*
    // Ensure the motif ID is provided
    if (targetMotifID.empty()) {
        std::cerr << "E: You must specify a motif ID with the -m option." << std::endl;
        printHelp(argv[0],genomeFile,pssmFile,targetMotifID,threshold,targetChromosome,targetFrom,targetTo,regionsFile,outdir,showSequence,scoreMode,scoreDistribution,distributionBinWidth);
        return 1;
    }
*/
    const std::string effectiveScoreMode = skipNormalization ? "raw_counts" : scoreMode;
    const double effectivePseudocount = skipNormalization ? 0.0 :
        (std::isfinite(pseudocount) ? pseudocount : defaultPseudocountForScoreMode(scoreMode));
    if (skipNormalization && scoreMode != "log2_relative_risk") {
        std::cerr << "W: --skip-normalization ignores --score-mode '" << scoreMode << "'; output ScoreMode will be raw_counts." << std::endl;
    }
    if (skipNormalization && std::isfinite(pseudocount)) {
        std::cerr << "W: --skip-normalization ignores --pseudocount " << pseudocount << "." << std::endl;
    }
    PSSM::debug = showDebug;

    if (scoreDistribution && !strandModeSet) {
        strandMode = StrandMode::Plus;
    }
    if (denseScores) {
        if (targetMotifID.empty()) {
            std::cerr << "E: --dense-scores requires a single --motif to avoid accidental all-motif dense output." << std::endl;
            return 1;
        }
        if (targetChromosome.empty()) {
            std::cerr << "E: --dense-scores requires a single --chr to avoid accidental whole-genome dense output." << std::endl;
            return 1;
        }
        if (thresholdSet) {
            std::cerr << "W: --dense-scores writes all alignment scores and ignores --threshold." << std::endl;
        }
        if (std::isfinite(minPwmRelativeScore) || std::isfinite(maxPwmRelativeScore)) {
            std::cerr << "W: --dense-scores writes all alignment scores and ignores PWM-relative score filters." << std::endl;
        }
        if (showSequence) {
            std::cerr << "W: --dense-scores ignores --show-sequence." << std::endl;
        }
    }
    const bool scanMinusStrandRequested = strandMode == StrandMode::Minus || strandMode == StrandMode::Both;

    if (!regionsFile.empty()) {
        regions = Region::parseRegionsFile(regionsFile);
        if (regions.empty()) {
            std::cerr << "E: No valid regions found in regions file." << std::endl;
            return 1;
        }
    }

    if (!regions.empty() && (!targetChromosome.empty() || targetFrom > 0L || targetTo > 0L)) {
        std::cerr << "E: Cannot specify both individual chromosome options and a regions file." << std::endl;
        return 1;
    }
    if (denseScores && !regions.empty()) {
        std::cerr << "E: --dense-scores does not support --regions; use --chr with optional --from/--to." << std::endl;
        return 1;
    }

    std::filesystem::path sparseHitOutputDirectory;
    if (!scoreDistribution && !denseScores) {
        const HitOutputOptions hitOutputOptions{
            effectiveScoreMode,
            effectivePseudocount,
            thresholdSet,
            threshold,
            minPwmRelativeScore,
            maxPwmRelativeScore,
            coordinateMode,
            showSequence,
            skipN
        };
        sparseHitOutputDirectory = hitOutputDirectory(outdir, hitOutputOptions);
        std::error_code directoryError;
        std::filesystem::create_directories(sparseHitOutputDirectory, directoryError);
        if (directoryError) {
            std::cerr << "E: Could not create hit output directory '"
                      << sparseHitOutputDirectory.string() << "': "
                      << directoryError.message() << std::endl;
            return 1;
        }
    }

    chromosome_set_type targetChromosomesForGenome;
    if (!targetChromosome.empty()) {
        targetChromosomesForGenome.insert(targetChromosome);
    } else {
        for (const auto& region : regions) {
            targetChromosomesForGenome.insert(region.chromosome);
        }
    }

    // Load the PSSM matrix from a JASPAR-like file
    pssm_list_type pssm_list;

    if ( PSSM::parsePSSMFile(pssmFile, pssm_list , targetMotifID , beVerbose + showDebug ) ) {
        std::cerr << "E: Error parsing PSSM file '" << pssmFile << "'" << std::endl;
        return 1;
    }
    if (showDebug) std::cerr << "D: Read " << pssm_list.size() << " PSSMs from file '" << pssmFile << "'" << std::endl;
    if (pssm_list.empty()) {
        if (!targetMotifID.empty()) {
            std::cerr << "E: Target motif " << targetMotifID << " was not found in PSSM file '"
                      << pssmFile << "'." << std::endl;
            const std::vector<std::string> availableVersions = findMotifVersionsInPSSMFile(pssmFile, targetMotifID);
            if (!availableVersions.empty()) {
                std::cerr << "I: Available motif version(s) with base ID "
                          << motifBaseID(targetMotifID) << ": "
                          << joinStrings(availableVersions, ", ") << std::endl;
            }
        } else {
            std::cerr << "E: PSSM file '" << pssmFile
                      << "' parsed successfully but no motifs were loaded." << std::endl;
        }
        return 1;
    }

    std::vector<std::string> motifIDs;
    motifIDs.reserve(pssm_list.size());
    for (const auto& motifEntry : pssm_list) {
        motifIDs.push_back(motifEntry.first);
    }
    std::sort(motifIDs.begin(), motifIDs.end());

    std::vector<PreparedMotif> preparedMotifs;
    preparedMotifs.reserve(motifIDs.size());
    for (const std::string& motifID : motifIDs) {
        const PSSM& parsedPssm = pssm_list.at(motifID);
        auto normalizedPssm = std::make_unique<PSSM>(parsedPssm);
        if (normalizedPssm->pssm.empty()) {
            std::cerr << "E: Failed to parse PSSM " << motifID << "." << std::endl;
            return 1;
        }
        if (motifID != normalizedPssm->motifID) {
            std::cerr << "E: Mismatch in motif ID for " << motifID << "." << std::endl;
            return 1;
        }

        const size_t inferredMotifLength = normalizedPssm->pssm.begin()->second.size();
        if (inferredMotifLength != static_cast<size_t>(normalizedPssm->motifLength)) {
            std::cerr << "W: Motif " << motifID << " length mismatch: first matrix row has "
                      << inferredMotifLength << " bp, parsed motif length is "
                      << normalizedPssm->motifLength << " bp." << std::endl;
        } else if (beVerbose) {
            std::cerr << "I: Motif " << motifID << " (" << normalizedPssm->motifName
                      << ") length: " << normalizedPssm->motifLength << " bp" << std::endl;
        }

        if (skipNormalization) {
            if (beVerbose) std::cerr << "I: Skipping normalization for " << motifID << "." << std::endl;
        } else {
            normalizedPssm->normalizePSSM(backgroundFrequencies, scoreMode, effectivePseudocount);
        }

        const ScoreRange scoreRange = scoreRangeForPSSM(*normalizedPssm);
        FlatPSSM flatPssm = flattenPSSM(*normalizedPssm, skipN);
        if (showDebug) {
            std::cerr << "PSSM after normalization:" << std::endl << *normalizedPssm;
        }

        preparedMotifs.push_back(PreparedMotif{
            motifID,
            std::move(normalizedPssm),
            std::move(flatPssm),
            scoreRange,
            sanitizeMotifNameForFile(parsedPssm.motifName),
            scoreDistribution ? std::make_unique<ScoreDistribution>(distributionBinWidth) : nullptr
        });
    }
    pssm_list.clear();

    const std::string effectiveFastaIndexFile = fastaIndexFile.empty()
        ? genomeFile + ".fai" : fastaIndexFile;
    fasta_index_type fastaIndex;
    if (readFastaIndexFile(effectiveFastaIndexFile, fastaIndex) != 0) {
        std::cerr << "E: Indexed scanning requires a readable FASTA index: "
                  << effectiveFastaIndexFile << std::endl;
        return 1;
    }

    std::vector<FastaIndexEntry> selectedFastaEntries;
    if (selectFastaIndexEntries(fastaIndex, targetChromosomesForGenome,
                                selectedFastaEntries) != 0) {
        std::cerr << "E: No indexed chromosomes selected for scanning." << std::endl;
        return 1;
    }

    const CompressedFileReader::FileType genomeFileType =
        CompressedFileReader::determineFileType(genomeFile);
    if (genomeFileType == CompressedFileReader::BZIP2) {
        std::cerr << "E: Indexed production scans do not support bzip2 FASTA; "
                  << "use plain FASTA or BGZF with .fai and .gzi indexes." << std::endl;
        return 1;
    }

    gzip_index_type gzipIndex;
    if (genomeFileType == CompressedFileReader::GZIP) {
        const std::string effectiveGzipIndexFile = gzipIndexFile.empty()
            ? genomeFile + ".gzi" : gzipIndexFile;
        if (readGzipIndexFile(effectiveGzipIndexFile, gzipIndex) != 0) {
            std::cerr << "E: Indexed gzip scans require BGZF and a readable .gzi index: "
                      << effectiveGzipIndexFile << std::endl;
            return 1;
        }
    }

    for (size_t chromosomeIndex = 0; chromosomeIndex < selectedFastaEntries.size();
         ++chromosomeIndex) {
        const FastaIndexEntry& fastaEntry = selectedFastaEntries[chromosomeIndex];
        const bool firstChromosome = chromosomeIndex == 0;
        const bool lastChromosome = chromosomeIndex + 1 == selectedFastaEntries.size();

        std::string chromosomeSequence;
        if (readIndexedFastaChromosome(genomeFile, genomeFileType, fastaEntry,
                                       gzipIndex, chromosomeSequence) != 0) {
            std::cerr << "E: Could not read indexed chromosome " << fastaEntry.name
                      << " from " << genomeFile << std::endl;
            return 1;
        }
        EncodedChromosome encodedChromosome =
            encodeChromosome(std::move(chromosomeSequence), scanMinusStrandRequested);
        const std::string& chromosome = fastaEntry.name;
        if (beVerbose) {
            std::cerr << "I: Loaded indexed chromosome " << chromosome
                      << " (" << encodedChromosome.sequence.size()
                      << " bp); scanning all selected motifs before release." << std::endl;
        }

        // Iterate over every motif while this is the only chromosome in memory.
        for (PreparedMotif& preparedMotif : preparedMotifs) {
            const std::string& motifID = preparedMotif.motifID;
            const PSSM& preparedPssm = *preparedMotif.pssm;
            const FlatPSSM& flatPssm = preparedMotif.flatPssm;
            const ScoreRange& scoreRange = preparedMotif.scoreRange;
            const std::string& motifNameForFile = preparedMotif.motifNameForFile;

            if (scoreDistribution) {
                const std::string chromosomeLabel = targetChromosome.empty() ? "all" : targetChromosome;
                const std::string distributionStrandLabel = strandModeName(strandMode);
                const std::string distributionStrandFileLabel = strandModeFileLabel(strandMode);
                std::string distributionBinWidthForFile = distributionBinWidth;
                std::replace(distributionBinWidthForFile.begin(), distributionBinWidthForFile.end(), '/', '-');
                std::replace(distributionBinWidthForFile.begin(), distributionBinWidthForFile.end(), '\\', '-');
                std::replace(distributionBinWidthForFile.begin(), distributionBinWidthForFile.end(), ':', '-');
                std::string outputFileNameDistribution = motifNameForFile + "_" + motifID + "_score_distribution_" + effectiveScoreMode + "_bins_" + distributionBinWidthForFile + "_" + distributionStrandFileLabel + "_" + chromosomeLabel;
                if (effectivePseudocount > 0.0) {
                    outputFileNameDistribution += "_pseudocount_" + formatDoubleForFileLabel(effectivePseudocount);
                }
                if (targetFrom > 0L) {
                    outputFileNameDistribution += "_" + std::to_string(targetFrom);
                }
                if (targetTo > 0L) {
                    outputFileNameDistribution += "-" + std::to_string(targetTo);
                }
                outputFileNameDistribution += ".tsv";

                std::filesystem::path outputFilePathDistribution = outdir;
                outputFilePathDistribution /= outputFileNameDistribution;
                outputFileNameDistribution = outputFilePathDistribution.string();

                ScoreDistribution& distribution = *preparedMotif.distribution;

                if (!regions.empty()) {
                    for (const auto& region : regions) {
                        if (region.chromosome != chromosome) {
                            continue;
                        }
                        const std::string& sequence = encodedChromosome.sequence;
                        const long seqLength = static_cast<long>(sequence.length());
                        const long regionFrom = std::max(0L, region.from);
                        const long regionTo = std::min(seqLength, region.to);
                        if (strandMode == StrandMode::Plus || strandMode == StrandMode::Both) {
                            if (scanScoreDistribution(chromosome, encodedChromosome, "+", preparedPssm,
                                                      flatPssm, regionFrom, regionTo, distribution) != 0) {
                                return 1;
                            }
                        }
                        if (strandMode == StrandMode::Minus || strandMode == StrandMode::Both) {
                            if (scanScoreDistribution(chromosome, encodedChromosome, "-", preparedPssm,
                                                      flatPssm, regionFrom, regionTo, distribution) != 0) {
                                return 1;
                            }
                        }
                    }
                } else {
                    if (strandMode == StrandMode::Plus || strandMode == StrandMode::Both) {
                        if (scanScoreDistribution(chromosome, encodedChromosome, "+", preparedPssm,
                                                  flatPssm, targetFrom, targetTo, distribution) != 0) {
                            return 1;
                        }
                    }
                    if (strandMode == StrandMode::Minus || strandMode == StrandMode::Both) {
                        if (scanScoreDistribution(chromosome, encodedChromosome, "-", preparedPssm,
                                                  flatPssm, targetFrom, targetTo, distribution) != 0) {
                            return 1;
                        }
                    }
                }

                if (lastChromosome) {
                    std::ofstream outFileDistribution(outputFileNameDistribution);
                    if (!outFileDistribution.is_open()) {
                        std::cerr << "E: Error opening output file '" << outputFileNameDistribution << "'." << std::endl;
                        return 1;
                    }
                    writeScoreDistribution(outFileDistribution, distribution, preparedPssm,
                                           chromosomeLabel, distributionStrandLabel,
                                           effectiveScoreMode, effectivePseudocount);
                    std::cout << "Score distribution saved to: " << outputFileNameDistribution << std::endl;
                }
                continue;
            }

            if (denseScores) {
                const bool scanPlusStrand = strandMode == StrandMode::Plus || strandMode == StrandMode::Both;
                const bool scanMinusStrand = strandMode == StrandMode::Minus || strandMode == StrandMode::Both;

                if (scanPlusStrand) {
                    const std::filesystem::path outputFilePath = denseScoreOutputPath(
                        outdir, preparedPssm, pssmFile, effectiveScoreMode, effectivePseudocount,
                        targetChromosome, "+", targetFrom, targetTo, skipN);
                    if (scanDenseScores(targetChromosome, encodedChromosome, "+", preparedPssm, flatPssm,
                                        targetFrom, targetTo, denseBlockSize, outputFilePath) != 0) {
                        return 1;
                    }
                    std::cout << "Dense " << denseScoreFormatName() << " scores saved to: " << outputFilePath.string() << std::endl;
                }

                if (scanMinusStrand) {
                    const std::filesystem::path outputFilePath = denseScoreOutputPath(
                        outdir, preparedPssm, pssmFile, effectiveScoreMode, effectivePseudocount,
                        targetChromosome, "-", targetFrom, targetTo, skipN);
                    if (scanDenseScores(targetChromosome, encodedChromosome, "-", preparedPssm, flatPssm,
                                        targetFrom, targetTo, denseBlockSize, outputFilePath) != 0) {
                        return 1;
                    }
                    std::cout << "Dense " << denseScoreFormatName() << " scores saved to: " << outputFilePath.string() << std::endl;
                }

                continue;
            }

            std::filesystem::path outputFilePathPositive = sparseHitOutputDirectory;
            outputFilePathPositive /= motifNameForFile + "_" + motifID + "_positive";
            std::string outputFileNamePositive = outputFilePathPositive.string();

            std::filesystem::path outputFilePathNegative = sparseHitOutputDirectory;
            outputFilePathNegative /= motifNameForFile + "_" + motifID + "_negative";
            std::string outputFileNameNegative = outputFilePathNegative.string();

            if (!targetChromosome.empty()) {
                outputFileNamePositive += "_"+targetChromosome;
                outputFileNameNegative += "_"+targetChromosome;
            }
            if (targetFrom > 0L) {
                outputFileNamePositive += "_"+std::to_string(targetFrom);
                outputFileNameNegative += "_"+std::to_string(targetFrom);
            }
            if (targetTo > 0L) {
                outputFileNamePositive += "-"+std::to_string(targetTo);
                outputFileNameNegative += "-"+std::to_string(targetTo);
            }
            outputFileNamePositive += ".bed";
            outputFileNameNegative += ".bed";

            const bool scanPlusStrand = strandMode == StrandMode::Plus || strandMode == StrandMode::Both;
            const bool scanMinusStrand = strandMode == StrandMode::Minus || strandMode == StrandMode::Both;

            std::ofstream outFilePositive;
            std::ofstream outFileNegative;
            const std::ios::openmode outputMode = firstChromosome
                ? (std::ios::out | std::ios::trunc)
                : (std::ios::out | std::ios::app);
            if (scanPlusStrand) {
                outFilePositive.open(outputFileNamePositive, outputMode);
            }
            if (scanMinusStrand) {
                outFileNegative.open(outputFileNameNegative, outputMode);
            }

            if (scanPlusStrand && !outFilePositive.is_open()) {
                std::cerr << "E: Error opening output file '" << outputFileNamePositive << "'." << std::endl;
                return 1;
            }

            if (scanMinusStrand && !outFileNegative.is_open()) {
                std::cerr << "E: Error opening output file '" << outputFileNameNegative << "'." << std::endl;
                return 1;
            }

            bool showHeader = firstChromosome;

            if (!regions.empty()) {
                // Scan specified regions
                for (const auto& region : regions) {
                    if (region.chromosome != chromosome) {
                        continue;
                    }

                    const std::string& sequence = encodedChromosome.sequence;
                    const long seqLength = static_cast<long>(sequence.length());
                    const long regionFrom = std::max(0L, region.from);
                    const long regionTo = std::min(seqLength, region.to);

                    if (regionFrom >= regionTo) {
                        std::cerr << "E: Invalid region: " << region.chromosome << " " << region.from << "-" << region.to << std::endl;
                        continue;
                    }

                    if (beVerbose && scanPlusStrand) std::cerr << "I: Scanning region: " << region.chromosome << ":" << regionFrom << "-" << regionTo << " (" << region.name << ") (+ strand)" << std::endl;

                    // Scan positive strand
                    if (scanPlusStrand) {
                        if (scanSequence(chromosome, encodedChromosome, "+", preparedPssm,
                                         flatPssm, outFilePositive, threshold, minPwmRelativeScore,
                                         maxPwmRelativeScore, regionFrom, regionTo, showHeader,
                                         showSequence, effectiveScoreMode, effectivePseudocount,
                                         scoreRange, coordinateMode) != 0) {
                            return 1;
                        }
                    }

                    // Scan negative strand
                    if (beVerbose && scanMinusStrand) std::cerr << "I: Scanning region: " << region.chromosome << ":" << regionFrom << "-" << regionTo << " (" << region.name << ") (- strand)" << std::endl;
                    if (scanMinusStrand) {
                        if (scanSequence(chromosome, encodedChromosome, "-", preparedPssm,
                                         flatPssm, outFileNegative, threshold, minPwmRelativeScore,
                                         maxPwmRelativeScore, regionFrom, regionTo, showHeader,
                                         showSequence, effectiveScoreMode, effectivePseudocount,
                                         scoreRange, coordinateMode) != 0) {
                            return 1;
                        }
                    }

                    showHeader = false;
                }

            } else {
                if (beVerbose && scanPlusStrand) std::cerr << "I: Scanning chromosome " << chromosome << " from " << targetFrom << " to " << targetTo << " (+ strand)" << std::endl;

                // Scan positive strand
                if (scanPlusStrand) {
                    if (scanSequence(chromosome, encodedChromosome, "+", preparedPssm,
                                     flatPssm, outFilePositive, threshold, minPwmRelativeScore,
                                     maxPwmRelativeScore, targetFrom, targetTo, showHeader,
                                     showSequence, effectiveScoreMode, effectivePseudocount,
                                     scoreRange, coordinateMode) != 0) {
                        return 1;
                    }
                }

                // Scan negative strand (reverse complement)
                if (beVerbose && scanMinusStrand) std::cerr << "I: Scanning chromosome " << chromosome << " from " << targetFrom << " to " << targetTo << " (- strand)" << std::endl;
                if (scanMinusStrand) {
                    if (scanSequence(chromosome, encodedChromosome, "-", preparedPssm,
                                     flatPssm, outFileNegative, threshold, minPwmRelativeScore,
                                     maxPwmRelativeScore, targetFrom, targetTo, showHeader,
                                     showSequence, effectiveScoreMode, effectivePseudocount,
                                     scoreRange, coordinateMode) != 0) {
                        return 1;
                    }
                }
            }

            if (scanPlusStrand) outFilePositive.close();
            if (scanMinusStrand) outFileNegative.close();

            if (!lastChromosome) {
                continue;
            }
            if (scanPlusStrand && scanMinusStrand) {
                std::cout << "Results saved to: " << outputFileNamePositive << " and " << outputFileNameNegative << std::endl;
            } else if (scanPlusStrand) {
                std::cout << "Results saved to: " << outputFileNamePositive << std::endl;
            } else {
                std::cout << "Results saved to: " << outputFileNameNegative << std::endl;
            }
        } // prepared motif
    } // indexed chromosome

    return 0;
}
