#include "pssm_scan_core.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

std::array<std::uint8_t, 256> buildBaseCodeTable() {
    std::array<std::uint8_t, 256> table{};
    table.fill(BASE_N);
    table[static_cast<unsigned char>('A')] = BASE_A;
    table[static_cast<unsigned char>('a')] = BASE_A;
    table[static_cast<unsigned char>('C')] = BASE_C;
    table[static_cast<unsigned char>('c')] = BASE_C;
    table[static_cast<unsigned char>('G')] = BASE_G;
    table[static_cast<unsigned char>('g')] = BASE_G;
    table[static_cast<unsigned char>('T')] = BASE_T;
    table[static_cast<unsigned char>('t')] = BASE_T;
    table[static_cast<unsigned char>('N')] = BASE_N;
    table[static_cast<unsigned char>('n')] = BASE_N;
    return table;
}

const std::array<std::uint8_t, 256>& baseCodeTable() {
    static const std::array<std::uint8_t, 256> table = buildBaseCodeTable();
    return table;
}

std::string sanitizePathComponent(const std::string& value) {
    std::string sanitized;
    sanitized.reserve(value.size());
    bool previousWasDash = false;
    for (const unsigned char c : value) {
        if (std::isalnum(c) || c == '.' || c == '_' || c == '-') {
            sanitized.push_back(static_cast<char>(c));
            previousWasDash = false;
        } else if (!previousWasDash) {
            sanitized.push_back('-');
            previousWasDash = true;
        }
    }
    while (!sanitized.empty() && sanitized.back() == '-') {
        sanitized.pop_back();
    }
    return sanitized.empty() ? "unknown" : sanitized;
}

std::string optionalScoreBoundLabel(const double value) {
    return std::isfinite(value) ? formatDoubleForFileLabel(value) : "none";
}

std::string coordinateModeFileLabel(const CoordinateMode coordinateMode) {
    return coordinateMode == CoordinateMode::Bed ? "bed" : "legacy";
}

} // namespace

char complementBase(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N';
    }
}

bool parseDoubleStrict(const char* text, double& value) {
    errno = 0;
    char* end = nullptr;
    const double parsed = std::strtod(text, &end);
    if (errno != 0 || end == text || (end != nullptr && *end != '\0') || !std::isfinite(parsed)) {
        return false;
    }
    value = parsed;
    return true;
}

std::string formatDoubleForFileLabel(const double value) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value < 0.0 ? "minf" : "inf";
    }

    std::ostringstream ss;
    ss << std::setprecision(12) << value;
    std::string label = ss.str();
    std::replace(label.begin(), label.end(), '+', 'p');
    std::replace(label.begin(), label.end(), '-', 'm');
    return label;
}

std::filesystem::path hitOutputDirectory(const std::filesystem::path& outdir,
                                         const HitOutputOptions& options) {
    std::filesystem::path outputPath = outdir;
    outputPath /= "hits";
    outputPath /= "score_mode-" + sanitizePathComponent(options.scoreMode);
    outputPath /= "pseudocount-" + formatDoubleForFileLabel(options.pseudocount);
    outputPath /= "threshold-" +
        (options.thresholdSet ? formatDoubleForFileLabel(options.threshold) : "none");
    outputPath /= "pwm_relative_min-" + optionalScoreBoundLabel(options.minPwmRelativeScore);
    outputPath /= "pwm_relative_max-" + optionalScoreBoundLabel(options.maxPwmRelativeScore);
    outputPath /= "coordinate_mode-" + coordinateModeFileLabel(options.coordinateMode);
    outputPath /= std::string("sequence-") + (options.showSequence ? "included" : "omitted");
    outputPath /= std::string("n_policy-") + (options.skipN ? "skip" : "neutral");
    return outputPath;
}

std::string motifDatasetLabelFromPSSMFile(const std::filesystem::path& pssmFile) {
    std::string filename = pssmFile.filename().string();
    std::transform(filename.begin(), filename.end(), filename.begin(), [](const unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });

    const size_t jasparPosition = filename.find("jaspar");
    if (jasparPosition != std::string::npos) {
        size_t versionStart = jasparPosition + std::string("jaspar").size();
        while (versionStart < filename.size() &&
               !std::isdigit(static_cast<unsigned char>(filename[versionStart]))) {
            ++versionStart;
        }
        size_t versionEnd = versionStart;
        while (versionEnd < filename.size() &&
               std::isdigit(static_cast<unsigned char>(filename[versionEnd]))) {
            ++versionEnd;
        }
        if (versionEnd > versionStart) {
            return "jaspar" + filename.substr(versionStart, versionEnd - versionStart);
        }
    }

    std::string stem = pssmFile.stem().string();
    std::transform(stem.begin(), stem.end(), stem.begin(), [](const unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return sanitizePathComponent(stem);
}

std::string denseScorePartFilename(const long from, const long to, const bool skipN,
                                   const std::string& extension) {
    const std::string fromLabel = from > 0L ? std::to_string(from) : "0";
    const std::string toLabel = to > 0L ? std::to_string(to) : "end";
    return "part-from=" + fromLabel + "-to=" + toLabel +
        "-n_policy=" + (skipN ? "skip" : "neutral") + "-000000" + extension;
}

std::uint8_t codeForBase(char base) {
    return baseCodeTable()[static_cast<unsigned char>(base)];
}

std::uint8_t complementCode(std::uint8_t code) {
    switch (code) {
        case BASE_A: return BASE_T;
        case BASE_C: return BASE_G;
        case BASE_G: return BASE_C;
        case BASE_T: return BASE_A;
        default: return BASE_N;
    }
}

bool isSkippedScore(double score) {
    return !std::isfinite(score) || score < SENTINEL_SCORE / 10.0;
}

double calculateScoreAt(const std::vector<std::uint8_t>& codes, size_t start, const FlatPSSM& pssm) {
    double score = 0.0;

    for (size_t i = 0; i < pssm.motifLength; ++i) {
        score += pssm.scores[i * BASE_CODE_COUNT + codes[start + i]];
    }

    return score;
}

double calculateScoreAtGenomicStart(const std::vector<std::uint8_t>& codes, size_t sequenceLength,
                                    bool reverseComplementWindow, size_t genomicStart,
                                    const FlatPSSM& pssm) {
    if (pssm.motifLength == 0 || genomicStart > sequenceLength ||
        pssm.motifLength > sequenceLength - genomicStart) {
        return SENTINEL_SCORE;
    }

    const size_t codeStart = reverseComplementWindow ?
        sequenceLength - genomicStart - pssm.motifLength :
        genomicStart;
    if (codeStart > codes.size() || pssm.motifLength > codes.size() - codeStart) {
        return SENTINEL_SCORE;
    }
    return calculateScoreAt(codes, codeStart, pssm);
}

ScoreBlock calculateScoreBlock(const std::vector<std::uint8_t>& codes, size_t sequenceLength,
                               bool reverseComplementWindow, size_t blockStart,
                               size_t windowCount, const FlatPSSM& pssm) {
    ScoreBlock block;
    block.blockStart = blockStart;
    block.scores.reserve(windowCount);

    for (size_t offset = 0; offset < windowCount; ++offset) {
        const double score = calculateScoreAtGenomicStart(
            codes, sequenceLength, reverseComplementWindow, blockStart + offset, pssm);
        block.scores.push_back(score);
        if (isSkippedScore(score)) {
            block.skippedWindows++;
        } else {
            block.validWindows++;
        }
    }

    return block;
}

double calculateScore(const std::string& window, const pssm_type& pssm, bool skipN) {
    FlatPSSM flatPssm;
    flatPssm.motifLength = window.size();
    flatPssm.scores.assign(flatPssm.motifLength * BASE_CODE_COUNT, skipN ? SENTINEL_SCORE : 0.0);
    for (size_t i = 0; i < flatPssm.motifLength; ++i) {
        if (pssm.find('A') != pssm.end()) flatPssm.scores[i * BASE_CODE_COUNT + BASE_A] = pssm.at('A')[i];
        if (pssm.find('C') != pssm.end()) flatPssm.scores[i * BASE_CODE_COUNT + BASE_C] = pssm.at('C')[i];
        if (pssm.find('G') != pssm.end()) flatPssm.scores[i * BASE_CODE_COUNT + BASE_G] = pssm.at('G')[i];
        if (pssm.find('T') != pssm.end()) flatPssm.scores[i * BASE_CODE_COUNT + BASE_T] = pssm.at('T')[i];
    }
    std::vector<std::uint8_t> codes;
    codes.reserve(window.size());
    for (char base : window) {
        codes.push_back(codeForBase(base));
    }
    return calculateScoreAt(codes, 0, flatPssm);
}

FlatPSSM flattenPSSM(const PSSM& pssm, bool skipN) {
    FlatPSSM flatPssm;
    flatPssm.motifLength = static_cast<size_t>(pssm.motifLength);
    flatPssm.scores.assign(flatPssm.motifLength * BASE_CODE_COUNT, skipN ? SENTINEL_SCORE : 0.0);

    const std::array<std::pair<char, std::uint8_t>, 4> bases{{
        {'A', BASE_A},
        {'C', BASE_C},
        {'G', BASE_G},
        {'T', BASE_T}
    }};
    for (const auto& [base, code] : bases) {
        auto it = pssm.pssm.find(base);
        if (it == pssm.pssm.end()) {
            continue;
        }
        for (size_t i = 0; i < flatPssm.motifLength && i < it->second.size(); ++i) {
            flatPssm.scores[i * BASE_CODE_COUNT + code] = it->second[i];
        }
    }

    return flatPssm;
}

std::string reverseComplement(const std::string& sequence) {
    std::string revComp = sequence;
    std::reverse(revComp.begin(), revComp.end());
    for (char& base : revComp) {
        base = complementBase(base);
    }
    return revComp;
}

ScoreRange scoreRangeForPSSM(const PSSM& pssm) {
    ScoreRange range;
    for (int i = 0; i < pssm.motifLength; ++i) {
        double columnMin = std::numeric_limits<double>::infinity();
        double columnMax = -std::numeric_limits<double>::infinity();
        for (const char base : {'A', 'C', 'G', 'T'}) {
            auto it = pssm.pssm.find(base);
            if (it == pssm.pssm.end() || static_cast<size_t>(i) >= it->second.size()) {
                continue;
            }
            const double score = it->second[i];
            if (!std::isfinite(score) || std::abs(score) >= SENTINEL_ABS_THRESHOLD) {
                continue;
            }
            columnMin = std::min(columnMin, score);
            columnMax = std::max(columnMax, score);
        }
        if (std::isfinite(columnMin)) range.min += columnMin;
        if (std::isfinite(columnMax)) range.max += columnMax;
    }
    return range;
}

double pwmRelativeScore(double score, const ScoreRange& scoreRange) {
    if (!std::isfinite(score) || score < -SENTINEL_ABS_THRESHOLD || scoreRange.max <= scoreRange.min) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return (score - scoreRange.min) / (scoreRange.max - scoreRange.min);
}

size_t outputStartForCoordinateMode(CoordinateMode coordinateMode, size_t windowStart) {
    return coordinateMode == CoordinateMode::Bed ? windowStart : windowStart + 1;
}

size_t outputEndForCoordinateMode(CoordinateMode coordinateMode, size_t windowStart, size_t motifLength) {
    return coordinateMode == CoordinateMode::Bed ? windowStart + motifLength : windowStart + 1 + motifLength;
}

bool ScoreBin::operator<(const ScoreBin& other) const {
    if (start != other.start) {
        return start < other.start;
    }
    return end < other.end;
}

ScoreDistribution::ScoreDistribution(const std::string& binWidthSpec)
    : minScore(std::numeric_limits<double>::infinity()),
      maxScore(-std::numeric_limits<double>::infinity()) {
    std::string normalizedSpec = binWidthSpec;
    std::transform(normalizedSpec.begin(), normalizedSpec.end(), normalizedSpec.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    if (normalizedSpec == "adaptive") {
        useAdaptiveBins = true;
        binScheme = "adaptive_log_ratio";
    } else {
        if (!parseDoubleStrict(binWidthSpec.c_str(), fixedBinWidth) || fixedBinWidth <= 0.0) {
            throw std::invalid_argument("Score distribution bin width must be 'adaptive' or greater than 0.");
        }
        binScheme = "fixed";
    }
}

double ScoreDistribution::binWidthForScore(double score) const {
    if (!useAdaptiveBins) {
        return fixedBinWidth;
    }
    if (score >= -10.0) {
        return 0.2;
    }
    if (score >= -50.0) {
        return 1.0;
    }
    if (score >= -250.0) {
        return 5.0;
    }
    if (score >= -1000.0) {
        return 10.0;
    }
    if (score >= -10000.0) {
        return 100.0;
    }
    return 500.0;
}

void ScoreDistribution::add(double score) {
    if (isSkippedScore(score)) {
        skippedWindows++;
        return;
    }
    const double binWidth = binWidthForScore(score);
    const double binStart = std::floor(score / binWidth) * binWidth;
    const ScoreBin bin{binStart, binStart + binWidth, binWidth};
    bins[bin]++;
    validWindows++;
    sumScore += score;
    if (score < minScore) minScore = score;
    if (score > maxScore) maxScore = score;
}

double ScoreDistribution::meanScore() const {
    if (validWindows == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return sumScore / validWindows;
}
