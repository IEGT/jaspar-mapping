#include "pssm_scan_core.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

char complementBase(const char& base) {
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

std::uint8_t codeForBase(const char& base) {
    return baseCodeTable()[static_cast<unsigned char>(base)];
}

std::uint8_t complementCode(const std::uint8_t& code) {
    switch (code) {
        case BASE_A: return BASE_T;
        case BASE_C: return BASE_G;
        case BASE_G: return BASE_C;
        case BASE_T: return BASE_A;
        default: return BASE_N;
    }
}

bool isSkippedScore(const double& score) {
    return !std::isfinite(score) || score < SENTINEL_SCORE / 10.0;
}

double calculateScoreAt(const std::vector<std::uint8_t>& codes, const size_t& start, const FlatPSSM& pssm) {
    double score = 0.0;

    for (size_t i = 0; i < pssm.motifLength; ++i) {
        score += pssm.scores[i * BASE_CODE_COUNT + codes[start + i]];
    }

    return score;
}

double calculateScoreAtGenomicStart(const std::vector<std::uint8_t>& codes, const size_t& sequenceLength,
                                    const bool reverseComplementWindow, const size_t& genomicStart,
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

ScoreBlock calculateScoreBlock(const std::vector<std::uint8_t>& codes, const size_t& sequenceLength,
                               const bool reverseComplementWindow, const size_t& blockStart,
                               const size_t& windowCount, const FlatPSSM& pssm) {
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
    for (const char& base : window) {
        codes.push_back(codeForBase(base));
    }
    return calculateScoreAt(codes, 0, flatPssm);
}

FlatPSSM flattenPSSM(const PSSM& pssm, const bool skipN) {
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

double pwmRelativeScore(const double& score, const ScoreRange& scoreRange) {
    if (!std::isfinite(score) || score < -SENTINEL_ABS_THRESHOLD || scoreRange.max <= scoreRange.min) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return (score - scoreRange.min) / (scoreRange.max - scoreRange.min);
}

size_t outputStartForCoordinateMode(const CoordinateMode& coordinateMode, const size_t& windowStart) {
    return coordinateMode == CoordinateMode::Bed ? windowStart : windowStart + 1;
}

size_t outputEndForCoordinateMode(const CoordinateMode& coordinateMode, const size_t& windowStart, const size_t& motifLength) {
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

double ScoreDistribution::binWidthForScore(const double& score) const {
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

void ScoreDistribution::add(const double& score) {
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
