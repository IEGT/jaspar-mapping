#ifndef PSSM_SCAN_CORE_H
#define PSSM_SCAN_CORE_H

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "pssm.h"

static constexpr std::uint8_t BASE_A = 0;
static constexpr std::uint8_t BASE_C = 1;
static constexpr std::uint8_t BASE_G = 2;
static constexpr std::uint8_t BASE_T = 3;
static constexpr std::uint8_t BASE_N = 4;
static constexpr size_t BASE_CODE_COUNT = 5;
static constexpr double SENTINEL_SCORE = -1e9;
static constexpr double SENTINEL_ABS_THRESHOLD = 1e8;

enum class CoordinateMode {
    Legacy,
    Bed
};

struct HitOutputOptions {
    std::string scoreMode;
    double pseudocount;
    bool thresholdSet;
    double threshold;
    double minPwmRelativeScore;
    double maxPwmRelativeScore;
    CoordinateMode coordinateMode;
    bool showSequence;
    bool skipN;
};

struct ScoreRange {
    double min = 0.0;
    double max = 0.0;
};

struct FlatPSSM {
    size_t motifLength = 0;
    std::vector<double> scores;
};

struct ScoreBlock {
    size_t blockStart = 0;
    std::vector<double> scores;
    // Separates unavailable sequence (for example N under --skip-N) from a
    // valid alignment whose unsmoothed motif score is negative infinity.
    std::vector<bool> sequenceValid;
    std::uint64_t validWindows = 0;
    std::uint64_t skippedWindows = 0;
};

struct ScoreBin {
    double start;
    double end;
    double width;

    bool operator<(const ScoreBin& other) const;
};

struct ScoreDistribution {
    explicit ScoreDistribution(const std::string& binWidthSpec);

    bool useAdaptiveBins = false;
    double fixedBinWidth = 1.0;
    std::string binScheme = "fixed";
    std::uint64_t validWindows = 0;
    std::uint64_t skippedWindows = 0;
    double minScore;
    double maxScore;
    double sumScore = 0.0;
    std::map<ScoreBin, std::uint64_t> bins;

    double binWidthForScore(double score) const;
    void add(double score);
    double meanScore() const;
};

char complementBase(char base);
bool parseDoubleStrict(const char* text, double& value);
std::string formatDoubleForFileLabel(double value);
std::filesystem::path hitOutputDirectory(const std::filesystem::path& outdir,
                                         const HitOutputOptions& options);
std::string motifDatasetLabelFromPSSMFile(const std::filesystem::path& pssmFile);
std::string denseScorePartFilename(long from, long to, bool skipN,
                                   const std::string& extension);
std::uint8_t codeForBase(char base);
std::uint8_t complementCode(std::uint8_t code);
bool isSkippedScore(double score);
double calculateScoreAt(const std::vector<std::uint8_t>& codes, size_t start, const FlatPSSM& pssm);
double calculateScoreAtGenomicStart(const std::vector<std::uint8_t>& codes, size_t sequenceLength,
                                    bool reverseComplementWindow, size_t genomicStart,
                                    const FlatPSSM& pssm);
ScoreBlock calculateScoreBlock(const std::vector<std::uint8_t>& codes, size_t sequenceLength,
                               bool reverseComplementWindow, size_t blockStart,
                               size_t windowCount, const FlatPSSM& pssm);
double calculateScore(const std::string& window, const pssm_type& pssm, bool skipN);
FlatPSSM flattenPSSM(const PSSM& pssm, bool skipN);
std::string reverseComplement(const std::string& sequence);
ScoreRange scoreRangeForPSSM(const PSSM& pssm);
double pwmRelativeScore(double score, const ScoreRange& scoreRange);
size_t outputStartForCoordinateMode(CoordinateMode coordinateMode, size_t windowStart);
size_t outputEndForCoordinateMode(CoordinateMode coordinateMode, size_t windowStart, size_t motifLength);

#endif
