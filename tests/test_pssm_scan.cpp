#include "../pssm_scan_core.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace {

int failures = 0;

const std::unordered_map<char, const double> testBackgroundFrequencies = {
    {'A', 0.25},
    {'C', 0.25},
    {'G', 0.25},
    {'T', 0.25},
    {'a', 0.25},
    {'c', 0.25},
    {'g', 0.25},
    {'t', 0.25},
    {'N', 1.0},
    {'n', 1.0}
};

void expectTrue(const bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << std::endl;
        failures++;
    }
}

void expectEqual(const std::string& actual, const std::string& expected, const std::string& message) {
    if (actual != expected) {
        std::cerr << "FAIL: " << message << " expected '" << expected
                  << "' got '" << actual << "'" << std::endl;
        failures++;
    }
}

template <typename T, typename U>
void expectEqual(const T& actual, const U& expected, const std::string& message) {
    if (actual != expected) {
        std::cerr << "FAIL: " << message << " expected " << expected
                  << " got " << actual << std::endl;
        failures++;
    }
}

void expectNear(const double actual, const double expected, const double tolerance, const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        std::cerr << "FAIL: " << message << " expected " << expected
                  << " got " << actual << std::endl;
        failures++;
    }
}

std::filesystem::path writeTinyPSSMFixture() {
    const auto stamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("jaspar_mapping_test_" + std::to_string(stamp) + ".pfm");
    std::ofstream out(path);
    out << ">MA0001.1 MotifOne\n";
    out << "A [ 1 0 3 ]\n";
    out << "C [ 0 2 0 ]\n";
    out << "G [ 1 1 0 ]\n";
    out << "T [ 0 0 1 ]\n";
    out << ">MA0002.1 MotifTwo\n";
    out << "A [ 2 2 ]\n";
    out << "C [ 0 0 ]\n";
    out << "G [ 0 0 ]\n";
    out << "T [ 2 2 ]\n";
    return path;
}

PSSM makeRangeTestPSSM() {
    pssm_type matrix;
    matrix['A'] = {2.0, 0.0};
    matrix['C'] = {1.0, 2.0};
    matrix['G'] = {0.0, 1.0};
    matrix['T'] = {1.0, 1.0};
    PSSM pssm(matrix, "TEST.1", "RangeTest", 2);
    pssm.normalizePSSM(testBackgroundFrequencies, "log2_relative_risk", 0.0);
    return pssm;
}

void testLogScores() {
    expectNear(PSSM::logRelativeRisk(0.5, 0.25), 1.0, 1e-12, "logRelativeRisk known value");
    expectNear(PSSM::logRelativeRiskACGT(0.25), 0.0, 1e-12, "logRelativeRiskACGT background value");
    expectNear(PSSM::logOddsRatioACGT(1.0, 4.0), 0.0, 1e-12, "logOddsRatioACGT background odds");
    expectNear(PSSM::logOddsRatioACGT(3.0, 4.0), std::log2(9.0), 1e-12, "logOddsRatioACGT enriched odds");
    expectEqual(PSSM::logRelativeRisk(0.0, 0.25), -1e9, "logRelativeRisk zero sentinel");
    expectEqual(PSSM::logRelativeRiskACGT(0.0), -1e9, "logRelativeRiskACGT zero sentinel");
    expectEqual(PSSM::logOddsRatioACGT(0.0, 4.0), -1e9, "logOddsRatioACGT zero sentinel");
    expectEqual(PSSM::logOddsRatioACGT(4.0, 4.0), 1e9, "logOddsRatioACGT fixed-nucleotide sentinel");
}

void testScoreModeAliases() {
    const std::vector<std::string> relativeRiskAliases = {
        "log2_relative_risk", "log_relative_risk", "relative_risk", "log_ratio",
        "log2_ratio", "log2-relative-risk", "LOG_RATIO"
    };
    for (const std::string& alias : relativeRiskAliases) {
        expectEqual(PSSM::canonicalScoreModeName(alias), "log2_relative_risk",
                    "relative-risk alias " + alias);
    }

    const std::vector<std::string> logOddsAliases = {
        "log_odds", "log2_odds", "log_odds_ratio", "odds_ratio",
        "log-odds", "LOG2_ODDS"
    };
    for (const std::string& alias : logOddsAliases) {
        expectEqual(PSSM::canonicalScoreModeName(alias), "log_odds",
                    "log-odds alias " + alias);
    }

    expectEqual(PSSM::canonicalScoreModeName("unknown"), "", "unknown score mode");
}

void testParsePSSMFile() {
    const std::filesystem::path fixture = writeTinyPSSMFixture();
    pssm_list_type motifs;
    expectEqual(PSSM::parsePSSMFile(fixture.string(), motifs, "", 0), 0, "parse tiny fixture");
    expectEqual(motifs.size(), static_cast<size_t>(2), "parsed motif count");
    const auto motifOne = motifs.find("MA0001.1");
    const auto motifTwo = motifs.find("MA0002.1");
    expectTrue(motifOne != motifs.end(), "first motif present");
    expectTrue(motifTwo != motifs.end(), "second motif present");
    if (motifOne != motifs.end() && motifTwo != motifs.end()) {
        expectEqual(motifOne->second.motifLength, 3, "first motif length");
        expectEqual(motifTwo->second.motifLength, 2, "second motif length");
        expectNear(motifOne->second.colsums[0], 2.0, 1e-12, "first colsum");
        expectNear(motifOne->second.colsums[1], 3.0, 1e-12, "second colsum");
        expectNear(motifOne->second.colsums[2], 4.0, 1e-12, "third colsum");
    }

    pssm_list_type targetMotifs;
    expectEqual(PSSM::parsePSSMFile(fixture.string(), targetMotifs, "MA0002.1", 0), 0,
                "parse target motif");
    expectEqual(targetMotifs.size(), static_cast<size_t>(1), "target motif count");
    expectTrue(targetMotifs.find("MA0002.1") != targetMotifs.end(), "target motif present");

    std::filesystem::remove(fixture);
}

void testNormalizePSSM() {
    pssm_type matrix;
    matrix['A'] = {1.0, 0.0};
    matrix['C'] = {0.0, 2.0};
    matrix['G'] = {1.0, 1.0};
    matrix['T'] = {0.0, 1.0};

    PSSM unsmoothed(matrix, "NORM.1", "Norm", 2);
    unsmoothed.normalizePSSM(testBackgroundFrequencies, "log2_relative_risk", 0.0);
    expectEqual(unsmoothed.pssm['A'][1], -1e9, "unsmoothed zero count sentinel");

    PSSM smoothed(matrix, "NORM.1", "Norm", 2);
    smoothed.normalizePSSM(testBackgroundFrequencies, "log2_relative_risk", 1.0);
    expectTrue(std::isfinite(smoothed.pssm['A'][1]), "smoothed zero count finite");
}

void testBaseCodesAndReverseComplement() {
    expectEqual(codeForBase('A'), BASE_A, "A code");
    expectEqual(codeForBase('c'), BASE_C, "lowercase c code");
    expectEqual(codeForBase('G'), BASE_G, "G code");
    expectEqual(codeForBase('t'), BASE_T, "lowercase t code");
    expectEqual(codeForBase('x'), BASE_N, "unknown base code");

    expectEqual(complementCode(BASE_A), BASE_T, "A complement code");
    expectEqual(complementCode(BASE_C), BASE_G, "C complement code");
    expectEqual(complementCode(BASE_G), BASE_C, "G complement code");
    expectEqual(complementCode(BASE_T), BASE_A, "T complement code");
    expectEqual(complementCode(BASE_N), BASE_N, "N complement code");

    expectEqual(reverseComplement("ACGTN"), "NACGT", "reverse complement");
    expectEqual(reverseComplement(reverseComplement("ACGTAC")), "ACGTAC",
                "reverse complement round trip");
}

void testCalculateScoreAt() {
    FlatPSSM flat;
    flat.motifLength = 3;
    flat.scores.assign(flat.motifLength * BASE_CODE_COUNT, 0.0);
    flat.scores[0 * BASE_CODE_COUNT + BASE_A] = 1.0;
    flat.scores[1 * BASE_CODE_COUNT + BASE_C] = 2.0;
    flat.scores[2 * BASE_CODE_COUNT + BASE_G] = 3.0;

    std::vector<std::uint8_t> codes = {
        codeForBase('T'), codeForBase('A'), codeForBase('C'), codeForBase('G'), codeForBase('T')
    };
    expectNear(calculateScoreAt(codes, 1, flat), 6.0, 1e-12, "calculateScoreAt offset window");
}

void testScoreRangeAndPWMRelativeScore() {
    PSSM pssm = makeRangeTestPSSM();
    const ScoreRange range = scoreRangeForPSSM(pssm);
    expectTrue(std::isfinite(range.min), "score range min finite");
    expectTrue(std::isfinite(range.max), "score range max finite");
    expectTrue(range.max > range.min, "score range non-empty");
    expectNear(range.min, 0.0, 1e-12, "score range ignores negative sentinels");
    expectNear(range.max, 2.0, 1e-12, "score range finite maximum");

    const double score = calculateScore("TC", pssm.pssm, true);
    const double relativeScore = pwmRelativeScore(score, range);
    expectTrue(std::isfinite(relativeScore), "relative score finite despite zero-count cells elsewhere");
    expectTrue(relativeScore >= 0.0 && relativeScore <= 1.0, "relative score in [0,1]");
    expectNear(relativeScore, 0.5, 1e-12, "relative score uses reachable finite range");
}

void testCoordinateModeOffsets() {
    const size_t windowStart = 10;
    const size_t motifLength = 16;
    expectEqual(outputStartForCoordinateMode(CoordinateMode::Legacy, windowStart),
                static_cast<size_t>(11), "legacy output start");
    expectEqual(outputEndForCoordinateMode(CoordinateMode::Legacy, windowStart, motifLength),
                static_cast<size_t>(27), "legacy output end");
    expectEqual(outputStartForCoordinateMode(CoordinateMode::Bed, windowStart),
                static_cast<size_t>(10), "BED output start");
    expectEqual(outputEndForCoordinateMode(CoordinateMode::Bed, windowStart, motifLength),
                static_cast<size_t>(26), "BED output end");
}

void testScoreDistribution() {
    ScoreDistribution distribution("adaptive");
    expectNear(distribution.binWidthForScore(-10.0), 0.2, 1e-12, "adaptive bin width at -10");
    expectNear(distribution.binWidthForScore(-10.1), 1.0, 1e-12, "adaptive bin width below -10");
    expectNear(distribution.binWidthForScore(-50.0), 1.0, 1e-12, "adaptive bin width at -50");
    expectNear(distribution.binWidthForScore(-50.1), 5.0, 1e-12, "adaptive bin width below -50");
    expectNear(distribution.binWidthForScore(-250.0), 5.0, 1e-12, "adaptive bin width at -250");
    expectNear(distribution.binWidthForScore(-250.1), 10.0, 1e-12, "adaptive bin width below -250");
    expectNear(distribution.binWidthForScore(-1000.0), 10.0, 1e-12, "adaptive bin width at -1000");
    expectNear(distribution.binWidthForScore(-1000.1), 100.0, 1e-12, "adaptive bin width below -1000");
    expectNear(distribution.binWidthForScore(-10000.0), 100.0, 1e-12, "adaptive bin width at -10000");
    expectNear(distribution.binWidthForScore(-10000.1), 500.0, 1e-12, "adaptive bin width below -10000");

    distribution.add(SENTINEL_SCORE);
    expectEqual(distribution.skippedWindows, static_cast<std::uint64_t>(1), "sentinel score skipped");
    expectEqual(distribution.validWindows, static_cast<std::uint64_t>(0), "sentinel score not valid");

    distribution.add(-9.95);
    expectEqual(distribution.validWindows, static_cast<std::uint64_t>(1), "normal score counted");
    const ScoreBin expectedBin{-10.0, -9.8, 0.2};
    expectTrue(distribution.bins.find(expectedBin) != distribution.bins.end(), "normal score adaptive bin");
    expectEqual(distribution.bins[expectedBin], static_cast<std::uint64_t>(1), "normal score bin count");
}

} // namespace

int main() {
    testLogScores();
    testScoreModeAliases();
    testParsePSSMFile();
    testNormalizePSSM();
    testBaseCodesAndReverseComplement();
    testCalculateScoreAt();
    testScoreRangeAndPWMRelativeScore();
    testCoordinateModeOffsets();
    testScoreDistribution();

    if (failures != 0) {
        std::cerr << failures << " test failure(s)" << std::endl;
        return 1;
    }

    std::cout << "All pssm_scan unit tests passed." << std::endl;
    return 0;
}
