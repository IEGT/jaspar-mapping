#include "../pssm_scan_core.h"

#include <cmath>
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

void testGenomicScoreBlocks() {
    FlatPSSM flat;
    flat.motifLength = 2;
    flat.scores.assign(flat.motifLength * BASE_CODE_COUNT, 0.0);
    flat.scores[0 * BASE_CODE_COUNT + BASE_A] = 1.0;
    flat.scores[0 * BASE_CODE_COUNT + BASE_C] = 2.0;
    flat.scores[0 * BASE_CODE_COUNT + BASE_G] = 3.0;
    flat.scores[0 * BASE_CODE_COUNT + BASE_T] = 4.0;
    flat.scores[0 * BASE_CODE_COUNT + BASE_N] = SENTINEL_SCORE;
    flat.scores[1 * BASE_CODE_COUNT + BASE_A] = 10.0;
    flat.scores[1 * BASE_CODE_COUNT + BASE_C] = 20.0;
    flat.scores[1 * BASE_CODE_COUNT + BASE_G] = 30.0;
    flat.scores[1 * BASE_CODE_COUNT + BASE_T] = 40.0;
    flat.scores[1 * BASE_CODE_COUNT + BASE_N] = SENTINEL_SCORE;

    std::vector<std::uint8_t> plusCodes = {
        codeForBase('A'), codeForBase('A'), codeForBase('C'), codeForBase('G')
    };
    std::vector<std::uint8_t> minusCodes(plusCodes.size());
    for (size_t i = 0; i < plusCodes.size(); ++i) {
        minusCodes[plusCodes.size() - 1 - i] = complementCode(plusCodes[i]);
    }

    expectNear(calculateScoreAtGenomicStart(plusCodes, plusCodes.size(), false, 1, flat),
               21.0, 1e-12, "plus genomic-start score");
    expectNear(calculateScoreAtGenomicStart(minusCodes, plusCodes.size(), true, 1, flat),
               43.0, 1e-12, "minus genomic-start score");

    const ScoreBlock plusBlock = calculateScoreBlock(plusCodes, plusCodes.size(), false, 0, 3, flat);
    expectEqual(plusBlock.blockStart, static_cast<size_t>(0), "score block start");
    expectEqual(plusBlock.scores.size(), static_cast<size_t>(3), "score block size");
    expectNear(plusBlock.scores[0], 11.0, 1e-12, "first dense block score");
    expectNear(plusBlock.scores[1], 21.0, 1e-12, "second dense block score");
    expectNear(plusBlock.scores[2], 32.0, 1e-12, "third dense block score");
    expectEqual(plusBlock.sequenceValid.size(), static_cast<size_t>(3), "score block validity size");
    expectTrue(plusBlock.sequenceValid[0] && plusBlock.sequenceValid[1] && plusBlock.sequenceValid[2],
               "ordinary dense windows have valid sequence");
    expectEqual(plusBlock.validWindows, static_cast<std::uint64_t>(3), "score block valid count");
    expectEqual(plusBlock.skippedWindows, static_cast<std::uint64_t>(0), "score block skipped count");

    std::vector<std::uint8_t> codesWithN = {codeForBase('A'), codeForBase('N')};
    const ScoreBlock skippedBlock = calculateScoreBlock(codesWithN, codesWithN.size(), false, 0, 1, flat);
    expectTrue(isSkippedScore(skippedBlock.scores[0]), "N-containing dense score skipped");
    expectTrue(!skippedBlock.sequenceValid[0], "N-containing dense sequence invalid");
    expectEqual(skippedBlock.validWindows, static_cast<std::uint64_t>(0), "skipped block valid count");
    expectEqual(skippedBlock.skippedWindows, static_cast<std::uint64_t>(1), "skipped block skipped count");

    FlatPSSM zeroCountFlat = flat;
    zeroCountFlat.scores[0 * BASE_CODE_COUNT + BASE_A] = SENTINEL_SCORE;
    const std::vector<std::uint8_t> zeroCountCodes = {codeForBase('A'), codeForBase('C')};
    const ScoreBlock zeroCountBlock = calculateScoreBlock(
        zeroCountCodes, zeroCountCodes.size(), false, 0, 1, zeroCountFlat);
    expectTrue(isSkippedScore(zeroCountBlock.scores[0]), "zero-count motif score is sentinel");
    expectTrue(zeroCountBlock.sequenceValid[0], "zero-count motif score retains valid DNA sequence");
    expectEqual(zeroCountBlock.validWindows, static_cast<std::uint64_t>(1),
                "zero-count motif score remains a dense candidate");
    expectEqual(zeroCountBlock.skippedWindows, static_cast<std::uint64_t>(0),
                "zero-count motif score is not a missing DNA window");
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

void testAnchorContextMaximum() {
    const ContextStartRange ordinary = anchorContextStartRange(
        100, 20, 25, 2, 4);
    expectTrue(!ordinary.empty, "ordinary anchor context is non-empty");
    expectEqual(ordinary.first, static_cast<size_t>(14),
                "context start includes upstream motif span");
    expectEqual(ordinary.last, static_cast<size_t>(27),
                "context start includes downstream flank");

    const ContextStartRange leftEdge = anchorContextStartRange(20, 0, 3, 5, 4);
    expectEqual(leftEdge.first, static_cast<size_t>(0),
                "context start clips at chromosome start");
    expectEqual(leftEdge.last, static_cast<size_t>(8),
                "left-edge context preserves downstream reach");

    const ContextStartRange rightEdge = anchorContextStartRange(20, 18, 20, 5, 4);
    expectEqual(rightEdge.first, static_cast<size_t>(9),
                "right-edge context preserves upstream reach");
    expectEqual(rightEdge.last, static_cast<size_t>(16),
                "context start clips at last complete motif window");
    expectTrue(anchorContextStartRange(20, 5, 5, 2, 4).empty,
               "empty anchor interval is rejected");

    FlatPSSM flat;
    flat.motifLength = 2;
    flat.scores.assign(flat.motifLength * BASE_CODE_COUNT, -3.0);
    flat.scores[0 * BASE_CODE_COUNT + BASE_T] = 5.0;
    flat.scores[1 * BASE_CODE_COUNT + BASE_T] = 7.0;
    flat.scores[0 * BASE_CODE_COUNT + BASE_N] = SENTINEL_SCORE;
    flat.scores[1 * BASE_CODE_COUNT + BASE_N] = SENTINEL_SCORE;

    const std::string sequence = "AAAATAAA";
    std::vector<std::uint8_t> plusCodes;
    plusCodes.reserve(sequence.size());
    for (const char base : sequence) {
        plusCodes.push_back(codeForBase(base));
    }
    std::vector<std::uint8_t> minusCodes(plusCodes.size());
    for (size_t i = 0; i < plusCodes.size(); ++i) {
        minusCodes[plusCodes.size() - 1 - i] = complementCode(plusCodes[i]);
    }

    const ContextMaximum maximum = maximumScoreInAnchorContext(
        plusCodes, minusCodes, sequence.size(), 3, 5, 0, flat,
        true, true, -20.0);
    expectTrue(maximum.available, "both-strand context maximum is available");
    expectNear(maximum.score, 12.0, 1e-12,
               "both-strand context maximum uses genomic coordinates");
    expectEqual(maximum.validWindows, static_cast<std::uint64_t>(10),
                "all in-context windows on both strands are inspected");
    expectEqual(maximum.skippedWindows, static_cast<std::uint64_t>(0),
                "ordinary context has no skipped windows");

    const ContextMaximum suppressed = maximumScoreInAnchorContext(
        plusCodes, minusCodes, sequence.size(), 3, 5, 0, flat,
        true, true, 13.0);
    expectTrue(!suppressed.available,
               "source floor suppresses lower context maxima");
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

void testOutputNaming() {
    const HitOutputOptions options{
        "log2_relative_risk",
        1.0,
        true,
        0.0,
        0.8,
        std::numeric_limits<double>::infinity(),
        CoordinateMode::Bed,
        true,
        true,
        "jaspar2026_core_nonredundant",
        "homo_sapiens_grch38_ensembl113_primary",
        "uniform_acgt_v1",
        "additive_per_base"
    };
    const std::filesystem::path outputDirectory = hitOutputDirectory("out", options);
    expectEqual(
        outputDirectory.generic_string(),
        "out/hits/score_mode-log2_relative_risk/pseudocount-1/threshold-0/"
        "pwm_relative_min-0.8/pwm_relative_max-none/coordinate_mode-bed/"
        "sequence-included/n_policy-skip",
        "hit output path records all content-changing options"
    );

    HitOutputOptions noThreshold = options;
    noThreshold.thresholdSet = false;
    expectTrue(hitOutputDirectory("out", noThreshold) != outputDirectory,
               "explicit threshold zero has a distinct output path");

    HitOutputOptions differentPseudocount = options;
    differentPseudocount.pseudocount = 0.0;
    expectTrue(hitOutputDirectory("out", differentPseudocount) != outputDirectory,
               "pseudocount changes the output path");

    HitOutputOptions legacyCoordinates = options;
    legacyCoordinates.coordinateMode = CoordinateMode::Legacy;
    expectTrue(hitOutputDirectory("out", legacyCoordinates) != outputDirectory,
               "coordinate mode changes the output path");

    HitOutputOptions noRelativeScoreFilter = options;
    noRelativeScoreFilter.minPwmRelativeScore =
        -std::numeric_limits<double>::infinity();
    expectTrue(hitOutputDirectory("out", noRelativeScoreFilter) != outputDirectory,
               "PWM-relative score bounds change the output path");

    expectEqual(
        sparseHitParquetOutputPath(
            "out", "MA0861.2", options, "1", "+", 100, 200).generic_string(),
        "out/tables/jaspar2026/motif_hit/"
        "motif_set_id=jaspar2026_core_nonredundant/"
        "genome_id=homo_sapiens_grch38_ensembl113_primary/"
        "motif_id=MA0861.2/score_mode=log2_relative_risk/pseudocount=1/"
        "background_model_id=uniform_acgt_v1/"
        "pseudocount_scheme=additive_per_base/minimum_score=0/"
        "minimum_pwm_relative_score=0.8/maximum_pwm_relative_score=none/"
        "chrom=1/strand=plus/n_policy=skip/matched_sequence=included/"
        "part-from=100-to=200-000000.parquet",
        "sparse Parquet path moves constant hit identity into partitions"
    );
    HitOutputOptions mouseGenome = options;
    mouseGenome.genomeID = "mus_musculus_grcm39_ensembl113_primary";
    expectTrue(
        sparseHitParquetOutputPath("out", "MA0861.2", mouseGenome,
                                   "1", "+", 100, 200) !=
        sparseHitParquetOutputPath("out", "MA0861.2", options,
                                   "1", "+", 100, 200),
        "genome identity changes the sparse Parquet path");
    expectEqual(denseScorePartFilename(100, 200, true, ".parquet"),
                "part-from=100-to=200-n_policy=skip-000000.parquet",
                "dense range and N policy in part filename");
    expectEqual(denseScorePartFilename(-1, -1, false, ".tsv"),
                "part-from=0-to=end-n_policy=neutral-000000.tsv",
                "dense full-range part filename");
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
    testGenomicScoreBlocks();
    testScoreRangeAndPWMRelativeScore();
    testAnchorContextMaximum();
    testCoordinateModeOffsets();
    testOutputNaming();
    testScoreDistribution();

    if (failures != 0) {
        std::cerr << failures << " test failure(s)" << std::endl;
        return 1;
    }

    std::cout << "All pssm_scan unit tests passed." << std::endl;
    return 0;
}
