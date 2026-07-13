#include "../context_core.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int failures = 0;

void expect(const bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

void expectNear(const double actual, const double expected,
                const std::string& message) {
    expect(std::abs(actual - expected) < 1e-12,
           message + " (expected " + std::to_string(expected)
               + ", got " + std::to_string(actual) + ")");
}

motif_context::BedEntry entry(const std::string& chrom, const std::int64_t start,
                              const std::int64_t end, const std::string& name,
                              const double score, const char strand) {
    motif_context::BedEntry result;
    result.hasName = true;
    result.chrom = chrom;
    result.start = start;
    result.end = end;
    result.name = name;
    result.score = score;
    result.strand = strand;
    return result;
}

class TemporaryFile {
public:
    explicit TemporaryFile(const std::string& content) {
        path_ = std::filesystem::temp_directory_path()
                / ("jaspar-context-test-" + std::to_string(std::rand()) + ".bed");
        std::ofstream output(path_);
        output << content;
        if (!output) throw std::runtime_error("failed to create context test fixture");
    }

    ~TemporaryFile() {
        std::error_code error;
        std::filesystem::remove(path_, error);
    }

    const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
};

void testParsingAndCenters() {
    expect(motif_context::isBedHeader("Chromosome\tFrom\tTo\tName\tScore\tStrand"),
           "Chromosome header recognized");
    expect(!motif_context::isBedHeader("Chr1\t10\t20\tM\t1\t+"),
           "Chr-prefixed data is not mistaken for a header");

    motif_context::BedEntry parsed;
    std::string error;
    expect(motif_context::parseBedLine("chr1 100 116 TP73 -Inf -", parsed, error),
           "BED6 parser accepts an infinite score");
    expect(parsed.start == 100 && parsed.end == 116 && parsed.strand == '-',
           "BED6 fields parsed correctly");
    expect(std::isinf(parsed.score) && parsed.score < 0,
           "negative infinity remains a score");

    const auto anchor = entry("1", 100, 116, "TP73", 5.0, '+');
    const auto neighbor = entry("1", 125, 135, "SP1", 2.0, '-');
    expect(motif_context::centerTwice(anchor) == 216, "twice-center is exact");
    expectNear(motif_context::genomicCenterShift(anchor, neighbor), 22.0,
               "center-to-center displacement");

    const auto reverseAnchor = entry("1", 100, 116, "TP73", 5.0, '-');
    expectNear(motif_context::anchorOrientedCenterShift(reverseAnchor, neighbor), -22.0,
               "minus anchor flips displacement");
}

void testSummary() {
    const auto anchor = entry("1", 100, 116, "TP73", 5.0, '+');
    const std::vector<motif_context::BedEntry> neighbors = {
        entry("1", 100, 116, "TP73", 99.0, '-'),
        entry("1", 120, 136, "TP73", 2.0, '+'),
        entry("1", 70, 106, "TP73", 3.0, '-'),
        entry("1", 120, 140, "SP1", 50.0, '+'),
    };

    const motif_context::ContextAnnotation annotation =
        motif_context::summarizeNeighbors(anchor, neighbors, 20, true);
    expect(annotation.count == 2,
           "only non-self centers within the flank are counted");
    expect(annotation.hasMatch, "eligible neighbor produces an annotation");
    expectNear(annotation.genomicCenterShift, -20.0,
               "highest-scoring eligible neighbor supplies genomic shift");
    expectNear(annotation.anchorOrientedCenterShift, -20.0,
               "plus anchor retains genomic direction");
    expect(!annotation.strandEqual, "relative orientation is retained");
    expectNear(annotation.score, 3.0, "best eligible score retained");
}

void testStreamingAnnotation() {
    const std::vector<motif_context::BedEntry> anchors = {
        entry("1", 100, 116, "TP73", 5.0, '+'),
        entry("1", 200, 216, "TP73", 5.0, '-'),
    };
    const TemporaryFile neighbors(
        "Chromosome\tFrom\tTo\tName\tScore\tStrand\n"
        "1\t100\t116\tTP73\t100\t-\n"
        "1\t120\t136\tTP73\t2\t+\n"
        "1\t180\t196\tTP73\t3\t+\n");

    const auto annotations = motif_context::annotateFromBedFile(
        anchors, neighbors.path().string(), 20, true, false);
    expect(annotations.size() == 2, "one annotation returned per anchor");
    expect(annotations[0].count == 1, "exact anchor locus excluded from tandem count");
    expectNear(annotations[0].anchorOrientedCenterShift, 20.0,
               "forward tandem displacement");
    expect(annotations[1].count == 1, "reverse anchor finds upstream genomic partner");
    expectNear(annotations[1].genomicCenterShift, -20.0,
               "reverse tandem genomic displacement");
    expectNear(annotations[1].anchorOrientedCenterShift, 20.0,
               "reverse tandem oriented displacement");
    expect(!annotations[1].strandEqual, "opposite tandem orientation retained");

    const TemporaryFile equalStartNeighbors(
        "1\t120\t150\tSP1\t1\t+\n"
        "1\t120\t136\tSP1\t2\t+\n");
    const auto equalStartAnnotations = motif_context::annotateFromBedFile(
        {anchors.front()}, equalStartNeighbors.path().string(), 30, true, false);
    expect(equalStartAnnotations.front().count == 2,
           "BED rows tied on start need no secondary end/strand ordering");
}

} // namespace

int main() {
    testParsingAndCenters();
    testSummary();
    testStreamingAnnotation();
    if (failures != 0) {
        std::cerr << failures << " context test(s) failed.\n";
        return 1;
    }
    std::cout << "Context core tests passed.\n";
    return 0;
}
