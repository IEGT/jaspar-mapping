#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "../gtf_file_region.h"

namespace {

int failures = 0;

void expect(const bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

void expectRegion(const GeneRegion& region, const size_t expectedStart,
                  const size_t expectedEnd, const std::string& message) {
    expect(region.start == expectedStart,
           message + " start (expected " + std::to_string(expectedStart) +
               ", got " + std::to_string(region.start) + ")");
    expect(region.end == expectedEnd,
           message + " end (expected " + std::to_string(expectedEnd) +
               ", got " + std::to_string(region.end) + ")");
    expect(region.end >= region.start &&
               region.end - region.start == expectedEnd - expectedStart,
           message + " width");
}

}  // namespace

int main() {
    const auto unique = std::chrono::high_resolution_clock::now()
                            .time_since_epoch()
                            .count();
    const auto fixture = std::filesystem::temp_directory_path() /
                         ("jaspar_mapping_gtf_" + std::to_string(unique) + ".gtf");

    {
        std::ofstream out(fixture);
        out << "1\ttest\ttranscript\t1000\t2000\t.\t+\t.\t"
               "gene_id \"ENSG_PLUS\"; gene_name \"PLUS\"; "
               "transcript_name \"PLUS-201\";\n";
        out << "2\ttest\ttranscript\t3000\t4000\t.\t-\t.\t"
               "gene_id \"ENSG_MINUS\"; gene_name \"MINUS\"; "
               "transcript_name \"MINUS-201\";\n";
    }

    const auto regions = parseGTFFile(fixture.string());
    std::filesystem::remove(fixture);

    const auto plusIt = regions.find("PLUS-201");
    const auto minusIt = regions.find("MINUS-201");
    expect(plusIt != regions.end(), "plus-strand transcript was parsed");
    expect(minusIt != regions.end(), "minus-strand transcript was parsed");

    if (plusIt != regions.end()) {
        const GeneRegion& plus = plusIt->second;
        expectRegion(plus, 999, 2000, "GTF-to-BED conversion on plus strand");
        expect(plus.toBedString() == "1\t999\t2000\tPLUS\t0\t+",
               "BED serialization uses converted coordinates");
        expectRegion(plus.relative_upstream(1, 500), 499, 999,
                     "500 bp plus-strand promoter");
        expectRegion(plus.relative_downstream(1, 500), 999, 1499,
                     "500 bp plus-strand downstream region");
        expectRegion(plus.relative_upstream(10, 20), 979, 990,
                     "inclusive plus-strand upstream distances");
    }

    if (minusIt != regions.end()) {
        const GeneRegion& minus = minusIt->second;
        expectRegion(minus, 2999, 4000, "GTF-to-BED conversion on minus strand");
        expect(minus.toBedString() == "2\t2999\t4000\tMINUS\t0\t-",
               "minus-strand BED serialization uses converted coordinates");
        expectRegion(minus.relative_upstream(1, 500), 4000, 4500,
                     "500 bp minus-strand promoter");
        expectRegion(minus.relative_downstream(1, 500), 3500, 4000,
                     "500 bp minus-strand downstream region");
        expectRegion(minus.relative_upstream(10, 20), 4009, 4020,
                     "inclusive minus-strand upstream distances");
    }

    const GeneRegion chromosomeStart{"1", 99, 200, "+", "ENSG_EDGE",
                                     "EDGE", "EDGE-201"};
    expectRegion(chromosomeStart.relative_upstream(1, 500), 0, 99,
                 "promoter clipped at chromosome start");

    if (failures != 0) {
        std::cerr << failures << " GTF region test(s) failed.\n";
        return 1;
    }

    std::cout << "All GTF region tests passed.\n";
    return 0;
}
