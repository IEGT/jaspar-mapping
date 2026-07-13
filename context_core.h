#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace motif_context {

// BED coordinates are 0-based and half-open. Distances below compare the
// centers of the scored sequence spans, not an inferred TF-DNA footprint.
struct BedEntry {
    bool hasName = false;
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::string name;
    double score = 0.0;
    std::string line;
    char strand = '.';
};

struct BedTable {
    std::string header;
    std::vector<BedEntry> entries;
};

struct ContextAnnotation {
    bool hasMatch = false;
    // Genomic shift is positive toward increasing reference coordinates;
    // anchor-oriented shift reverses that sign for a minus-strand anchor.
    double genomicCenterShift = 0.0;
    double anchorOrientedCenterShift = 0.0;
    double score = 0.0;
    bool strandEqual = false;
    std::uint64_t count = 0;
};

bool isBedHeader(const std::string& line);
bool parseBedLine(const std::string& line, BedEntry& entry, std::string& error);
// Twice-center arithmetic preserves half-base centers without rounding.
std::int64_t centerTwice(const BedEntry& entry);
double genomicCenterShift(const BedEntry& anchor, const BedEntry& neighbor);
double anchorOrientedCenterShift(const BedEntry& anchor, const BedEntry& neighbor);
bool isSameNamedLocus(const BedEntry& anchor, const BedEntry& neighbor);

ContextAnnotation summarizeNeighbors(const BedEntry& anchor,
                                     const std::vector<BedEntry>& neighbors,
                                     std::int64_t flankBp,
                                     bool excludeSameNamedLocus = true);

BedTable readBedTable(const std::string& filename, bool verbose = false);

std::vector<ContextAnnotation> annotateFromBedFile(
    const std::vector<BedEntry>& anchors,
    const std::string& neighborFilename,
    std::int64_t flankBp,
    bool excludeSameNamedLocus = true,
    bool verbose = false);

std::string basenameWithoutContextSuffix(const std::string& filePath);

} // namespace motif_context
