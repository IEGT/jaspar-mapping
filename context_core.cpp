#include "context_core.h"

#include "compressed_file_reader.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace motif_context {
namespace {

std::string firstField(const std::string& line) {
    std::istringstream input(line);
    std::string field;
    input >> field;
    return field;
}

std::string normalizedChromosome(std::string chromosome) {
    if (chromosome.size() > 3
        && (chromosome.compare(0, 3, "chr") == 0
            || chromosome.compare(0, 3, "CHR") == 0)) {
        chromosome.erase(0, 3);
    }
    return chromosome;
}

struct ChromosomeKey {
    int group = 1;
    int rank = 0;
    std::string text;
};

ChromosomeKey chromosomeKey(const std::string& chromosome) {
    const std::string normalized = normalizedChromosome(chromosome);
    if (!normalized.empty()
        && std::all_of(normalized.begin(), normalized.end(), [](const unsigned char value) {
               return value >= '0' && value <= '9';
           })) {
        try {
            return {0, std::stoi(normalized), {}};
        } catch (const std::exception&) {
            // Fall through to deterministic lexical ordering.
        }
    }
    if (normalized == "X") return {0, 23, {}};
    if (normalized == "Y") return {0, 24, {}};
    if (normalized == "M" || normalized == "MT") return {0, 25, {}};
    return {1, 0, normalized};
}

int compareChromosomes(const std::string& left, const std::string& right) {
    const ChromosomeKey leftKey = chromosomeKey(left);
    const ChromosomeKey rightKey = chromosomeKey(right);
    if (leftKey.group != rightKey.group) return leftKey.group < rightKey.group ? -1 : 1;
    if (leftKey.rank != rightKey.rank) return leftKey.rank < rightKey.rank ? -1 : 1;
    if (leftKey.text != rightKey.text) return leftKey.text < rightKey.text ? -1 : 1;
    return 0;
}

bool sameChromosome(const std::string& left, const std::string& right) {
    return normalizedChromosome(left) == normalizedChromosome(right);
}

bool comesBefore(const BedEntry& left, const BedEntry& right) {
    const int chromosomeComparison = compareChromosomes(left.chrom, right.chrom);
    if (chromosomeComparison != 0) return chromosomeComparison < 0;
    return left.start < right.start;
}

bool isBetterMatch(const BedEntry& candidate,
                   const ContextAnnotation& current,
                   const BedEntry& anchor) {
    if (!current.hasMatch || candidate.score > current.score) return true;
    if (candidate.score < current.score) return false;

    const double candidateShift = anchorOrientedCenterShift(anchor, candidate);
    const double candidateAbsoluteShift = std::abs(candidateShift);
    const double currentAbsoluteShift = std::abs(current.anchorOrientedCenterShift);
    if (candidateAbsoluteShift != currentAbsoluteShift) {
        return candidateAbsoluteShift < currentAbsoluteShift;
    }
    return candidateShift < current.anchorOrientedCenterShift;
}

std::int64_t flankTwiceChecked(const std::int64_t flankBp) {
    if (flankBp < 0) throw std::invalid_argument("context flank must be non-negative");
    if (flankBp > std::numeric_limits<std::int64_t>::max() / 2) {
        throw std::overflow_error("context flank is too large");
    }
    return flankBp * 2;
}

class BedStream {
public:
    BedStream(const std::string& filename, const bool verbose)
        : filename_(filename), reader_(filename, verbose) {}

    bool next(BedEntry& entry) {
        std::string line;
        while (reader_.getline(line)) {
            ++lineNumber_;
            if (line.empty() || line.front() == '#') continue;
            if (isBedHeader(line)) continue;

            std::string error;
            if (!parseBedLine(line, entry, error)) {
                throw std::runtime_error(filename_ + ":" + std::to_string(lineNumber_)
                                         + ": " + error);
            }
            if (previous_ && comesBefore(entry, *previous_)) {
                throw std::runtime_error(
                    filename_ + ":" + std::to_string(lineNumber_)
                    + ": input is not sorted by chromosome and start");
            }
            previous_ = entry;
            return true;
        }
        return false;
    }

private:
    std::string filename_;
    CompressedFileReader reader_;
    std::size_t lineNumber_ = 0;
    std::optional<BedEntry> previous_;
};

} // namespace

bool isBedHeader(const std::string& line) {
    const std::string field = firstField(line);
    return field == "Chr" || field == "Chromosome" || field == "chrom"
           || field == "chromosome";
}

bool parseBedLine(const std::string& line, BedEntry& entry, std::string& error) {
    std::istringstream input(line);
    std::string startText;
    std::string endText;
    std::string scoreText;
    std::string strandText;
    BedEntry parsed;

    if (!(input >> parsed.chrom >> startText >> endText >> parsed.name
          >> scoreText >> strandText)) {
        error = "expected at least six BED fields: chrom start end name score strand";
        return false;
    }

    try {
        std::size_t consumed = 0;
        parsed.start = std::stoll(startText, &consumed);
        if (consumed != startText.size()) throw std::invalid_argument("start");
        parsed.end = std::stoll(endText, &consumed);
        if (consumed != endText.size()) throw std::invalid_argument("end");
        parsed.score = std::stod(scoreText, &consumed);
        if (consumed != scoreText.size()) throw std::invalid_argument("score");
    } catch (const std::exception&) {
        error = "invalid start, end, or score field";
        return false;
    }

    if (parsed.start < 0 || parsed.end <= parsed.start) {
        error = "BED coordinates must satisfy 0 <= start < end";
        return false;
    }
    if (strandText.size() != 1
        || (strandText.front() != '+' && strandText.front() != '-'
            && strandText.front() != '.')) {
        error = "strand must be '+', '-', or '.'";
        return false;
    }

    parsed.strand = strandText.front();
    parsed.hasName = true;
    parsed.line = line;
    entry = std::move(parsed);
    return true;
}

std::int64_t centerTwice(const BedEntry& entry) {
    if (entry.start > std::numeric_limits<std::int64_t>::max() - entry.end) {
        throw std::overflow_error("BED interval center overflows int64");
    }
    return entry.start + entry.end;
}

double genomicCenterShift(const BedEntry& anchor, const BedEntry& neighbor) {
    return static_cast<double>(centerTwice(neighbor) - centerTwice(anchor)) / 2.0;
}

double anchorOrientedCenterShift(const BedEntry& anchor, const BedEntry& neighbor) {
    const double genomicShift = genomicCenterShift(anchor, neighbor);
    return anchor.strand == '-' ? -genomicShift : genomicShift;
}

bool isSameNamedLocus(const BedEntry& anchor, const BedEntry& neighbor) {
    return sameChromosome(anchor.chrom, neighbor.chrom)
           && anchor.start == neighbor.start
           && anchor.end == neighbor.end
           && anchor.name == neighbor.name;
}

ContextAnnotation summarizeNeighbors(const BedEntry& anchor,
                                     const std::vector<BedEntry>& neighbors,
                                     const std::int64_t flankBp,
                                     const bool excludeSameNamedLocus) {
    ContextAnnotation annotation;
    const std::int64_t flankTwice = flankTwiceChecked(flankBp);
    for (const BedEntry& neighbor : neighbors) {
        if (!sameChromosome(anchor.chrom, neighbor.chrom)) continue;
        if (std::llabs(centerTwice(neighbor) - centerTwice(anchor)) > flankTwice) continue;
        if (excludeSameNamedLocus && isSameNamedLocus(anchor, neighbor)) continue;

        ++annotation.count;
        if (!isBetterMatch(neighbor, annotation, anchor)) continue;

        annotation.hasMatch = true;
        annotation.genomicCenterShift = genomicCenterShift(anchor, neighbor);
        annotation.anchorOrientedCenterShift =
            anchorOrientedCenterShift(anchor, neighbor);
        annotation.score = neighbor.score;
        annotation.strandEqual = anchor.strand == neighbor.strand;
    }
    return annotation;
}

BedTable readBedTable(const std::string& filename, const bool verbose) {
    CompressedFileReader reader(filename, verbose);
    BedTable table;
    std::optional<BedEntry> previous;
    std::string line;
    std::size_t lineNumber = 0;
    while (reader.getline(line)) {
        ++lineNumber;
        if (line.empty() || line.front() == '#') continue;
        if (isBedHeader(line)) {
            if (table.header.empty()) table.header = line;
            continue;
        }

        BedEntry entry;
        std::string error;
        if (!parseBedLine(line, entry, error)) {
            throw std::runtime_error(filename + ":" + std::to_string(lineNumber)
                                     + ": " + error);
        }
        if (previous && comesBefore(entry, *previous)) {
            throw std::runtime_error(filename + ":" + std::to_string(lineNumber)
                                     + ": input is not sorted by chromosome and start");
        }
        if (previous && sameChromosome(previous->chrom, entry.chrom)
            && centerTwice(entry) < centerTwice(*previous)) {
            throw std::runtime_error(
                filename + ":" + std::to_string(lineNumber)
                + ": motif centers are not sorted; split variable-length anchors first");
        }
        previous = entry;
        table.entries.push_back(std::move(entry));
    }
    return table;
}

std::vector<ContextAnnotation> annotateFromBedFile(
    const std::vector<BedEntry>& anchors,
    const std::string& neighborFilename,
    const std::int64_t flankBp,
    const bool excludeSameNamedLocus,
    const bool verbose) {
    BedStream stream(neighborFilename, verbose);
    std::vector<ContextAnnotation> annotations;
    annotations.reserve(anchors.size());
    std::deque<BedEntry> active;
    std::optional<BedEntry> pending;
    const std::int64_t flankTwice = flankTwiceChecked(flankBp);

    for (const BedEntry& anchor : anchors) {
        active.erase(
            std::remove_if(active.begin(), active.end(), [&](const BedEntry& candidate) {
                const int chromosomeComparison = compareChromosomes(candidate.chrom, anchor.chrom);
                return chromosomeComparison < 0
                       || (chromosomeComparison == 0
                           && centerTwice(candidate) < centerTwice(anchor) - flankTwice);
            }),
            active.end());

        while (true) {
            BedEntry candidate;
            if (pending) {
                candidate = std::move(*pending);
                pending.reset();
            } else if (!stream.next(candidate)) {
                break;
            }

            const int chromosomeComparison = compareChromosomes(candidate.chrom, anchor.chrom);
            if (chromosomeComparison < 0) continue;
            const long double candidateStartTwice =
                static_cast<long double>(candidate.start) * 2.0L;
            const long double maximumStartTwice =
                static_cast<long double>(centerTwice(anchor)) + flankTwice;
            if (chromosomeComparison > 0
                || candidateStartTwice > maximumStartTwice) {
                pending = std::move(candidate);
                break;
            }
            active.push_back(std::move(candidate));
        }

        std::vector<BedEntry> candidates(active.begin(), active.end());
        annotations.push_back(
            summarizeNeighbors(anchor, candidates, flankBp, excludeSameNamedLocus));
    }
    return annotations;
}

std::string basenameWithoutContextSuffix(const std::string& filePath) {
    const std::size_t lastSlash = filePath.find_last_of("/\\");
    const std::string base = lastSlash == std::string::npos
                                 ? filePath
                                 : filePath.substr(lastSlash + 1);
    const std::size_t firstUnderscore = base.find('_');
    if (firstUnderscore == std::string::npos) return base;
    const std::size_t secondUnderscore = base.find('_', firstUnderscore + 1);
    return secondUnderscore == std::string::npos ? base : base.substr(0, secondUnderscore);
}

} // namespace motif_context
