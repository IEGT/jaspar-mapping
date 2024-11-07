#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <tuple>
#include <limits>
#include <cmath>

/** \brief Representation of single line in .bed file */
struct BedEntry {
    std::string chrom;
    int start;
    int end;
    double score;
    std::string line;
    char strand;

    BedEntry(const std::string& chrom, int start, int end, double score, const char strand = '.', const std::string& line = "")
        : chrom(chrom), start(start), end(end), score(score), strand(strand), line(line) {}
};

/** \brief Splits line of .bed file and retrieves values of relevance 
 * @param line - line from .bed file
 * @param entry - BedEntry object to store values
 */
bool parseBedLine(const std::string& line, BedEntry& entry) {
    std::istringstream iss(line);
    std::string chrom;
    int start, end;
    double score;
    char strand;
    if (iss >> chrom >> start >> end >> score >> strand) {
        entry = BedEntry(chrom, start, end, score, strand, line);
        return true;
    }
    return false;
}

/** \brief Process main .bed file and list of other .bed files
 * @param mainBedFile - path to main .bed file
 * @param otherBedFiles - list of paths to other .bed files
 */
void processBedFiles(const std::string& mainBedFile, const std::vector<std::string>& otherBedFiles) {
    std::ifstream mainFile(mainBedFile);
    if (!mainFile.is_open()) {
        std::cerr << "Error opening main .bed file: " << mainBedFile << std::endl;
        return;
    }

    /** open all other .bed files and add to otherFiles vector */
    std::vector<std::ifstream> otherFiles;
    for (const auto& file : otherBedFiles) {
        otherFiles.emplace_back(file);
        if (!otherFiles.back().is_open()) {
            std::cerr << "Error opening .bed file: " << file << std::endl;
            return;
        }
    }

    /** Create all the queues that buffer the lines in the .bed files. */
    std::vector<std::queue<BedEntry>> queues(otherFiles.size());

    /** Iterate over all lines in the main .bed file and annotate it. */
    std::string line;
    while (std::getline(mainFile, line)) {

        /** Read next entry of main .bed file. */
        BedEntry mainEntry("", 0, 0, 0.0, '.', line);
        if (!parseBedLine(line, mainEntry)) {
            std::cerr << "Error parsing line in main .bed file: " << line << std::endl;
            continue;
        }

        /** Store the best annotation for each other file. */
        std::vector<std::tuple<int, double, std::string>> annotations(otherFiles.size(), std::make_tuple(std::numeric_limits<int>::max(), -std::numeric_limits<double>::max(), ""));

        for (size_t i = 0; i < otherFiles.size(); ++i) {

            /** Ensure that all queue elements beyond the 100 bp limit are removed from queue. */
            while (!queues[i].empty()
                && (   (queues[i].front().chrom < mainEntry.chrom)
                    || (queues[i].front().chrom == mainEntry.chrom && queues[i].front().end < mainEntry.start - 100))) {
                queues[i].pop();
            }

            /** Ensure that all queue elements within the 100 bp limit are added to queue. */
            while (true) {
                if (!std::getline(otherFiles[i], line)) break;
                BedEntry otherEntry("", 0, 0, 0.0, '.', line);
                if (!parseBedLine(line, otherEntry)) {
                    std::cerr << "Error parsing line in .bed file: " << line << std::endl;
                    continue;
                }

                if (otherEntry.chrom > mainEntry.chrom || (otherEntry.chrom == mainEntry.chrom && otherEntry.start > mainEntry.end + 100)) {
                    // Unread the line if it's too far downstream or on the next chromosome
                    const std::streampos prevPos = otherFiles[i].tellg();
                    if (!std::getline(otherFiles[i], line)) break;
                    BedEntry otherEntry("", 0, 0, 0.0, '.', line);
                    if (!parseBedLine(line, otherEntry)) {
                        std::cerr << "Error parsing line in .bed file: " << line << std::endl;
                        continue;
                    }

                    if (otherEntry.chrom > mainEntry.chrom || (otherEntry.chrom == mainEntry.chrom && otherEntry.start > mainEntry.end + 100)) {
                        otherFiles[i].seekg(prevPos);
                        break;
                    }

                    if (otherEntry.chrom < mainEntry.chrom || (otherEntry.chrom == mainEntry.chrom && otherEntry.end < mainEntry.start - 100)) {
                        // Discard the entry if it's too far upstream or on the previous chromosome
                        continue;
                    }
                }

                queues[i].push(otherEntry);
            } // while true

            /** Find the best annotation within the 100 bp limit. */
            std::queue<BedEntry> tempQueue = queues[i];
            while (!tempQueue.empty()) {
                const auto& entry = tempQueue.front();
                if (entry.chrom == mainEntry.chrom && std::abs(entry.start - mainEntry.start) <= 100) {
                    if (entry.score > std::get<1>(annotations[i])) {
                        annotations[i] = std::make_tuple(entry.start, entry.score, entry.line);
                    }
                }
                tempQueue.pop();
            }
        } // for each other file

        /** Write to stdout, first the complete line of the main .bed, then the extra columns for the other .bed files */
        std::cout << mainEntry.line;
        for (const auto& annotation : annotations) {
            std::cout << "\t" << std::get<2>(annotation);
        }
        std::cout << std::endl;
    }
}

/** \brief Main entry point.
 * Handling arguments and invoking "processBedFiles"
 */
int main( int argc, char* argv[] ) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <main.bed> <other1.bed> [<other2.bed> ...]" << std::endl;
        return 1;
    }

    std::string mainBedFile = argv[1];
    std::vector<std::string> otherBedFiles;
    for (int i = 2; i < argc; ++i) {
        otherBedFiles.push_back(argv[i]);
    }

    processBedFiles(mainBedFile, otherBedFiles);

    return 0;
}