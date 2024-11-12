#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <tuple>
#include <limits>
#include <cmath>

// Function to get the basename of a file path
std::string basename(const std::string& filePath) {
    const size_t lastSlash = filePath.find_last_of("/\\");
    if (lastSlash == std::string::npos) {
        return filePath; // No slash found, return the whole path
    }
    return filePath.substr(lastSlash + 1);
}

/** \brief Representation of single line in .bed file */
struct BedEntry {
    bool hasName=false;
    std::string chrom;
    int start;
    int end;
    std::string name;
    double score;
    std::string line;
    char strand;

    BedEntry(const std::string& chrom, int start, int end,
            const std::string& name, const double& score,
            const char strand = '.', const std::string& line = "", const bool hasName=false)
        : chrom(chrom), start(start), end(end), name(name), score(score), strand(strand), line(line), hasName(hasName) {}
};

/** \brief Splits line of .bed file and retrieves values of relevance
 * @param line - line from .bed file
 * @param entry - BedEntry object to store values
 */
bool parseBedLine(const std::string& line, BedEntry& entry, const bool hasName=false) {
    std::istringstream iss(line);
    std::string chrom;
    int start, end;
    std::string name;
    double score;
    char strand;

    if ( hasName && iss >> chrom >> start >> end >> name >> score >> strand) {
        entry = BedEntry(chrom, start, end, name, score, strand, line, true);
        return true;
    } else if (iss >> chrom >> start >> end >> score >> strand) {
        entry = BedEntry(chrom, start, end, "", score, strand, line, false);
        return true;
    }
    std::cerr << "Error parsing line (hasName="<<hasName<<"): " << line << std::endl;
    return false;
}

/** \brief Process main .bed file and list of other .bed files
 * @param mainBedFile - path to main .bed file
 * @param otherBedFiles - list of paths to other .bed files
 */
void processBedFiles(const std::string& mainBedFile,
                     const std::vector<std::string>& otherBedFiles, const bool hasName=false) {
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
    std::vector<std::string> shortName(otherFiles.size());
    for (size_t i = 0; i < shortName.size(); ++i) {
        shortName[i] = basename(otherBedFiles[i]);
    }

    /** Iterate over all lines in the main .bed file and annotate it. */
    std::string line;
    size_t lineCount = 0;
    while (std::getline(mainFile, line)) {

        lineCount++;
        bool readingHeader=0;

        /** Read next entry of main .bed file. */
        BedEntry mainEntry("", 0, 0, "", 0.0, '.', line, true);

        if (line.starts_with('#')) continue;
        if (line.starts_with("Chr")) {
            std::cout << line;
            readingHeader=1;
        }

        if (!readingHeader && !parseBedLine(line, mainEntry, true)) {
            std::cerr << "Main file: Error parsing line " << lineCount << " in main .bed file:"
                      << std::endl << line << std::endl;
            continue;
        }

        /** Store the best annotation for each other file, see header of other files for details. */
        std::vector<std::tuple<int, double, bool>> annotations(otherFiles.size(),
                        std::make_tuple(std::numeric_limits<int>::max(), -std::numeric_limits<double>::max(), 0));
        std::vector<std::size_t> otherLineCount(otherFiles.size(), 0);

        std::cerr << "Processing line " << lineCount << " : ";

        for (size_t i = 0; i < otherFiles.size(); ++i) {
            std::cerr << ".";

            /** Ensure that all queue elements beyond the 100 bp limit are removed from queue. */
            while (!queues[i].empty()
                && (   (queues[i].front().chrom < mainEntry.chrom)
                    || (queues[i].front().chrom == mainEntry.chrom && queues[i].front().end < mainEntry.start - 100))) {
                queues[i].pop();
            }

            /** Ensure that all queue elements within the 100 bp limit are added to queue. */
            while (true) {
                if (!std::getline(otherFiles[i], line)) break;
                otherLineCount[i]++;

                if (line.starts_with('#')) continue;
                if (line.starts_with("Chr")) {
                    if (!readingHeader) {
                        std::cerr << "Other file '" << otherBedFiles[i] << "': reading header when expecting a regular line." << std::endl;
                        abort();
                    }
                    std::cout << "\t" << shortName[i]<<".Shift" << "\t" << shortName[i]<<".Score" << shortName[i]<<".Strand.Equal";
                    continue;
                }

                BedEntry otherEntry("", 0, 0, "", 0.0, '.', line, true);
                if (!parseBedLine(line, otherEntry, true)) {
                    std::cerr << "Other file '" << otherBedFiles[i] << "': Error parsing line '" << otherLineCount[i] << "' in .bed file: " << line << std::endl;
                    continue;
                }

                if (    otherEntry.chrom > mainEntry.chrom
                    || (otherEntry.chrom == mainEntry.chrom && otherEntry.start > mainEntry.end + 100)) {
                    // Unread the line if it's too far downstream or on the next chromosome
                    const std::streampos prevPos = otherFiles[i].tellg();
                    if (!std::getline(otherFiles[i], line)) break;
                    BedEntry otherEntry("", 0, 0, "", 0.0, '.', line, true);
                    if (!parseBedLine(line, otherEntry, true)) {
                        std::cerr << "Error parsing line in .bed file: " << line << std::endl;
                        continue;
                    }

                    if (    otherEntry.chrom > mainEntry.chrom
                        || (otherEntry.chrom == mainEntry.chrom && otherEntry.start > mainEntry.end + 100)) {
                        otherFiles[i].seekg(prevPos);
                        break;
                    }

                    if (    otherEntry.chrom < mainEntry.chrom
                        || (otherEntry.chrom == mainEntry.chrom && otherEntry.end < mainEntry.start - 100)) {
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
                        annotations[i] = std::make_tuple(entry.start - mainEntry.start, entry.score, entry.strand == mainEntry.strand );
                    }
                }
                tempQueue.pop();
            }
        } // for each other file
        std::cerr << std::endl; // end of line processing

        if (readingHeader) {
            std::cout << std::endl;
            continue;
        }

        /** Write to stdout, first the complete line of the main .bed, then the extra columns for the other .bed files */
        std::cout << mainEntry.line;
        for (const auto& annotation : annotations) {
            std::cout << "\t" << std::get<0>(annotation) << "\t" << std::get<1>(annotation) << "\t" << std::get<2>(annotation);
        }
        std::cout << std::endl;
    }
}

/** \brief Main entry point.
 * Handling arguments and invoking "processBedFiles"
 */
int main( int argc, const char* const argv[] ) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <main.bed> <other1.bed> [<other2.bed> ...]" << std::endl;
        std::cerr << argv[0] << " - Annotate a .bed file with annotations from other .bed files." << std::endl;
        std::cerr << std::endl;
        std::cerr << "The program reads the main .bed file and annotates each line with the highest-scoring annotation from the other .bed files within a 100 bp window. "
                     "The output is written to stdout." << std::endl;

        std::cerr << "The .bed files must be sorted by chromosome and start position, "
                     "chromosome names be withot a 'chr' prefix." << std::endl;

        std::cerr << "To start from a UNIX shell, you may use the following command:" << std::endl
                  << argv[0] << " main.bed somefolder/*positive.bed" << std::endl;

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
