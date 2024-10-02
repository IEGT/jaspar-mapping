#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <filesystem> // test and creating for directories
#include <algorithm> // For removing brackets, replace
#include <iomanip>   // For formatted output / setting precision
#include <chrono>    // For progress timing
#include <getopt.h>  // For GNU Getopt
#include <math.h>
#include "progress.h" // progress indicator

typedef std::unordered_map<std::string, std::string>   genome_type;  //< Map of chromosome ID to sequence
typedef std::unordered_map<char, std::vector<double> > pssm_type;    //< Map of nucleotide to counts

int beVerbose = 0;
int showDebug = 1;

// Function to calculate log-odds score for a given nucleotide at a position
double logOddsScore(const double& frequency, const double& background) {
    if (frequency == 0) {
        //return -1e9;  // Prevent log(0) by returning a large negative score for unobserved nucleotides
        return -2;  // Symmetry
    }
    return std::max(log2(frequency / background),-2.0);  // Log-odds ratio
}

class pssm_class;

/**
 * A class representing a PSSM (Position-Specific Scoring Matrix).
 * Stores nucleotide counts and associated motif information.
 */
class pssm_class {
    public:
        pssm_type pssm;
        std::vector<double> colsums;
        std::string motifID;
        std::string motifName;
        int motifLength;
        /**
         * Default constructor that initializes empty PSSM and motif information.
         */
        pssm_class() {
            pssm.clear();
            motifID.clear();
            motifName.clear();
            motifLength = 0;
            colsums.clear();
        }

        pssm_class(const pssm_class& c) {
            std::cerr << "D: invoked copy constructor" << std::endl;
            this->motifID = c.motifID;
            this->motifName = c.motifName;
            this->motifLength = c.motifLength;
            this->colsums = c.colsums;
            this->pssm = c.pssm;
        }
        /** \brief Constructor that initializes the PSSM and motif information.
         *
         * @param pssm The PSSM matrix
         * @param motifID The motif's ID
         * @param motifName The motif's name
         * @param motifLength The motif's length
         */
        pssm_class(const pssm_type& pssm, const std::string& motifID, const std::string& motifName, const int& motifLength) {
            this->pssm = pssm;
            this->motifID = motifID;
            this->motifName = motifName;
            this->motifLength = motifLength;
            this->colsums.clear();
            for (int i = 0; i< motifLength; i++) {
                this->colsums.push_back(0.0);
            }
            for (auto& [nucleotide, counts] : this->pssm) {
                for (int i = 0; i< motifLength; i++) {
                    this->colsums[i] += counts[i];
                }
            }
            std::cerr << "ColSums: " << std::endl;
            for (int i = 0; i< motifLength; i++) {
                std::cerr << "\t" << this->colsums[i];
            }
            std::cerr << std::endl;
        }
 
        ~pssm_class() {
            std::cerr << "pssm_class - deleting colsums[]" << std::endl;
        }

        static int parsePSSMFile(const std::string& pssmFile, std::unordered_map<std::string, pssm_class>& pssm_list, const std::string& targetMotifID);

        // Normalize the PSSM by converting counts to log-odds scores
        void normalizePSSM(const std::unordered_map<char, const double>& backgroundFrequencies) {

            std::cerr << "I: NormalizePSSM " << std::endl;
            std::cerr << *this;

            for (auto& [nucleotide, counts] : this->pssm) {
                const auto background = backgroundFrequencies.at(nucleotide);

                //int i = 0;
                std::cerr << "Counts.size= " << counts.size() << std::endl;
                //for (unsigned int i=0; auto& count : counts) {
                for (unsigned int i=0; i < counts.size(); i++) {
                    std::cerr << "D: i=" << i << std::endl;
                    if (0 == this->colsums[i]) {
                        //count = 0;
                        counts[i] = 0;
                        std::cerr << "W: Unexpected - column " << i << " has no scores assigned." << std::endl;
                        continue;
                    }
                    if (0.0 == counts[i]) { // an unobserved nucleotide
                        counts[i] = -2;
                        //count = -1e9;
                        continue;
                    }
                    std::cerr << "I: Normalization of " << counts[i] << " counts at position " << i << " with column sum " << this->colsums[i] << " to ";
                    counts[i] = logOddsScore(counts[i]/this->colsums[i], background);
                    std::cerr << counts[i] << std::endl;
                }
            }
        }

        // Overload the << operator for pssm_class
        friend std::ostream& operator<<(std::ostream& os, const pssm_class& pssmObj) {
            os << "Motif ID: " << pssmObj.motifID << std::endl;
            os << "Motif Name: " << pssmObj.motifName << std::endl;
            os << "Motif Length: " << pssmObj.motifLength << std::endl;
            os << "PSSM:" << std::endl;

            // Iterate over the PSSM and print each base and its associated counts
            for (const auto& [base, counts] : pssmObj.pssm) {
                os << base << ": ";
                for (double count : counts) {
                    os << count << " ";
                }
                os << std::endl;
            }

            os << "ColSums:\t" ;
            for (int i = 0; i< pssmObj.motifLength; i++) {
                os << "\t" << pssmObj.colsums[i];
            }
            std::cerr << std::endl;

            return os;
        }
};

typedef std::unordered_map<std::string, pssm_class> pssm_list_type;

// global variable to control verbosity

// Function to trim whitespace from strings
std::string trim(const std::string& str) {
    const std::string whitespace = " \t\n\r";
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos) return "";  // No content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

/** \brief Function to parse a JASPAR PSSM file format
 * The function reads all PSSMs from the file and stores them in a map, that is unless targetMotifID is provided.
 * @param pssmFile - path to JASPAR PSSM file to parse
 * @param pssm - the PSSM structure to read that PSSM file into
 * @param targetMotifID - the motif ID to search for in the PSSM file
 * @returns 0 upon success, else -1.
 */
int pssm_class::parsePSSMFile(const std::string& pssmFile, pssm_list_type& pssm_list, const std::string& targetMotifID) {

    std::ifstream inFile(pssmFile);
    if (!inFile.is_open()) {
        std::cerr << "E: Error opening file: '" << pssmFile << "'" << std::endl;
        return -1;
    }

    std::string line;
    pssm_type pssm;
    bool readingPSSM = false;
    char currentBase = '\0';
    std::string motifID, motifName;

    int length=0;
    while (getline(inFile, line)) {
        line = trim(line);
        // Check for the FASTA-like header, indicating a new motif
        if (!line.empty() && line[0] == '>') {
            // Save the previous PSSM if there is one
            if (   !motifID.empty() // this is not the very first line read,
                                    // so a complete PSSM has already been seen
                && ( targetMotifID.empty() || motifID == targetMotifID )
                                    // no motif was specified or the desired motif was found
                ) {

                pssm_class pssm_object(pssm, motifID, motifName, length);
                pssm_list[motifID] = pssm_object;
                if (showDebug) std::cerr << "I: Got the following PSSM: " << std::endl << pssm_object;
                if (beVerbose) {
                    std::cerr << "I: PSSM for motif " << motifID << " (" << motifName << ") saved with length " << length << "." << std::endl;
                }
                length=0; // redundant, just for clarity
            }

            // We iterate over all provides PSSMs, so clear the current one
            pssm.clear();

            if (!targetMotifID.empty() && motifID == targetMotifID) {
                // Found the target motif, no need to read further
                if (beVerbose) {
                    std::cerr << "I: Target motif " << motifID << " (" << motifName << ") found." << std::endl;
                    std::cerr << "   Parsing PSSM for motif: " << motifID << " (" << motifName << ")" << std::endl;
                }
                break;
            }
            // Parse new motif ID and name from the line read from the file
            std::stringstream ss(line.substr(1));  // Skip the '>'
            ss >> motifID >> motifName;
        } else if (line[0] == 'A' || line[0] == 'C' || line[0] == 'G' || line[0] == 'T'
                // unlikely to see in a PSSM file, but just in case
                || line[0] == 'a' || line[0] == 'c' || line[0] == 'g' || line[0] == 't') {
            // Read nucleotide counts from the JASPAR format
            std::stringstream ss(line);
            ss >> currentBase;  // The nucleotide character (A, C, G, T)

            //if (showDebug) std::cerr << "I: " << currentBase << " : ";
            std::string countsStr;
            std::getline(ss,countsStr);  // Read the bracketed counts (e.g., [ 4 19 0 0 ])
            //if (showDebug) std::cerr << countsStr << std::endl;

            // Remove brackets and split counts
            countsStr.erase(std::remove(countsStr.begin(), countsStr.end(), '['), countsStr.end());
            countsStr.erase(std::remove(countsStr.begin(), countsStr.end(), ']'), countsStr.end());

            std::stringstream countStream(countsStr);
            double countValue;
            length=0;
            while (countStream >> countValue) {
                length++;
                pssm[currentBase].push_back(countValue);
                //if (showDebug) std::cerr << " " << countValue;
                //if (showDebug) std::cerr << " l:" << length;
            }
            //if (showDebug) std::cerr << std::endl;
        } else {
            // Skip other lines
            if (beVerbose) {
                std::cerr << "I: Skipping line: " << line << std::endl;
            }
        }
    }

    inFile.close();
    return 0;
}

std::unordered_map<char, const double> backgroundFrequencies = {
    {'A', 0.25},  // Assuming equal background probabilities; adjust as needed
    {'C', 0.25},
    {'G', 0.25},
    {'T', 0.25},
    {'a', 0.25},
    {'c', 0.25},
    {'g', 0.25},
    {'t', 0.25},
    {'N', 1},
    {'n', 1}
};

// Region constraining the scan for motifs
class Region {
    public:
        std::string chromosome;
        long from;
        long to;
        std::string name;
        
        static std::vector<Region> parseRegionsFile(const std::string& regionsFile);
};

std::vector<Region> Region::parseRegionsFile(const std::string& regionsFile) {
    std::ifstream inFile(regionsFile);
    std::vector<Region> regions;

    if (!inFile.is_open()) {
        std::cerr << "E: Error opening regions file: " << regionsFile << std::endl;
        return regions;
    }

    std::string line;
    bool headerSkipped = false;
    while (getline(inFile, line)) {
        line = trim(line);
        if (line.empty()) continue;

        // Optionally skip a header
        if (!headerSkipped && (line.find("Chromosome") != std::string::npos || line.find("chromosome") != std::string::npos)) {
            headerSkipped = true;
            continue;
        }

        std::stringstream ss(line);
        Region region;

        if (!(ss >> region.chromosome >> region.from >> region.to >> region.name)) {
            std::cerr << "Error parsing line in regions file: " << line << std::endl;
            continue;
        }

        regions.push_back(region);
    }

    inFile.close();
    return regions;
}


// Function to calculate PSSM score for a given window of DNA sequence
double calculateScore(const std::string& window, const std::unordered_map<char, std::vector<double>>& pssm, bool skipN) {
    double score = 0.0;

    for (size_t i = 0; i < window.size(); ++i) {
        char nucleotide = window[i];
        if (nucleotide == 'N') {
            if (skipN) {
                return -1e9;  // Special flag to indicate an invalid window if skipping N
            } else {
                continue;  // Neutral impact on the score if treating N as neutral
            }
        }
        if (pssm.find(nucleotide) != pssm.end()) {
            score += pssm.at(nucleotide)[i];
        } else {
            std::cerr << "Unknown nucleotide: " << nucleotide << std::endl;
            return -1e9;  // Return an extremely low score for invalid nucleotides
        }
    }

    return score;
}

/** \brief Function to reverse complement a DNA sequence
 * Lowercase characters are transposed to upper case in the process
 * @param sequence - the DNA sequence to reverse complement
 * @return the reverse complement of the input sequence
 */
std::string reverseComplement(const std::string& sequence) {
    std::string revComp = sequence;
    std::reverse(revComp.begin(), revComp.end());
    for (char& base : revComp) {
        switch (base) {
            case 'A': base = 'T'; break;
            case 'T': base = 'A'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
            case 'a': base = 'T'; break; // lower case
            case 't': base = 'A'; break;
            case 'c': base = 'G'; break;
            case 'g': base = 'C'; break;
            default: break; // no change required, also N = reverse(N)
        }
    }
    return revComp;
}

/** \brief Function to slide the PSSM across the DNA sequence and calculate scores
 * @param chromosome - the chromosome ID
 * @param sequence - the DNA sequence to scan
 * @param strand - the strand of the sequence
 * @param pssm - the PSSM to use for scoring
 * @param outFile - the output file stream to write results to
 * @param skipN - whether to skip windows containing 'N'
 * @param threshold - the minimal score to achieve to print
 * @param from - the minimal position on the chromosome to consider
 * @param to - the maximal position on the chromosome to consider
 * @param name - the name of the PSSM to be used in the output
 * @param showHeader - whether to show the header in the output file
 * @return 0 if successful, else -1
 */
int scanSequence(const std::string& chromosome, const std::string& sequence, const std::string& strand, const pssm_class& pssm,
                 std::ofstream& outFile, const bool skipN, const float& threshold, const long& from, const long& to, const bool& showHeader, const bool& showSequence) {
    //size_t motifLength = pssm.begin()->second.size();
    size_t motifLength = pssm.motifLength;
    size_t sequenceLength = sequence.size();
    size_t reportInterval = sequenceLength / 100;  // Update progress every 1%

    if (showHeader) {
        outFile << "Chromosome\tFrom\tTo\tName\tScore\tStrand";
        if (showSequence) {
            outFile << "\tMatch";
        }
        outFile << std::endl;
    }

    auto start = std::chrono::high_resolution_clock::now();

    size_t posStart=0L;
    if (from>0L) {
        posStart = (size_t) from;
        std::cerr << "I: Scanning chromosome " << chromosome << " from position " << from << " for motif " << pssm.motifName << std::endl;
    }

    size_t posEnd=sequenceLength;
    if (to>0L && to<sequenceLength) {
        posEnd = (size_t) to;
        std::cerr << "I: Scanning chromosome " << chromosome << " up to position " << to << " for motif " << pssm.motifName << std::endl;
    }

    // Slide the window over the sequence
    if (beVerbose) {
        std::cerr << "I: Scanning chromosome " << chromosome << " for motif " << pssm.motifName << " from " << posStart << " to " << posEnd << " - " << motifLength << std::endl;
    }
    for (size_t i = posStart; i <= posEnd - motifLength; ++i) {
        const std::string window = sequence.substr(i, motifLength);
        const double score = calculateScore(window, pssm.pssm, skipN);
        // Skip output if the window contained 'N' or invalid nucleotides
        std::cerr << "D: Score: " << score << std::endl;
        if (score == -1e9) {
            continue;
        }
        if (score < threshold) {
            continue;
        }
        outFile << chromosome << "\t" << (i+1) << "\t" << (i+1 + motifLength) << "\t" << pssm.motifName << "\t" << std::fixed << std::setprecision(3) << score << "\t" << strand;
        if (showSequence) {
            outFile << "\t" << window;
        }
        outFile << std::endl;

        // Progress indicator
        if (i % reportInterval == 0) {
            double progress = (double)i / sequenceLength * 100;
            auto now = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
            std::cout << "Progress: " << std::fixed << std::setprecision(2) << progress << "% - Elapsed time: " << elapsed << " seconds\n";
        }
    }
    return 0;
}

/** \brief Function to read a genome sequence from a FASTA file (extract only chromosome name)
 * No check is performed if the genome representation passed is empty or already filled with some chromosomes.
 * Chromosomes added with the same name would substitute the entry passed.
 * @param fastaFile - path to FASTA genome representation to read
 * @genome - the genome structure to read that genome file into
 * @return -1 on error, else 0.
 */
int readFastaFile(const std::string& fastaFile, genome_type& genome) {

    if (beVerbose) std::cerr << "I: Reading genome from FASTA file: " << fastaFile << std::endl;

    if (! genome.empty() ) {
        std::cerr << "W: readFastaFile: genome passed is not empty." << std::endl;
    }
    std::ifstream inFile(fastaFile);

    if (!inFile.is_open()) {
        std::cerr << "E: Error opening FASTA file: " << fastaFile << std::endl;
        return -1;
    }
    inFile.seekg (0, inFile.end);
    std::streampos genomeFileSize = inFile.tellg();
    inFile.seekg (0, inFile.beg);
    if (genomeFileSize <= 0) {
        std::cerr << "E: Invalid file size for genome: " << genomeFileSize << std::endl;
        return -1;
    }


    std::string line, sequence, currentChromosome;
    std::streamsize bytesRead = 0;
    size_t lineNo=0;
    while (getline(inFile, line)) {
        line = trim(line);

        if (line.empty()) continue;

        // Check for the header
        if (line[0] == '>') {
            // Save the previous sequence if there is one
            if (!currentChromosome.empty() && !sequence.empty()) {
                genome[currentChromosome] = sequence;
                sequence.clear();
            }

            // Get the chromosome ID (everything after the '>' and before the first space)
            currentChromosome = line.substr(1, line.find(' ') - 1);

            if (!genome[currentChromosome].empty()) {
                std::cerr << "W: readFastaFile: Overwriting current entry for chromosmoe '" << currentChromosome << "'" << std::endl;
            }
        } else {
            // Append the line to the current sequence (FASTA may split sequences across lines)
            sequence += line;
        }

        // Update the number of bytes read
        bytesRead += line.size() + 1;  // +1 for the newline character
        // Calculate the progress as a float between 0 and 1
        // Display the progress bar
        if (lineNo++ % 10000 == 0) {
            float progress = static_cast<float>(bytesRead) / genomeFileSize;
            displayProgressBar(progress);
        }

    }

    // Add the last chromosome and sequence
    if (!currentChromosome.empty() && !sequence.empty()) {
        genome[currentChromosome] = sequence;
    }

    inFile.close();
    std::cout << std::endl;
    return 0;
}

void printHelp(const std::string& programName, const std::string& genomeFile, const std::string& pssmFile, const std::string& targetMotifID, const float& threshold, const std::string& chromosome, const size_t& from, const size_t& to, const std::string& regions, const std::string& outdir, const bool& showSequence) {
    std::cout << "Usage: " << programName << " [-v] [-c chromosome] [-t toBp] [-f fromBp] [-g genome_file] [-p pssm_file] [-m motif_id] [--skip-N | --neutral-N] [--skip-normalization]" << std::endl;
    std::cout << " -v, --verbose        Allow verbose output (set to " << beVerbose << ")" << std::endl;
    std::cout << " -g, --genome         Path to genome FASTA file (set to '" << genomeFile  << "')" << std::endl;
    std::cout << " -p, --pssm           Path to JASPAR PSSM file (set to '" << pssmFile << "')" << std::endl;
    std::cout << " -m, --motif          Target motif ID from JASPAR file (set to '" << targetMotifID << "')" << std::endl;
    std::cout << " -l, --threshold      Minimal score to achieve to print (set to " << threshold << ")" << std::endl;
    std::cout << " -c, --chr            Single chromosome to consider (set to '" << chromosome << "')" << std::endl;
    std::cout << " -f, --from           Minimal position on chromosome to consider (set to " << from << ")" << std::endl;
    std::cout << " -t, --to             Maximal position on chromosome to consider (set to " << to << ")" << std::endl;
    std::cout << " -r, --regions        Regions constraint file (set to '" << regions << "')" << std::endl;
    std::cout << " -o, --outdir         Directory to create output files in (set to '" << outdir << "')" << std::endl;
    std::cout << " -s, --show-sequence  Adds the sequence matched by the motif as another field in the output (set to " << showSequence << ")" << std::endl;
    std::cout << " --skip-N             Skip windows containing 'N'" << std::endl;
    std::cout << " --neutral-N          Treat 'N' as neutral (contribute 0 to the score)" << std::endl;
    std::cout << " -N, --skip-normalization Skip log-normalisation, will affect scoring." << std::endl;
    std::cout << " -h, --help           Display this help message" << std::endl;
}

int main(int argc, char* argv[]) {
    // Default file names
    std::string genomeFile = "Homo_sapiens.GRCh38.dna.primary_assembly.fasta";
    std::string pssmFile = "JASPAR2022_CORE_non-redundant_pfms_jaspar.txt";
    std::string regionsFile = "";
    std::string targetMotifID;
    float threshold = - 1e9; // Very small number, i.e. effectively no threshold
    bool showSequence = false;
    bool skipN = true;  // Default to skipping N
    bool skipNormalization = false;  // Default to not skip normalization
    bool showHelp = false; // Do not show help after all variables have been assigned
    long targetFrom = -1L, targetTo = -1L ;
    std::string targetChromosome ;
    std::vector<Region> regions;
    std::string outdir = ".";

    // Option flags and variables for getopt
    int option;
    static struct option long_options[] = {
        {"genome", required_argument, 0, 'g'},
        {"pssm", required_argument, 0, 'p'},
        {"motif", required_argument, 0, 'm'},
        {"threshold", required_argument, 0, 'l'},
        {"chr", required_argument, 0, 'c'},
        {"from", required_argument, 0, 'f'},
        {"to", required_argument, 0, 't'},
        {"regions", required_argument, 0, 'r'},  // regions file
        {"outdir", required_argument, 0, 'o'},  // output directory
        {"show-sequence", no_argument, 0, 's'},  // output directory
        {"skip-N", no_argument, 0, 0},
        {"neutral-N", no_argument, 0, 0},
        {"skip-normalization", no_argument, 0, 'N'},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    // Parse command line arguments using getopt
    int option_index = 0;
    while ((option = getopt_long(argc, argv, "g:p:m:l:c:f:t:o:Nsvh", long_options, &option_index)) != -1) {
        switch (option) {
            case 'g':
                genomeFile = optarg;
                break;
            case 'p':
                pssmFile = optarg;
                break;
            case 'm':
                targetMotifID = optarg;
                std::cerr << "I: Setting motif ID to " << targetMotifID << std::endl;
                break;
            case 'r':
                regionsFile = optarg;
                break;
            case 'c':
                targetChromosome = optarg;
                std::cerr << "I: Only showing matches on chromosome " << targetChromosome << "." << std::endl;
                break;
            case 'f':
                targetFrom = atol(optarg);
                std::cerr << "I: Only showing matches with position downstream of " << targetFrom << "." << std::endl;
                break;
            case 't':
                targetTo = atol(optarg);
                std::cerr << "I: Only showing matches with position up to " << targetTo << "." << std::endl;
                break;
            case 'l':
                threshold = atof(optarg);
                std::cerr << "I: Only showing matches with score > " << threshold << "." << std::endl;
                break;
            case 'o':
                outdir = optarg;
                if (!std::filesystem::is_directory(outdir)) {
                    if (std::filesystem::create_directory(outdir)) {
                        std::cerr << "I: Created output directory '" << outdir << "'." << std::endl;
                    } else {
                        std::cerr << "E: Output directory '" << outdir << "' not existing and could not be created." << std::endl;
                    }
                }
                break;
            case 0:
                if (std::string(long_options[option_index].name) == "skip-N") {
                    skipN = true;
                } else if (std::string(long_options[option_index].name) == "neutral-N") {
                    skipN = false;
                }
                break;
            case 's':
                showSequence=true;
                break;
            case 'N':
                skipNormalization=true;
                break;
            case 'v':
                beVerbose = 1;
                break;
            default:
                showHelp = 1;
                break;
        }

        if (showHelp) {
            printHelp(argv[0],genomeFile,pssmFile,targetMotifID,threshold,targetChromosome,targetFrom,targetTo,regionsFile,outdir,showSequence);
            return 1;
        }
    }

/*
    // Ensure the motif ID is provided
    if (targetMotifID.empty()) {
        std::cerr << "E: You must specify a motif ID with the -m option." << std::endl;
        printHelp(argv[0],genomeFile,pssmFile,targetMotifID,threshold,targetChromosome,targetFrom,targetTo,regionsFile,outdir,showSequence);
        return 1;
    }
*/

    // Load the PSSM matrix from a JASPAR-like file
    pssm_list_type pssm_list;

    if ( pssm_class::parsePSSMFile(pssmFile, pssm_list , targetMotifID ) ) {
        std::cerr << "E: Error parsing PSSM file '" << pssmFile << "'" << std::endl;
        return 1;
    } else if (pssm_list.empty()) {
        std::cerr << "E: PSSM sucessfully parsed but nonetheless empty." << std::endl;
        return 1;
    }
    std::cerr << "I: Read " << pssm_list.size() << " PSSMs from file '" << pssmFile << "'" << std::endl;

    // Load the genome from a FASTA file
    genome_type genome ;
    if ( 0 != readFastaFile(genomeFile, genome) ) {
        std::cerr << "E: Error in function reading the genome." << std::endl;
        return 1;
    } else if (genome.empty()) {
        std::cerr << "E: Genome apparently successfully parsed but nonetheless empty." << std::endl;
        return 1;
    }


    // iterating over PSSMs in list
    for(auto const& [motifID, pssm_object] : pssm_list) {

        // ensure we can modify the object, e.g. for normalization
        pssm_class pssm_object_copy = pssm_object;

        // Perform PSSM scanning on each chromosome in the genome
        if (pssm_object_copy.pssm.empty() ) {
            std::cerr << "E: Failed to parse PSSM." << std::endl;
            return 1;
        }

        pssm_type pssm = pssm_object_copy.pssm;
        std::cerr << "I: Found motif length to be " << pssm.begin()->second.size() << " but maybe better is " << pssm_object.motifLength << std::endl;

        if ( motifID != pssm_object.motifID ) {
            std::cerr << "E: Mismatch in motif ID." << std::endl;
            return 1;
        }
        std::string motifName = pssm_object_copy.motifName;

        if (skipNormalization) {
            std::cerr << "I: Skipping normalization." << std::endl;
        } else {
            pssm_object_copy.normalizePSSM(backgroundFrequencies);
        }

        std::cerr << "PSSM after normalization:" << std::endl << pssm_object_copy;

        if (!regionsFile.empty()) {
            regions = Region::parseRegionsFile(regionsFile);
            if (regions.empty()) {
                std::cerr << "E: No valid regions found in regions file." << std::endl;
                return 1;
            }
        }

        if (!regions.empty() && (!targetChromosome.empty() || targetFrom > 0L || targetTo > 0L)) {
            std::cerr << "E: Cannot specify both individual chromosome options and a regions file." << std::endl;
            return 1;
        }

        // Output filenames based on motif ID and name
        //
        std::string motifNameForFile = motifName;
        std::replace(motifNameForFile.begin(), motifNameForFile.end(), ':', '_'); // Replace colon with underscore
        std::replace(motifNameForFile.begin(), motifNameForFile.end(), '/', '_'); // Replace colon with underscore
        std::string outputFileNamePositive = outdir + std::filesystem::path::preferred_separator + motifNameForFile + "_" + motifID + "_positive";
        std::string outputFileNameNegative = outdir + std::filesystem::path::preferred_separator + motifNameForFile + "_" + motifID + "_negative";
        if (!targetChromosome.empty()) {
            outputFileNamePositive += "_"+targetChromosome;
            outputFileNameNegative += "_"+targetChromosome;
        }
        if (targetFrom > 0L) {
            outputFileNamePositive += "_"+std::to_string(targetFrom);
            outputFileNameNegative += "_"+std::to_string(targetFrom);
        }
        if (targetTo > 0L) {
            outputFileNamePositive += "-"+std::to_string(targetTo);
            outputFileNameNegative += "-"+std::to_string(targetTo);
        }
        if (threshold) {
            outputFileNamePositive += "_thresh_";
            outputFileNamePositive += std::to_string(threshold);
            outputFileNameNegative += "_thresh_";
            outputFileNameNegative += std::to_string(threshold);
        }
        outputFileNamePositive += ".bed";
        outputFileNameNegative += ".bed";

        std::ofstream outFilePositive(outputFileNamePositive);
        std::ofstream outFileNegative(outputFileNameNegative);

        if (!outFilePositive.is_open()) {
            std::cerr << "E: Error opening output file '" << outputFileNamePositive << "'." << std::endl;
            return 1;
        }

        if (!outFileNegative.is_open()) {
            std::cerr << "E: Error opening output file '" << outputFileNameNegative << "'." << std::endl;
            return 1;
        }

        bool showHeader=true;

        if (!regions.empty()) {
            // Scan specified regions
            for (const auto& region : regions) {
                auto it = genome.find(region.chromosome);
                if (it == genome.end()) {
                    std::cerr << "Chromosome " << region.chromosome << " not found in genome." << std::endl;
                    continue;
                }

                const std::string& sequence = it->second;
                size_t seqLength = sequence.length();
                size_t regionFrom = static_cast<size_t>(std::max(0L, region.from));
                size_t regionTo = static_cast<size_t>(std::min(static_cast<long>(seqLength), region.to));

                if (regionFrom >= regionTo) {
                    std::cerr << "E: Invalid region: " << region.chromosome << " " << region.from << "-" << region.to << std::endl;
                    continue;
                }

                std::cerr << "I: Scanning region: " << region.chromosome << ":" << regionFrom << "-" << regionTo << " (" << region.name << ") (+ strand)" << std::endl;

                // Scan positive strand
                scanSequence(region.chromosome, sequence, "+", pssm_object_copy, outFilePositive, skipN, threshold, regionFrom, regionTo, showHeader, showSequence);

                // Scan negative strand
                std::cerr << "Scanning region: " << region.chromosome << ":" << regionFrom << "-" << regionTo << " (" << region.name << ") (- strand)" << std::endl;
                std::string revCompSequence = reverseComplement(sequence);
                scanSequence(region.chromosome, revCompSequence, "-", pssm_object_copy, outFileNegative, skipN, threshold, regionFrom, regionTo, showHeader, showSequence);

                showHeader = false;
            }

        } else {
            // Original scanning logic for full genome or chromosome ranges
            for (const auto& [chromosome, sequence] : genome) {
                bool chromosomeFound = false;
                if ( chromosome.empty() ) {
                    std::cerr << "E: Stored data on empty chromosome - weird, please check." << std::endl;
                    outFilePositive.close();
                    outFileNegative.close();
                    return 1;
                }
                if ( !targetChromosome.empty() ) {
                    if ( 0 != chromosome.compare(targetChromosome) ) {
                        if (showDebug) {
                            std::cerr << "I: Skipping chromosome " << chromosome << " (+ strand)" << std::endl;
                        }
                        continue;
                    } else {
                        chromosomeFound = true;
                    }
                }

                std::cerr << "I: Scanning chromosome " << chromosome << " from " << targetFrom << " to " << targetTo << " (+ strand)" << std::endl;

                // Scan positive strand
                scanSequence(chromosome, sequence, "+", pssm_object_copy, outFilePositive, skipN, threshold, targetFrom, targetTo, showHeader, showSequence);

                // Scan negative strand (reverse complement)
                std::cout << "I: Scanning chromosome " << chromosome << " from " << targetFrom << " to " << targetTo << " (- strand)" << std::endl;
                std::string revCompSequence = reverseComplement(sequence);
                scanSequence(chromosome, revCompSequence, "-", pssm_object_copy, outFileNegative, skipN, threshold, targetFrom, targetTo, showHeader, showSequence);

                if (chromosomeFound) {
                    // Found single chromosome of interest, completed its search - loop shall end here.
                    break;
                }
                showHeader=false;
            }
        }

        outFilePositive.close();
        outFileNegative.close();

        std::cout << "Results saved to: " << outputFileNamePositive << " and " << outputFileNameNegative << std::endl;
    }


    return 0;
}

