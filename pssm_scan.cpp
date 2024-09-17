#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <algorithm> // For removing brackets
#include <iomanip>   // For formatted output
#include <chrono>    // For progress timing
#include <getopt.h>  // For GNU Getopt

// Function to trim whitespace from strings
std::string trim(const std::string& str) {
    const std::string whitespace = " \t\n\r";
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos) return "";  // No content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

// Function to parse a JASPAR PSSM file format
std::unordered_map<char, std::vector<double>> parsePSSMFile(const std::string& pssmFile, const std::string& targetMotifID, std::string& motifID, std::string& motifName) {
    std::ifstream inFile(pssmFile);
    std::unordered_map<char, std::vector<double>> pssm;

    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << pssmFile << std::endl;
        return pssm;
    }

    std::string line;
    bool readingPSSM = false;
    char currentBase = '\0';

    while (getline(inFile, line)) {
        line = trim(line);

        // Check for the FASTA-like header, indicating a new motif
        if (!line.empty() && line[0] == '>') {
            if ( readingPSSM ) {
                break;
            }
            // Parse motif ID and name
            std::stringstream ss(line.substr(1));  // Skip the '>'
            ss >> motifID >> motifName;
            readingPSSM = (motifID == targetMotifID);  // Only read if the motif ID matches the target ID
            if (readingPSSM) {
                std::cout << "Parsing PSSM for motif: " << motifID << " (" << motifName << ")" << std::endl;
            } else {
                // Do not clear motifID unnecessarily; ensure it is only cleared when invalid.
            }
        } else if (readingPSSM && (line[0] == 'A' || line[0] == 'C' || line[0] == 'G' || line[0] == 'T')) {
            // Read nucleotide counts from the JASPAR format
            std::stringstream ss(line);
            ss >> currentBase;  // The nucleotide character (A, C, G, T)

            std::cerr << "I: currentBase: '" << currentBase << "'" << std::endl;

            std::string countsStr;
            std::getline(ss,countsStr);  // Read the bracketed counts (e.g., [ 4 19 0 0 ])
            std::cerr << "I: countsStr: \"" << countsStr << "\"" << std::endl;

            // Remove brackets and split counts
            countsStr.erase(std::remove(countsStr.begin(), countsStr.end(), '['), countsStr.end());
            countsStr.erase(std::remove(countsStr.begin(), countsStr.end(), ']'), countsStr.end());

            std::stringstream countStream(countsStr);
            double countValue;
            while (countStream >> countValue) {
                pssm[currentBase].push_back(countValue);
            }
        }
    }

    inFile.close();
    return pssm;
}


// Function to calculate log-odds score for a given nucleotide at a position
double logOddsScore(double frequency, double background) {
    if (frequency == 0) {
        return -1e9;  // Prevent log(0) by returning a large negative score for unobserved nucleotides
    }
    return log2(frequency / background);  // Log-odds ratio
}

// Normalize the PSSM by converting counts to log-odds scores
void normalizePSSM(std::unordered_map<char, std::vector<double>>& pssm, const std::unordered_map<char, double>& backgroundFrequencies) {
    for (auto& [nucleotide, counts] : pssm) {
        double background = backgroundFrequencies.at(nucleotide);
        for (double& count : counts) {
            count = logOddsScore(count, background);
        }
    }
}

std::unordered_map<char, double> backgroundFrequencies = {
    {'A', 0.25},  // Assuming equal background probabilities; adjust as needed
    {'C', 0.25},
    {'G', 0.25},
    {'T', 0.25}
};

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

// Function to reverse complement a DNA sequence
std::string reverseComplement(const std::string& sequence) {
    std::string revComp = sequence;
    std::reverse(revComp.begin(), revComp.end());
    for (char& base : revComp) {
        switch (base) {
            case 'A': base = 'T'; break;
            case 'T': base = 'A'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
        }
    }
    return revComp;
}

// Function to slide the PSSM across the DNA sequence and calculate scores
void scanSequence(const std::string& chromosome, const std::string& sequence, const std::string& strand, const std::unordered_map<char, std::vector<double>>& pssm, std::ofstream& outFile, const bool skipN, const float& threshold) {
    size_t motifLength = pssm.begin()->second.size();
    size_t sequenceLength = sequence.size();
    size_t reportInterval = sequenceLength / 100;  // Update progress every 1%

    outFile << "Chromosome\tFrom\tTo\tScore\tStrand\n";

    auto start = std::chrono::high_resolution_clock::now();

    // Slide the window over the sequence
    for (size_t i = 0; i <= sequenceLength - motifLength; ++i) {
        std::string window = sequence.substr(i, motifLength);
        double score = calculateScore(window, pssm, skipN);
        // Skip output if the window contained 'N' or invalid nucleotides
        if (score == -1e9) {
            continue;
        }
        if (score < threshold) {
            continue;
        }
        outFile << chromosome << "\t" << i << "\t" << (i + motifLength) << "\t" << std::fixed << std::setprecision(3) << score << "\t" << strand << "\n";

        // Progress indicator
        if (i % reportInterval == 0) {
            double progress = (double)i / sequenceLength * 100;
            auto now = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
            std::cout << "Progress: " << std::fixed << std::setprecision(2) << progress << "% - Elapsed time: " << elapsed << " seconds\n";
        }
    }
}

// Function to read a genome sequence from a FASTA file (extract only chromosome name)
std::unordered_map<std::string, std::string> readFastaFile(const std::string& fastaFile) {
    std::ifstream inFile(fastaFile);
    std::unordered_map<std::string, std::string> genome;  // Map of chromosome ID to sequence

    if (!inFile.is_open()) {
        std::cerr << "Error opening FASTA file: " << fastaFile << std::endl;
        return genome;
    }

    std::string line, sequence, currentChromosome;
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
        } else {
            // Append the line to the current sequence (FASTA may split sequences across lines)
            sequence += line;
        }
    }

    // Add the last chromosome and sequence
    if (!currentChromosome.empty() && !sequence.empty()) {
        genome[currentChromosome] = sequence;
    }

    inFile.close();
    return genome;
}

void printHelp(const std::string& programName, const std::string& genomeFile, const std::string& pssmFile, const std::string& targetMotifID, const float& threshold) {
    std::cout << "Usage: " << programName << " [-g genome_file] [-p pssm_file] [-m motif_id] [--skip-N | --neutral-N]" << std::endl;
    std::cout << "  -g, --genome     Path to genome FASTA file (set to " << genomeFile  << ")" << std::endl;
    std::cout << "  -p, --pssm       Path to JASPAR PSSM file (set to " << pssmFile << ")" << std::endl;
    std::cout << "  -m, --motif      Target motif ID from JASPAR file (set to " << targetMotifID << ")" << std::endl;
    std::cout << "  -t, --threshold  Minimal score to achieve to print (set to " << threshold << ")" << std::endl;
    std::cout << "  --skip-N         Skip windows containing 'N'" << std::endl;
    std::cout << "  --neutral-N      Treat 'N' as neutral (contribute 0 to the score)" << std::endl;
    std::cout << "  -h, --help       Display this help message" << std::endl;
}

int main(int argc, char* argv[]) {
    // Default file names
    std::string genomeFile = "Homo_sapiens.GRCh38.dna.primary_assembly.fasta";
    std::string pssmFile = "JASPAR2022_CORE_non-redundant_pfms_jaspar.txt";
    std::string targetMotifID;
    float threshold = 0.0f;
    bool skipN = true;  // Default to skipping N
    bool showHelp = false; // Do not show help after all variables have been assigned

    // Option flags and variables for getopt
    int option;
    static struct option long_options[] = {
        {"genome", required_argument, 0, 'g'},
        {"pssm", required_argument, 0, 'p'},
        {"motif", required_argument, 0, 'm'},
        {"threshold", required_argument, 0, 't'},
        {"skip-N", no_argument, 0, 0},
        {"neutral-N", no_argument, 0, 0},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    // Parse command line arguments using getopt
    while ((option = getopt_long(argc, argv, "g:p:m:t:h", long_options, nullptr)) != -1) {
        switch (option) {
            case 'g':
                genomeFile = optarg;
                break;
            case 'p':
                pssmFile = optarg;
                break;
            case 'm':
                targetMotifID = optarg;
                break;
            case 't':
                threshold = static_cast<float>(*optarg);
                std::cerr << "I: Only showing matches with score > " << threshold << "." << std::endl;
            case 0:
                if (std::string(long_options[optind].name) == "skip-N") {
                    skipN = true;
                } else if (std::string(long_options[optind].name) == "neutral-N") {
                    skipN = false;
                }
                break;
            default:
                showHelp = 1;
                break;
        }

        if (showHelp) {
            printHelp(argv[0],genomeFile,pssmFile,targetMotifID,threshold);
            return 1;
        }
    }

    // Ensure the motif ID is provided
    if (targetMotifID.empty()) {
        std::cerr << "Error: You must specify a motif ID with the -m option." << std::endl;
        printHelp(argv[0],genomeFile,pssmFile,targetMotifID,threshold);
        return 1;
    }

    // Load the PSSM matrix from a JASPAR-like file
    std::string motifID, motifName;
    std::unordered_map<char, std::vector<double>> pssm = parsePSSMFile(pssmFile, targetMotifID, motifID, motifName);

    if (motifID.empty()) {
        std::cerr << "E: PSSM with ID '" << targetMotifID << "' not found in file." << std::endl;
        return 1;
    }

    if (pssm.empty()) {
        std::cerr << "E: Failed to parse PSSM." << std::endl;
        return 1;
    }

    // Output filenames based on motif ID and name
    std::string outputFileNamePositive = motifName + "_" + motifID + "_positive.txt";
    std::string outputFileNameNegative = motifName + "_" + motifID + "_negative.txt";

    std::ofstream outFilePositive(outputFileNamePositive);
    std::ofstream outFileNegative(outputFileNameNegative);

    if (!outFilePositive.is_open() || !outFileNegative.is_open()) {
        std::cerr << "E: Error opening output files." << std::endl;
        return 1;
    }

    // Load the genome from a FASTA file
    std::unordered_map<std::string, std::string> genome = readFastaFile(genomeFile);
    if (genome.empty()) {
        std::cerr << "E: Failed to parse genome." << std::endl;
        return 1;
    }

    // Perform PSSM scanning on each chromosome in the genome
    if (pssm.empty() || genome.empty()) {
        std::cerr << "Failed to parse PSSM or genome." << std::endl;
        return 1;
    }

    for (const auto& [chromosome, sequence] : genome) {
        std::cout << "Scanning chromosome " << chromosome << " (+ strand)" << std::endl;

        // Scan positive strand
        scanSequence(chromosome, sequence, "+", pssm, outFilePositive, skipN, threshold);

        // Scan negative strand (reverse complement)
        std::cout << "Scanning chromosome " << chromosome << " (- strand)" << std::endl;
        std::string revCompSequence = reverseComplement(sequence);
        scanSequence(chromosome, revCompSequence, "-", pssm, outFileNegative, skipN, threshold);
    }

    outFilePositive.close();
    outFileNegative.close();

    std::cout << "Results saved to: " << outputFileNamePositive << " and " << outputFileNameNegative << std::endl;

    return 0;
}

