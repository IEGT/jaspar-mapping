#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <getopt.h>  // For GNU Getopt
#include "progress.h" // progress indicator

// Structure to store gene region information
struct GeneRegion {
    std::string chromosome;
    size_t start;
    size_t end;
    std::string id;
};

bool beVerbose=false;
bool showHelp=false;

void storeGeneReference(std::unordered_map<std::string, GeneRegion>& geneRegions, const std::string geneIdentifier, std::string chrom, size_t start, size_t end) {

        // Check if the gene already exists in geneRegions
        if (geneRegions.find(geneIdentifier) != geneRegions.end()) {
           //std::cerr << "D: Gene " << geneIdentifier << " already exists in geneRegions." << std::endl;
           GeneRegion existing = geneRegions[geneIdentifier];
           // Bail out if two genes with the same name are found to exits on different chromosomes.
           if (existing.chromosome != chrom) {
               std::cerr << "E: inconsistency for gene " << geneIdentifier << " - existing.chromosome = " << existing.chromosome << " != " << chrom << " = chrom" << std::endl;
               return;
           }
           // Ignore the shorter version
           if (abs(existing.end - existing.start) > abs(end-start)) {
               //std::cerr << "I: Ignoring alternative shorter annotation for gene " << geneIdentifier << std::endl;
               return;
           }
        }

        // Store the gene region
        geneRegions[geneIdentifier] = {chrom, start, end, geneIdentifier};
}


// Function to parse the GTF file and extract regions for specific gene names
std::unordered_map<std::string, GeneRegion> parseGTFFile(const std::string& gtfFile) {
    std::unordered_map<std::string, GeneRegion> geneRegions;

    std::ifstream inFile(gtfFile, std::ifstream::ate | std::ifstream::binary);
    if (!inFile.is_open()) {
        std::cerr << "Error opening GTF file: " << gtfFile << std::endl;
        return geneRegions;
    }
    std::streampos gtfFileSize = inFile.tellg();
    std::cerr << "gtfFileSize == " << gtfFileSize << std::endl;

    if (gtfFileSize <= 0) {
        std::cerr << "Error: File size is 0 or file is empty: " << gtfFile << std::endl;
        return geneRegions;
    }

    inFile.seekg(0);

    std::string line;
    size_t bytesRead = 0L;
    size_t lineNo = 0L;
    while (getline(inFile, line)) {
        lineNo++;
        std::istringstream ss(line);
        bytesRead += line.size();

                     //1   havana	gene    2581560   2584533  .      +       .      gene_id "ENSG00000228037"; gene_version "1"; gene_source "havana"; gene_biotype "lncRNA";
        std::string chrom, source, feature, startStr, endStr, score, strand, frame, attributes;
        if (!(ss >> chrom >> source >> feature >> startStr >> endStr >> score >> strand >> frame)) {
            continue;  // Skip lines that don't match the expected format
        }
        if (chrom.starts_with('#')) {
            continue;  // Skip header / comment
        }
        std::getline(ss,attributes);

        if (feature != "gene"  && feature != "transcript") continue;  // We only care about "gene" features

        long start = std::stol(startStr);
        long end = std::stol(endStr);

        std::string geneName, geneID;

        // Extract the gene name from the attributes field
        size_t geneNamePos = attributes.find("gene_name");
        if (geneNamePos != std::string::npos) {
            size_t startQuote = attributes.find("\"", geneNamePos);
            if (std::string::npos != startQuote) {
                size_t endQuote = attributes.find("\"", startQuote + 1);
                geneName = attributes.substr(startQuote + 1, endQuote - startQuote - 1);
                if ("TP73" == geneName) {
                    std::cerr << "I: found gene name '" << geneName << "' at line " << lineNo << " from " << startQuote << " to " << endQuote << "." << std::endl;
                }
            } else {
                //std::cerr << "D: endQuote: " << endQuote << std::endl;
                geneName.clear();
            }
        }

        size_t geneIDPos = attributes.find("gene_id");
        if (geneIDPos != std::string::npos) {
            size_t startQuote = attributes.find("\"", geneIDPos);
            if (std::string::npos == startQuote) {
                std::cerr << "E: Line " << lineNo << ": Found entry 'gene_id' at position " << geneIDPos << " of string '" << attributes << "' but did not find subsequent double quote." << std::endl;
                exit(1);
            }
            size_t endQuote = attributes.find("\"", startQuote + 1);
            if (std::string::npos == endQuote) {
                std::cerr << "E: Line " << lineNo << ": Found entry 'gene_id' at position " << geneIDPos << " of string '" << attributes << "' but did not second subsequent double quote." << std::endl;
                exit(1);
            }
            if (endQuote > startQuote) {
                geneID = attributes.substr(startQuote + 1, endQuote - startQuote - 1);
                //std::cerr << "D: found gene ID '" << geneID << "' at line " << lineNo << " from " << startQuote << " to " << endQuote << "." << std::endl;
            } else {
                //std::cerr << "D: endQuote: " << endQuote << std::endl;
                geneID.clear();
            }
        }

        if (!geneName.empty()) storeGeneReference(geneRegions, geneName, chrom, start, end);
        if (!geneID.empty()) storeGeneReference(geneRegions, geneID, chrom, start, end);

        if (lineNo % 100 == 0) {
            float progress = static_cast<float>(bytesRead) / gtfFileSize;
            displayProgressBar(progress);
        }

    }
    inFile.close();
    return geneRegions;
}

void printHelp(const std::string& gtfFile, const std::vector<std::string>& geneNames) {
    std::cerr << "Usage: gtf_file_region_retrieval [-g gtf_file] [--gene gene_name1 --gene gene_name2 ...]\n";
    exit(1);
}

// Function to read gene names from stdin
std::vector<std::string> readGenesFromStdin() {
    std::vector<std::string> genes;
    std::string gene;
    std::cout << "Enter gene names (one per line). Press Ctrl+D to finish:\n";
    while (std::getline(std::cin, gene)) {
        if (!gene.empty()) {
            genes.push_back(gene);
        }
    }
    return genes;
}

int main(int argc, char* argv[]) {
    std::string gtfFile = "Homo_sapiens.GRCh38.112.gtf";  // Default GTF file
    std::vector<std::string> geneNames;

    // Option flags and variables for getopt
    int option;
    static struct option long_options[] = {
        {"gtf", required_argument, 0, 'g'},   // Path to GTF file
        {"gene", required_argument, 0, 'n'},  // Specify genes of interest
        {"verbose", required_argument, 0, 'v'},  // Instructs to generate help info
        {"help", no_argument, 0, 'h'}, // Instructs to be chatty about what is happening
        {0, 0, 0, 0}
    };

    // Parse command line arguments using getopt
    while ((option = getopt_long(argc, argv, "g:n:vh", long_options, nullptr)) != -1) {
        switch (option) {
            case 'g':
                gtfFile = optarg;
                break;
            case 'n':
                geneNames.push_back(optarg);  // Collect the gene names
                break;
            case 'h':
                showHelp=true;
                break;
            case 'v':
                beVerbose=true;
            default:
                return 1;
        }
    }

    if (showHelp) {
        printHelp(gtfFile,geneNames);
    }

    // If no gene names were provided, read from stdin
    if (geneNames.empty()) {
        geneNames = readGenesFromStdin();
    }

    if (geneNames.empty()) {
        std::cerr << "Error: No genes provided." << std::endl;
        return 1;
    }

    // Parse the GTF file to get regions
    std::unordered_map<std::string, GeneRegion> geneRegions = parseGTFFile(gtfFile);

    // Print the regions for the genes
    for (const std::string& geneName : geneNames) {
        std::cerr << "D: testing for relevance of gene name '" << geneName << "'" << std::endl;
        if (geneRegions.find(geneName) != geneRegions.end()) {
            const GeneRegion& region = geneRegions[geneName];
            std::cout << "Gene " << geneName << " is located on chromosome "
                      << region.chromosome << " from " << region.start << " to " << region.end << std::endl;
        } else {
            std::cerr << "Gene " << geneName << " not found in the GTF file." << std::endl;
            std::cerr << "Got the following to choose from:" << std::endl;
            for (const auto & [geneName,geneRegion] : geneRegions) {
                //std::cerr << geneName << std::endl;
            }
        }
    }



    return 0;
}
