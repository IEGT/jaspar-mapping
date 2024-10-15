#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <getopt.h>  // For GNU Getopt
#include "progress.h" // progress indicator
#include "gtf_file_region.h"  // GeneRegion struct

static bool beVerbose=false;
static bool showHelp=false;

void printHelp(const std::string& gtfFile, const std::vector<std::string>& geneNames) {
    std::cerr << "Usage: gtf_file_region_retrieval [-g gtf_file] [--gene gene_name1 --gene gene_name2 ...]\n";
    exit(1);
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
