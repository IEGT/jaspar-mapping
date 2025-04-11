/*
 * This file, gtf_file_region_retrieval.cpp, is a command-line utility designed to extract specific gene regions from a GTF (Gene Transfer Format) file.
 * 
 * GTF files are used in bioinformatics to store information about gene structure, including the locations of exons, coding sequences, and other genomic features.
 * 
 * Key Components of the Code:
 * 
 * 1. Command-Line Argument Parsing:
 *    - The program uses `getopt_long` to parse command-line arguments. It supports options for specifying the GTF file (`-g` or `--gtf`), gene names (`-n` or `--gene`), verbosity (`-v` or `--verbose`), and help (`-h` or `--help`).
 * 
 * 2. Default GTF File:
 *    - If no GTF file is specified, it defaults to `Homo_sapiens.GRCh38.112.gtf`.
 * 
 * 3. Help Function:
 *    - The `printHelp` function provides usage instructions and exits the program.
 * 
 * 4. Gene Name Input:
 *    - Gene names can be provided via command-line arguments. If no gene names are provided, the program reads them from standard input.
 *    - Multiple gene names can be passed together, separated by a comma, a blank, or a combination of the two.
 *    The gene names are stored in a vector for to serve as a filter for the output.
 * 
 * 5. GTF File Parsing:
 *    - The `parseGTFFile` function (assumed to be defined elsewhere) reads the GTF file and extracts gene regions into an `unordered_map` where the key is the gene name and the value is a `GeneRegion` struct.
 * 
 * 6. Output:
 *    - For each gene name provided, the program checks if it exists in the GTF file and prints its location (chromosome, start, and end positions). If a gene is not found, it lists the available gene names from the GTF file.
 *    - Output is filtered for gene names selected.
 *    - Coordinates can be constrained on promoter regions from 1 to 500 bp upstream (considering the strand information) of the transcription start site.
 * 
 * Purpose of the Tool:
 * 
 * This tool is essential for bioinformatics workflows where specific gene regions need to be extracted from a GTF file for further analysis. For example, these regions can be intersected with `.bed` files, which are commonly used to represent genomic intervals in various analyses, such as identifying overlaps with other genomic features or regions of interest.
 * 
 * Example Usage:
 * 
 * ./gtf_file_region_retrieval -g path/to/gtf_file.gtf --gene BRCA1 --gene TP53
 * ./gtf_file_region_retrieval -g path/to/gtf_file.gtf --gene "BRCA1, TP53"
 * 
 * This command would extract the regions for the genes `BRCA1` and `TP53` from the specified GTF file and print their locations.
 * 
 * Why This Tool is Needed:
 * 
 * - Efficiency: Automates the extraction of gene regions, saving time and reducing the potential for manual errors.
 * - Integration: Facilitates the integration of gene region data with other genomic data formats like `.bed` files.
 * - Flexibility: Allows users to specify genes of interest directly via command-line arguments or standard input, making it adaptable to various workflows.
 * 
 * Overall, this tool streamlines the process of extracting and utilizing gene region information from GTF files, which is a common requirement in genomic research and bioinformatics analyses.
 */

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
    std::cerr << "Usage: gtf_file_region_retrieval [-g gtf_file] [-c head_count] [-f \"promoter\"] [--gene gene_name1 --gene gene_name2 ...]\n";
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
    std::string filter; // Modification of GeneRegion to perform
    std::vector<std::string> geneNames;
    int head = 0; // number of top hits to show

    // Option flags and variables for getopt
    int option;
    static struct option long_options[] = {
        {"gtf", required_argument, 0, 'g'},   // Path to GTF file
        {"gene", required_argument, 0, 'n'},  // Specify genes of interest
        {"filter", required_argument, 0, 'f'},  // Specify genes of interest
        {"head", required_argument, 0, 'c'}, // Head - show only specified number of first hits
        {"verbose", required_argument, 0, 'v'},  // Instructs to generate help info
        {"help", no_argument, 0, 'h'}, // Instructs to be chatty about what is happening
        {0, 0, 0, 0}
    };

    // Parse command line arguments using getopt
    while ((option = getopt_long(argc, argv, "c:f:g:n:vh", long_options, nullptr)) != -1) {
        switch (option) {
            case 'c':
                head = std::stoi(optarg);
                break;
            case 'f':
                filter = optarg;
                break;
            case 'g':
                gtfFile = optarg;
                break;
            case 'n': {
                std::istringstream iss(optarg);
                std::string gene;
                while (std::getline(iss, gene, ',')) {
                    std::istringstream iss2(gene);
                    std::string subGene;
                    while (iss2 >> subGene) {
                        std::cerr << "I: Adding gene name '" << subGene << "'" << std::endl;
                        geneNames.push_back(subGene);
                    }
                }
                break;
            }
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

    /*
    // If no gene names were provided, read from stdin
    if (geneNames.empty()) {
        geneNames = readGenesFromStdin();
    }

    if (geneNames.empty()) {
        std::cerr << "Error: No genes provided." << std::endl;
        return 1;
    }
    */
    
    // Parse the GTF file to get regions
    std::unordered_map<std::string, GeneRegion> geneRegions = parseGTFFile(gtfFile);
    std::cerr << "I: Found " << geneRegions.size() << " gene regions in the GTF file." << std::endl;

    // Create a map from gene name to vector of GeneRegion references
    std::unordered_map<std::string, std::vector<std::reference_wrapper<const GeneRegion>>> geneNameToRegions;

    // Populate the geneNameToRegions map
    for (const auto& [geneName, geneRegion] : geneRegions) {
        geneNameToRegions[geneRegion.geneName].push_back(geneRegion);
    }

    if (0 == geneNames.size()) {
        std::cerr << "I: No gene names provided, assigning all existing." << std::endl;
        for (const auto& element : geneNameToRegions) {
            geneNames.push_back(element.first);
        }
    }
    // Sort the gene names alphabetically
    std::sort(geneNames.begin(), geneNames.end());

    // Print the regions for the genes
    for (const std::string& geneName : geneNames) {
        //std::cerr << "D: testing for relevance of gene name '" << geneName << "'" << std::endl;
        int c = 0;
        // Check if the gene name exists in the map
        if (geneNameToRegions.find(geneName) != geneNameToRegions.end()) {
            std::vector<std::reference_wrapper<const GeneRegion>>& regions = geneNameToRegions[geneName];
            //std::cout << "# Genename: " << geneName << std::endl;
            // Print the regions for the gene
            for (const GeneRegion& region : regions) {
                if ("promoter" == filter) {
                    std::cout << region.relative_upstream(1,500).toBedString() << std::endl;
                }
                else {
                    std::cout << region.toBedString() << std::endl;
                }
                c++;
                if (head > 0 && c >= head) {
                    break; // Stop after printing the specified number of regions
                }
            }
        } else {
            std::cerr << "Gene " << geneName << " not found in the GTF file." << std::endl;
            std::cerr << "Got the following to choose from:" << std::endl;
            for (const auto& regionName : geneRegions) {
                std::cerr << regionName.first << std::endl;
            }
        }
    }

    return 0;
}
