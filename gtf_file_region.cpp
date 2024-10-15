
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <getopt.h>  // For GNU Getopt
#include "progress.h" // progress indicator
#include "gtf_file_region.h"  // GeneRegion struct

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
           const size_t diffExisting = existing.end > existing.start ? existing.end - existing.start : existing.start - existing.end;
           const size_t diff         =          end >          start ?          end -          start :          start -          end;
           if ( diffExisting > diff ) {
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