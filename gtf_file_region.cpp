#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <getopt.h>  // For GNU Getopt
#include <cassert>
#include "progress.h" // progress indicator
#include "gtf_file_region.h"  // GeneRegion struct


void storeGeneReference(std::unordered_map<std::string, GeneRegion>& geneRegions,
    const std::string& geneId,
    const std::string& geneName,
    const std::string& transcriptName,
    std::string& chrom, const size_t& start, const size_t& end, const std::string& strand,
    // source file line number
    const size_t line) {
        // If the gene name is empty, use the transcript name

        const std::string& geneIdentifier = transcriptName;

        // Check if the gene already exists in geneRegions
        if (geneRegions.find(geneIdentifier) != geneRegions.end()) {
           //std::cerr << "D: Gene " << geneIdentifier << " already exists in geneRegions." << std::endl;
           GeneRegion existing = geneRegions[geneIdentifier];
           // Bail out if two genes with the same name are found to exits on different chromosomes.
           if (existing.chromosome != chrom) {
               std::cerr << "E: Inconsistency for gene " << geneIdentifier << " - existing.chromosome = " << existing.chromosome << " != " << chrom << " = chrom" << " in line " << line << "." << std::endl;
               std::cerr << "   Existing:  " << existing << std::endl;
               std::cerr << "   Candidate: " << GeneRegion{chrom, start, end, strand, geneId, geneName, transcriptName} << std::endl;
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
        geneRegions[geneIdentifier] = {chrom, start, end, strand, geneId, geneName, transcriptName};
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

    // The following are expected to be set in every line of the GTF file.
    // If not, especially for the gene name, it was set in the line above and is carried over.
    std::string geneName; 

    while (getline(inFile, line)) {
        std::string geneID, transcriptName;

        lineNo++;
        std::istringstream ss(line);
        bytesRead += line.size();

                     //1   havana	gene    2581560   2584533  .      +       .      gene_id "ENSG00000228037"; gene_version "1"; gene_source "havana"; gene_biotype "lncRNA";
        std::string chrom, source, feature, startStr, endStr, score, strand, frame, attributes;
        if (!(ss >> chrom >> source >> feature >> startStr >> endStr >> score >> strand >> frame)) {
            continue;  // Skip lines that don't match the expected format
        }
        if (chrom.starts_with('#') || chrom.starts_with("KI") || chrom.starts_with("GL")) {
            continue;  // Skip header / comment / unwanted chromosomes
        }
        std::getline(ss,attributes);

        if (
            feature != "gene" && feature != "transcript"
        ) {
            continue; // Skip lines that are not gene or transcript entries
        }

        long start = std::stol(startStr);
        long end = std::stol(endStr);

        // Extract the gene name from the attributes field
        size_t biotypePos = attributes.find("gene_biotype");
        std::string biotype;
        if (biotypePos != std::string::npos) {
            size_t startQuote = attributes.find("\"", biotypePos);
            if (std::string::npos != startQuote) {
                size_t endQuote = attributes.find("\"", startQuote + 1);
                biotype = attributes.substr(startQuote + 1, endQuote - startQuote - 1);
            } else {
                biotype.clear();
            }
        }

        //std::cerr << "I: biotype = " << biotype << ", feature = " << feature << std::endl;
        
        if ("gene" == feature) {
            if ("unprocessed_pseudogene"==biotype || "processed_pseudogene"==biotype || "TEC" == biotype || "lncRNA" == biotype || "misc_RNA" == biotype || "snoRNA" == biotype) {
                geneName="lncRNA_"+chrom+"_"+startStr+"_"+endStr;
                continue; // transcript is redundant
            } else {
                //geneName.clear(); // not necessarily set in gene line
            }
        }

        if (feature != "transcript") {
            continue; // Skip lines that are not gene or transcript entries
        }

        if ("unprocessed_pseudogene"==biotype || "processed_pseudogene"==biotype || "TEC" == biotype || "lncRNA" == biotype || "misc_RNA" == biotype) {
            geneName="lncRNA_"+chrom+"_"+startStr+"_"+endStr;
            continue; // transcript is redundant
        } else {
            // Extract the gene name from the attributes field
            size_t geneNamePos = attributes.find("gene_name");
            if (geneNamePos != std::string::npos) {
                size_t startQuote = attributes.find("\"", geneNamePos);
                if (std::string::npos != startQuote) {
                    size_t endQuote = attributes.find("\"", startQuote + 1);
                    geneName = attributes.substr(startQuote + 1, endQuote - startQuote - 1);
                    //if ("TP73" == geneName) {
                    //    std::cerr << "I: found gene name '" << geneName << "' at line " << lineNo << " from " << startQuote << " to " << endQuote << "." << std::endl;
                    //}
                } else {
                    //std::cerr << "D: endQuote: " << endQuote << std::endl;
                    // Not clearing the name, a transcript may reference only the ID.
                    geneName.clear();
                }
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

        // Extract the gene name from the attributes field
        size_t transcriptNamePos = attributes.find("transcript_name");
        if (transcriptNamePos != std::string::npos) {
            size_t startQuote = attributes.find("\"", transcriptNamePos);
            if (std::string::npos != startQuote) {
                size_t endQuote = attributes.find("\"", startQuote + 1);
                transcriptName = attributes.substr(startQuote + 1, endQuote - startQuote - 1);
            } else {
                //std::cerr << "D: endQuote: " << endQuote << std::endl;
                transcriptName.clear();
            }
        }

        if (geneID.empty()) {
            std::cerr << "E: Line " << lineNo << ": No gene ID found in line " << lineNo << ": '" << line << "'" << std::endl;
            exit(1);
        }

        if (transcriptName.empty()) {
            //std::cerr << "D: transcriptName is empty, using geneName for line " << lineNo << ": '" << line << "'" << std::endl;
            transcriptName = geneName;
        }

        if (geneName.empty()) {
            std::cerr << "E: Line " << lineNo << ": No gene name found in line '" << line << "'" << std::endl;
            std::cerr << "   GeneName was expected to be carried over from the line before." << std::endl;
            exit(1);
        }

        if (transcriptName.empty()) {
            storeGeneReference(geneRegions, geneID, geneName, geneName, chrom, start, end, strand, lineNo);
        }
        else {
            storeGeneReference(geneRegions, geneID, geneName, transcriptName, chrom, start, end, strand, lineNo);
        }

        if (lineNo % 100 == 0) {
            float progress = static_cast<float>(bytesRead) / gtfFileSize;
            displayProgressBar(progress);
        }

    }
    inFile.close();
    return geneRegions;
}

// Return a new GeneRegion object that represents the region upstream of the current region
GeneRegion GeneRegion::relative_upstream(size_t minDist=1, size_t maxDist=500) const {
    assert(end>start);
    if (minDist > maxDist) {
        std::cerr << "E GeneRegion::relative_upstream: minDist > maxDist" << std::endl;
        abort();
    }
    if ("+" == strand) {
        size_t newStart = start - maxDist;
        size_t newEnd = start - minDist;
        return {chromosome, newStart, newEnd, strand, geneId, geneName, transcriptName};
    }
    else if ("-" == strand) {
        size_t newStart = end + minDist;
        size_t newEnd = end + maxDist;
        return {chromosome, newStart, newEnd, strand, geneId, geneName, transcriptName};
    }
    std::cerr << "E GeneRegion::relative_upstream: unknown strand: '" << strand << "'" << std::endl;
    abort();
}

std::string GeneRegion::toBedString() const {
    return chromosome + "\t" + std::to_string(start) + "\t" + std::to_string(end) + "\t" + geneName + "\t" + "0" + "\t" + strand;
}

// Overload the << operator for GeneRegion
std::ostream& operator<<(std::ostream& os, const GeneRegion& region)
{
    os << "GeneRegion(chromosome: " << region.chromosome
       << ", start: " << region.start
       << ", end: " << region.end
       << ", strand: " << region.strand
       << ", geneId: " << region.geneId
       << ", geneName: " << region.geneName
       << ", transcriptName: " << region.transcriptName
       << ")";
    return os;
}