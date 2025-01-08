/** \file Implementation of PSSM class
*/

#include "pssm.h"

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
#include <math.h>

/** \brief Function to calculate log-odds score for a given nucleotide at a position
 */
double PSSM::logRelativeRisk(const double& frequency, const double& background) {
    if (frequency == 0) {
        return -1e9;  // Prevent log(0) by returning a large negative score for unobserved nucleotides
        //return -log2(1024*16);  // assume one of the next upcoming tests would have found that residue to avoid -inf
    }
    return log2(frequency / background);  // log Risk Ratio
}

/** \brief Function to calculate log-odds score for a given nucleotide at a position
 */
double PSSM::logRelativeRiskACGT(const double& frequency) {
    static double log2Background= -2 ; // == log2(0.25);
    if (frequency == 0) {
        return -1e9;  // Prevent log(0) by returning a large negative score for unobserved nucleotides
        //return -log2(1024*16);  // assume one of the next upcoming tests would have found that residue to avoid -inf
    }
    return log2(frequency) - log2Background;  // log Risk Ratio
}

double PSSM::logOddsRatioACGT(const double& count, const double& colsum) {

    static double log2BackgroundOdd = log2(1.0/(4.0-1.0));

    if (count == 0.0) {
        // not expected to happen because of numbers "smoothed" by 1 added to all counts
        return -1e9;  // Prevent log(0) by returning a large negative score for unobserved nucleotides
    }
    double remainer = colsum-count;
    if (remainer<0.0) {
        std::cerr << "E: PSSM::logOddsRatioACGT: remainer<0.0 for count " << count << " and colsum " << colsum << std::endl;
        exit(-1);
    }
    if (remainer == 0.0) {
        // not expected to happen because of numbers "smoothed" by 1 added to all counts
        return +1e9;  // Prevent frequency/0 division by returning a large negative score for unobserved nucleotides
    }
    return log2(count / remainer) - log2BackgroundOdd;  // Log-odds ratio, assuming equal distribution as background
}

/** \brief Function to trim whitespace from strings
 */
std::string PSSM::trim(const std::string& str) {
    const std::string whitespace = " \t\n\r";
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos) return "";  // No content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

/** \brief Function to parse a JASPAR PSSM file format
 * The function reads all PSSMs from the file and stores them in a map, that is unless targetMotifID is provided.
 * @param constpssmFile - path to JASPAR PSSM file to parse
 * @param pssm_list - the string->pssm mapping structure to read that PSSM file into
 * @param targetMotifID - the motif ID to search for in the PSSM file, may be empty
 * @returns 0 upon success, else -1.
 */
int PSSM::parsePSSMFile(const std::string& pssmFile, pssm_list_type& pssm_list, const std::string& targetMotifID, const int& beVerbose) {

    std::cerr << "D: Parsing file '" << pssmFile << "' aiming at target MotifID '" << targetMotifID << "' with verbosity level " << beVerbose << std::endl;

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
    int numPSSMsRead=0;
    while (getline(inFile, line)) {
        line = PSSM::trim(line);
        // Check for the FASTA-like header, indicating a new motif
        if ( ( !line.empty() && line[0] == '>')  || inFile.eof() ) {
            numPSSMsRead++;
            // Save the previous PSSM if there is one
            if (   !motifID.empty() // this is not the very first line read,
                                    // so a complete PSSM has already been seen
                && ( targetMotifID.empty() || motifID == targetMotifID )
                                    // no motif was specified or the desired motif was found
                ) {

                PSSM pssm_object(pssm, motifID, motifName, length);
                pssm_list[motifID] = pssm_object;
                if (beVerbose>0) std::cerr << "I: Got the following PSSM: " << std::endl << pssm_object;
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
            } else if (beVerbose>1 || debug) {
                std::cerr << "D: Skipping motif '" << motifID << "' that is != '" <<  targetMotifID << "'" << std:: endl;
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

            //if (beVerbose>1) std::cerr << "I: " << currentBase << " : ";
            std::string countsStr;
            std::getline(ss,countsStr);  // Read the bracketed counts (e.g., [ 4 19 0 0 ])
            //if (beVerbose>1) std::cerr << countsStr << std::endl;

            // Remove brackets and split counts
            countsStr.erase(std::remove(countsStr.begin(), countsStr.end(), '['), countsStr.end());
            countsStr.erase(std::remove(countsStr.begin(), countsStr.end(), ']'), countsStr.end());

            std::stringstream countStream(countsStr);
            double countValue;
            length=0;
            while (countStream >> countValue) {
                length++;
                pssm[currentBase].push_back(countValue);
                //if (beVerbose>1) std::cerr << " " << countValue;
                //if (beVerbose>1) std::cerr << " l:" << length;
            }
            //if (beVerbose>1) std::cerr << std::endl;
        } else {
            // Skip other lines
            if (beVerbose) {
                std::cerr << "I: Skipping line: " << line << std::endl;
            }
        }
    }

    // also inspecting last file, terminated by EOF, not by new ">" line
    if (   length>0
           && !motifID.empty() // this is not the very first line read,
                            // so a complete PSSM has already been seen
            && ( targetMotifID.empty() || motifID == targetMotifID )
                            // no motif was specified or the desired motif was found
       ) {
         PSSM pssm_object(pssm, motifID, motifName, length);
         pssm_list[motifID] = pssm_object;
         if (beVerbose>0) std::cerr << "I: EOF handling: Got the following PSSM << '" << motifID << "' for target '" << targetMotifID << "': " << std::endl << pssm_object;
         if (beVerbose) {
              std::cerr << "I: EOF handling: PSSM for motif " << motifID << " (" << motifName << ") saved with length " << length << "." << std::endl;
         }
    }

    std::cerr << "I: Read " << numPSSMsRead << " PSSMs from file '" << pssmFile << "'" << std::endl;

    inFile.close();
    return 0;
}

/** \brief Normalize the PSSM by converting counts to log-odds scores
 */
void PSSM::normalizePSSM(const std::unordered_map<char, const double>& backgroundFrequencies) {

        std::cerr << "I: NormalizePSSM " << std::endl;
        std::cerr << *this;

        for (auto& [nucleotide, counts] : this->pssm) {
                const auto background = backgroundFrequencies.at(nucleotide);

                //int i = 0;
                if (debug) std::cerr << "D: Counts.size= " << counts.size() << std::endl;
                //for (unsigned int i=0; auto& count : counts) {
                for (unsigned int i=0; i < counts.size(); i++) {
                        if (debug) std::cerr << "D: i=" << i << std::endl;
                        if (0 == this->colsums[i]) {
                                //count = 0;
                                counts[i] = 0;
                                std::cerr << "W: Unexpected - column " << i << " has no scores assigned." << std::endl;
                                continue;
                        }
/*
                        if (0.0 == counts[i]) { // an unobserved nucleotide
                                counts[i] = -2;
                                //count = -1e9;
                                continue;
                        }
*/
                        if (PSSM::debug) std::cerr << "I: Normalization of " << counts[i] << " counts at position " << i << " with column sum " << this->colsums[i] << " to ";
            // log odds ratio, smoothed by avoiding complete coverage
                        counts[i] = PSSM::logOddsRatioACGT(counts[i]+1, this->colsums[i]+4);
                        if (PSSM::debug) std::cerr << counts[i] << std::endl;
                }
        }
}

/** \brief debug flag
 * set to value > 0 to enable debug output
 */
int PSSM::debug = 0;
