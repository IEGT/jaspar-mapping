#ifndef PSSM_H
#define PSSM_H

#include <string>
#include <vector>
#include <unordered_map>
#include <ostream> // std::endl
#include <iostream> // std::cerr

class PSSM;
typedef std::unordered_map<char, std::vector<double> > pssm_type;    //< Map of nucleotide to counts

/**
 * A class representing a PSSM (Position-Specific Scoring Matrix).
 * Stores nucleotide counts and associated motif information.
 */
class PSSM {
    public:
        pssm_type pssm;
        std::vector<double> colsums;
        std::string motifID;
        std::string motifName;
        int motifLength;
        /**
         * Default constructor that initializes empty PSSM and motif information.
         */
        PSSM() {
            pssm.clear();
            motifID.clear();
            motifName.clear();
            motifLength = 0;
            colsums.clear();
        }

        PSSM(const PSSM& c) {
            std::cerr << "D: PSSM:invoked copy constructor" << std::endl;
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
        PSSM(const pssm_type& pssm, const std::string& motifID, const std::string& motifName, const int& motifLength) {
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
 
        ~PSSM() {
            std::cerr << "PSSM - deleting colsums[]" << std::endl;
        }

        static int parsePSSMFile(const std::string& pssmFile, std::unordered_map<std::string, PSSM>& pssm_list, const std::string& targetMotifID, const int& beVerbose=0);
        static std::string trim(const std::string& str);
        static inline double logOddsScore(const double& frequency, const double& background);

        // Normalize the PSSM by converting counts to log-odds scores
        void normalizePSSM(const std::unordered_map<char, const double>& backgroundFrequencies);

        // Overload the << operator for PSSM
        friend std::ostream& operator<<(std::ostream& os, const PSSM& pssmObj) {
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

typedef std::unordered_map<std::string, PSSM> pssm_list_type;

#endif
