#include <string>
#include <iostream>

// Structure to store gene region information
struct GeneRegion {
    std::string chromosome;
    size_t start;
    size_t end;
    std::string strand;
    std::string geneId;
    std::string geneName;
    std::string transcriptName;

    // Return content as a line that can be written to a BED file
    std::string toBedString() const;

    // Return a new GeneRegion object that represents the region upstream of the current region
    GeneRegion relative_upstream(size_t min_upstream, size_t max_upstream) const;

    // Return a new GeneRegion object that represents the region downstream of the current region
    GeneRegion relative_downstream(size_t min_downstream, size_t max_downstream) const;

    // Overload the << operator for GeneRegion
    friend std::ostream& operator<<(std::ostream& os, const GeneRegion& region);
};

/**
 * @brief Parses a GTF file to extract regions for specific gene names.
 *
 * This function reads a GTF (Gene Transfer Format) file and extracts gene regions
 * for specific gene names. It stores the gene regions in an unordered map where
 * the key is the gene identifier (either gene name or gene ID) and the value is
 * a GeneRegion structure containing the chromosome, start, end, and gene identifier.
 *
 * @param gtfFile The path to the GTF file to be parsed.
 * @return An unordered map where the key is the gene identifier (either gene name or gene ID)
 *         and the value is a GeneRegion structure containing the chromosome, start, end,
 *         and gene identifier.
 */
std::unordered_map<std::string, GeneRegion> parseGTFFile(const std::string& gtfFile);

// Function to read gene names from stdin
std::vector<std::string> readGenesFromStdin();