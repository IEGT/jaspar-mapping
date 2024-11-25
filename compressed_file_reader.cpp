#include "compressed_file_reader.h"

/** \brief reads a line from the file
 * Is aware of single line that may hav been unread in earlier processing.
 * @param line - string to store the line
 * @return true if line was read, false if end of file
 */
bool CompressedFileReader::getline(std::string& line) {
    // If there's an unread line, return it first
    if (unread_line_) {
        line = *unread_line_;
        unread_line_.reset();
        return true;
    }

    char buffer[4096];
    if (type_ == GZIP) {
        if (!gzgets(gz_file_, buffer, sizeof(buffer))) return false;
    } else if (type_ == BZIP2) {
        int bytesRead = BZ2_bzread(bz_file_, buffer, sizeof(buffer) - 1);
        if (bytesRead <= 0) return false;
        buffer[bytesRead] = '\0'; // Ensure null-terminated string
    } else {
        if (!std::getline(*file_, line)) return false;
        return true;
    }

    line = buffer;
    // Remove trailing newline characters for consistency
    if (!line.empty() && line.back() == '\n') line.pop_back();
    if (!line.empty() && line.back() == '\r') line.pop_back();

    return true;
}

/** \brief Unread a single line, to be read again on the next getline call
 * No check on line separator or previous content performed.
 * @param line - string to unread
 * @throws runtime_error if trying to unread multiple lines
 */
void CompressedFileReader::unread(const std::string& line) {
    if (unread_line_) {
        throw std::runtime_error("Cannot unread multiple lines at once.");
    }
    unread_line_ = line; // Save the line for the next getline call
}

/** \brief Returns enum to distinguish plain file formats from gzip and bzip2
 * Static member function.
 * @param filename - path to file
 * @return FileType enum
 */
CompressedFileReader::FileType CompressedFileReader::determineFileType(const std::string& filename) {
    if (filename.size() >= 3 && filename.compare(filename.size() - 3, 3, ".gz") == 0) {
        return GZIP;
    } else if (filename.size() >= 4 && filename.compare(filename.size() - 4, 4, ".bz2") == 0) {
        return BZIP2;
    } else {
        return PLAIN;
    }
}