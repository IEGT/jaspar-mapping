#include "compressed_file_reader.h"

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

void CompressedFileReader::unread(const std::string& line) {
    if (unread_line_) {
        throw std::runtime_error("Cannot unread multiple lines at once.");
    }
    unread_line_ = line; // Save the line for the next getline call
}
