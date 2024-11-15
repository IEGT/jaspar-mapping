#include <iostream>
#include <fstream>
#include <zlib.h>
#include <bzlib.h>
#include <memory>
#include <string>
#include <optional> // For handling unread lines

class CompressedFileReader {
public:
    enum FileType { PLAIN, GZIP, BZIP2 };

    explicit CompressedFileReader(const std::string& filename, FileType type)
        : type_(type), file_(nullptr), unread_line_(std::nullopt) {
        if (type == GZIP) {
            gz_file_ = gzopen(filename.c_str(), "rb");
            if (!gz_file_) throw std::runtime_error("Failed to open gzip file: " + filename);
        } else if (type == BZIP2) {
            bz_file_ = BZ2_bzopen(filename.c_str(), "rb");
            if (!bz_file_) throw std::runtime_error("Failed to open bzip2 file: " + filename);
        } else {
            file_ = std::make_unique<std::ifstream>(filename);
            if (!file_->is_open()) throw std::runtime_error("Failed to open file: " + filename);
        }
    }

    ~CompressedFileReader() {
        if (type_ == GZIP && gz_file_) gzclose(gz_file_);
        if (type_ == BZIP2 && bz_file_) BZ2_bzclose(bz_file_);
    }

    bool getline(std::string& line) {
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

    void unread(const std::string& line) {
        if (unread_line_) {
            throw std::runtime_error("Cannot unread multiple lines at once.");
        }
        unread_line_ = line; // Save the line for the next getline call
    }

private:
    FileType type_;
    std::unique_ptr<std::ifstream> file_;
    gzFile gz_file_ = nullptr;
    BZFILE* bz_file_ = nullptr;
    std::optional<std::string> unread_line_; // Buffer for the last read line
};
