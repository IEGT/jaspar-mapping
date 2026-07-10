#include <cstddef>
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <bzlib.h>
#include <memory>
#include <optional> // For handling unread lines
#include <stdexcept>
#include <string>

class CompressedFileReader {
public:
    enum FileType { PLAIN, GZIP, BZIP2 };

    explicit CompressedFileReader(const std::string& filename, const int verbose=true)
        : CompressedFileReader(filename, determineFileType(filename), verbose) {}

    explicit CompressedFileReader(const std::string& filename, FileType type, const int verbose=true)
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
        if (verbose) std::cerr << "I: Opened " << filename << std::endl;
    }

    // Delete copy constructor and copy assignment operator
    CompressedFileReader(const CompressedFileReader&) = delete;
    CompressedFileReader& operator=(const CompressedFileReader&) = delete;

    CompressedFileReader(CompressedFileReader&& other) noexcept;
    CompressedFileReader& operator=(CompressedFileReader&& other) noexcept;

    static CompressedFileReader::FileType determineFileType(const std::string& filename);

    /** Fixme: not implemented for bzip2 */
    std::streampos getFileSize() const {
        if (type_ == PLAIN && file_) {
            std::streampos current_pos = file_->tellg();
            file_->seekg(0, std::ios::end);
            std::streampos size = file_->tellg();
            file_->seekg(current_pos);
            return size;
        } else if (type_ == GZIP && gz_file_) {
            // Gzip file size is not straightforward to determine
            // This is a workaround to get the uncompressed size
            gzseek(gz_file_, 0, SEEK_END);
            std::streampos size = gztell(gz_file_);
            gzrewind(gz_file_);
            return size;
        } else if (type_ == BZIP2 && bz_file_) {
            // Bzip2 file size is not straightforward to determine
            // This is a workaround to get the uncompressed size
           /* BZ2_bzseek(bz_file_, 0, SEEK_END);
            std::streampos size = BZ2_bztell(bz_file_);
            BZ2_bzrewind(bz_file_);
            return size;*/
            return -1;
        }
        return -1;
    }

    /** Fixme: not implemented for bzip2 */
    std::streampos getBytesRead() const {
        if (type_ == PLAIN && file_) {
            return file_->tellg();
        } else if (type_ == GZIP && gz_file_) {
            return gztell(gz_file_);
        } else if (type_ == BZIP2 && bz_file_) {
            //return BZ2_bztell(bz_file_);
            return -1;
        }
        return -1;
    }

    ~CompressedFileReader();

    /** Read a line from the file, return false if EOF */
    bool getline(std::string& line);

    /** Unread a line, to be read again on the next getline call */
    void unread(const std::string& line);

private:
    bool getlineCompressed(std::string& line);
    void close() noexcept;
    void moveFrom(CompressedFileReader&& other) noexcept;

    FileType type_ = PLAIN;
    std::unique_ptr<std::ifstream> file_;
    gzFile gz_file_ = nullptr;
    BZFILE* bz_file_ = nullptr;
    std::optional<std::string> unread_line_; // Buffer for the last read line
    std::string compressed_buffer_;
    size_t compressed_buffer_offset_ = 0;
    bool compressed_eof_ = false;
};
