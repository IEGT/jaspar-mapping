#include "compressed_file_reader.h"

#include <utility>

namespace {

constexpr size_t COMPRESSED_READ_BUFFER_SIZE = 4096;

void removeTrailingCarriageReturn(std::string& line) {
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
}

} // namespace

CompressedFileReader::CompressedFileReader(CompressedFileReader&& other) noexcept {
    moveFrom(std::move(other));
}

CompressedFileReader& CompressedFileReader::operator=(CompressedFileReader&& other) noexcept {
    if (this != &other) {
        close();
        moveFrom(std::move(other));
    }
    return *this;
}

CompressedFileReader::~CompressedFileReader() {
    close();
}

void CompressedFileReader::close() noexcept {
    if (gz_file_) {
        gzclose(gz_file_);
        gz_file_ = nullptr;
    }
    if (bz_file_) {
        BZ2_bzclose(bz_file_);
        bz_file_ = nullptr;
    }
    if (file_) {
        file_->close();
        file_.reset();
    }
    unread_line_.reset();
    compressed_buffer_.clear();
    compressed_buffer_offset_ = 0;
    compressed_eof_ = true;
}

void CompressedFileReader::moveFrom(CompressedFileReader&& other) noexcept {
    type_ = other.type_;
    file_ = std::move(other.file_);
    gz_file_ = other.gz_file_;
    bz_file_ = other.bz_file_;
    unread_line_ = std::move(other.unread_line_);
    compressed_buffer_ = std::move(other.compressed_buffer_);
    compressed_buffer_offset_ = other.compressed_buffer_offset_;
    compressed_eof_ = other.compressed_eof_;

    other.type_ = PLAIN;
    other.gz_file_ = nullptr;
    other.bz_file_ = nullptr;
    other.unread_line_.reset();
    other.compressed_buffer_.clear();
    other.compressed_buffer_offset_ = 0;
    other.compressed_eof_ = true;
}

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

    if (type_ == PLAIN) {
        if (!file_) return false;
        if (!std::getline(*file_, line)) return false;
        return true;
    }

    return getlineCompressed(line);
}

bool CompressedFileReader::getlineCompressed(std::string& line) {
    while (true) {
        const size_t newline = compressed_buffer_.find('\n', compressed_buffer_offset_);
        if (newline != std::string::npos) {
            line.assign(compressed_buffer_, compressed_buffer_offset_,
                        newline - compressed_buffer_offset_);
            compressed_buffer_offset_ = newline + 1;
            removeTrailingCarriageReturn(line);
            return true;
        }

        if (compressed_eof_) {
            if (compressed_buffer_offset_ == compressed_buffer_.size()) {
                compressed_buffer_.clear();
                compressed_buffer_offset_ = 0;
                return false;
            }
            line.assign(compressed_buffer_, compressed_buffer_offset_, std::string::npos);
            compressed_buffer_.clear();
            compressed_buffer_offset_ = 0;
            removeTrailingCarriageReturn(line);
            return true;
        }

        if (compressed_buffer_offset_ != 0) {
            compressed_buffer_.erase(0, compressed_buffer_offset_);
            compressed_buffer_offset_ = 0;
        }

        char buffer[COMPRESSED_READ_BUFFER_SIZE];
        int bytesRead = 0;
        if (type_ == GZIP) {
            bytesRead = gzread(gz_file_, buffer, sizeof(buffer));
            if (bytesRead < 0) {
                int errorCode = Z_OK;
                const char* errorMessage = gzerror(gz_file_, &errorCode);
                throw std::runtime_error(
                    std::string("Failed to read gzip stream: ")
                    + (errorMessage ? errorMessage : "unknown zlib error"));
            }
            compressed_eof_ = bytesRead == 0;
        } else {
            bytesRead = BZ2_bzread(bz_file_, buffer, sizeof(buffer));
            int errorCode = BZ_OK;
            const char* errorMessage = BZ2_bzerror(bz_file_, &errorCode);
            if (errorCode != BZ_OK && errorCode != BZ_STREAM_END) {
                throw std::runtime_error(
                    std::string("Failed to read bzip2 stream: ")
                    + (errorMessage ? errorMessage : "unknown bzip2 error"));
            }
            compressed_eof_ = errorCode == BZ_STREAM_END || bytesRead == 0;
        }

        if (bytesRead > 0) {
            compressed_buffer_.append(buffer, static_cast<size_t>(bytesRead));
        }
    }
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
