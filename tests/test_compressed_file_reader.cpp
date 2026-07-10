#include "../compressed_file_reader.h"

#include <bzlib.h>
#include <zlib.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

void require(const bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

struct TemporaryDirectory {
    std::filesystem::path path;

    ~TemporaryDirectory() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

std::string joinLines(const std::vector<std::string>& lines) {
    std::string content;
    for (size_t i = 0; i < lines.size(); ++i) {
        if (i != 0) {
            content.push_back('\n');
        }
        content += lines[i];
    }
    return content;
}

void writePlainFile(const std::filesystem::path& path, const std::string& content) {
    std::ofstream output(path, std::ios::binary);
    require(output.is_open(), "Could not create plain-text fixture");
    output.write(content.data(), static_cast<std::streamsize>(content.size()));
    require(output.good(), "Could not write plain-text fixture");
}

void writeGzipFile(const std::filesystem::path& path, const std::string& content) {
    gzFile output = gzopen(path.string().c_str(), "wb");
    require(output != nullptr, "Could not create gzip fixture");

    const int bytesWritten = gzwrite(
        output, content.data(), static_cast<unsigned int>(content.size()));
    const int closeResult = gzclose(output);
    require(bytesWritten == static_cast<int>(content.size()), "Could not write complete gzip fixture");
    require(closeResult == Z_OK, "Could not close gzip fixture");
}

void writeBzip2File(const std::filesystem::path& path, const std::string& content) {
    unsigned int compressedSize = static_cast<unsigned int>(content.size() + content.size() / 100 + 601);
    std::vector<char> compressed(compressedSize);
    std::vector<char> source(content.begin(), content.end());
    const int result = BZ2_bzBuffToBuffCompress(
        compressed.data(), &compressedSize, source.data(),
        static_cast<unsigned int>(content.size()), 9, 0, 30);
    require(result == BZ_OK, "Could not compress bzip2 fixture");

    std::ofstream output(path, std::ios::binary);
    require(output.is_open(), "Could not create bzip2 fixture");
    output.write(compressed.data(), static_cast<std::streamsize>(compressedSize));
    require(output.good(), "Could not write bzip2 fixture");
}

void expectLines(CompressedFileReader& reader, const std::vector<std::string>& expected) {
    std::string line;
    for (size_t i = 0; i < expected.size(); ++i) {
        require(reader.getline(line), "Reader reached EOF before line " + std::to_string(i));
        require(line == expected[i],
                "Unexpected content for line " + std::to_string(i)
                + ": expected length " + std::to_string(expected[i].size())
                + ", got " + std::to_string(line.size()));
    }
    require(!reader.getline(line), "Reader returned an extra line after EOF");
}

void expectFileLines(const std::filesystem::path& path,
                     const std::vector<std::string>& expected) {
    CompressedFileReader reader(path.string(), false);
    expectLines(reader, expected);
}

void expectMoveConstruction(const std::filesystem::path& path,
                            const std::vector<std::string>& expected) {
    CompressedFileReader source(path.string(), false);
    std::string firstLine;
    require(source.getline(firstLine), "Move-construction source had no first line");
    source.unread(firstLine);

    CompressedFileReader destination(std::move(source));
    expectLines(destination, expected);
}

void expectMoveAssignment(const std::filesystem::path& sourcePath,
                          const std::filesystem::path& destinationPath,
                          const std::vector<std::string>& expected) {
    CompressedFileReader source(sourcePath.string(), false);
    std::string firstLine;
    require(source.getline(firstLine), "Move-assignment source had no first line");
    source.unread(firstLine);

    {
        CompressedFileReader destination(destinationPath.string(), false);
        destination = std::move(source);
        expectLines(destination, expected);
    }
}

} // namespace

int main() {
    try {
        const auto suffix = std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count());
        TemporaryDirectory temporary{
            std::filesystem::temp_directory_path() / ("jaspar-compressed-reader-" + suffix)};
        require(std::filesystem::create_directory(temporary.path),
                "Could not create temporary test directory");

        const std::vector<std::string> expected = {
            "alpha",
            std::string(10000, 'X'),
            "after-long-line",
            "",
            "final-line-without-newline"
        };
        const std::string content = joinLines(expected);
        const std::filesystem::path plainPath = temporary.path / "fixture.txt";
        const std::filesystem::path gzipPath = temporary.path / "fixture.txt.gz";
        const std::filesystem::path bzip2Path = temporary.path / "fixture.txt.bz2";

        writePlainFile(plainPath, content);
        writeGzipFile(gzipPath, content);
        writeBzip2File(bzip2Path, content);

        expectFileLines(plainPath, expected);
        expectFileLines(gzipPath, expected);
        expectFileLines(bzip2Path, expected);
        expectMoveConstruction(gzipPath, expected);
        expectMoveAssignment(gzipPath, bzip2Path, expected);
        expectMoveAssignment(bzip2Path, gzipPath, expected);

        std::cout << "All compressed-file reader tests passed." << std::endl;
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Compressed-file reader test failed: " << error.what() << std::endl;
        return 1;
    }
}
