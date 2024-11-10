#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

// Method 1: Buffered reading with custom buffer
void read_buffered_lines(const std::string& filename) {
    const size_t buffer_size = 1024 * 1024; // 1MB buffer
    std::ifstream file(filename);
    std::vector<char> buffer(buffer_size);

    while (file.read(buffer.data(), buffer_size) || file.gcount() > 0) {
        size_t buffer_end = file.gcount();
        size_t line_start = 0;

        for (size_t i = 0; i < buffer_end; ++i) {
            if (buffer[i] == '\n') {
                std::string line(buffer.begin() + line_start, buffer.begin() + i);
                std::cout << line << std::endl;  // Print each line
                line_start = i + 1; // Start next line
            }
        }

        if (line_start < buffer_end) {
            std::string remaining_line(buffer.begin() + line_start, buffer.begin() + buffer_end);
            std::cout << remaining_line << std::endl;  // Print the remaining line
        }
    }
}

// Method 2: Memory-mapped file for line parsing
void read_mmap_lines(const std::string& filename) {
    int fd = open(filename.c_str(), O_RDONLY);
    if (fd == -1) return;

    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        close(fd);
        return;
    }

    char* mapped = (char*)mmap(nullptr, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        close(fd);
        return;
    }

    char* start = mapped;
    char* end = mapped + sb.st_size;
    char* line_start = start;

    while (start < end) {
        if (*start == '\n') {
            std::string line(line_start, start - line_start);
            std::cout << line << std::endl;  // Print each line
            line_start = start + 1;
        }
        ++start;
    }

    if (line_start < end) {
        std::string remaining_line(line_start, end - line_start);
        std::cout << remaining_line << std::endl;  // Print the last line if no newline at the end
    }

    munmap(mapped, sb.st_size);
    close(fd);
}

// Method 3: Using std::istreambuf_iterator for line-by-line parsing
void read_using_iterator(const std::string& filename) {
    std::ifstream file(filename);
    std::istreambuf_iterator<char> start(file), end;

    std::string line;
    for (; start != end; ++start) {
        if (*start == '\n') {
            std::cout << line << std::endl;  // Print each line
            line.clear(); // Clear the current line
        } else {
            line.push_back(*start); // Add the character to the current line
        }
    }

    if (!line.empty()) {
        std::cout << line << std::endl;  // Print any remaining text
    }
}

// Method 4: Standard getline approach
void read_with_getline(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::cout << line << std::endl;  // Print each line
    }
}

// Benchmark each method
void benchmark(const std::string& filename) {
    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds duration;

    std::cout << "Benchmarking buffered reading...\n";
    start = std::chrono::high_resolution_clock::now();
    read_buffered_lines(filename);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now() - start);
    std::cout << "Buffered reading duration: " << duration.count() << " microseconds\n";

    std::cout << "Benchmarking memory-mapped reading...\n";
    start = std::chrono::high_resolution_clock::now();
    read_mmap_lines(filename);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now() - start);
    std::cout << "Memory-mapped reading duration: " << duration.count() << " microseconds\n";

    std::cout << "Benchmarking iterator reading...\n";
    start = std::chrono::high_resolution_clock::now();
    read_using_iterator(filename);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now() - start);
    std::cout << "Iterator reading duration: " << duration.count() << " microseconds\n";

    std::cout << "Benchmarking std::getline reading...\n";
    start = std::chrono::high_resolution_clock::now();
    read_with_getline(filename);
    duration = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now() - start);
    std::cout << "std::getline reading duration: " << duration.count() << " microseconds\n";
}

int main() {
    const std::string filename = "small_vcf.vcf";  // Replace with your file name
    benchmark(filename);
    return 0;
}
