#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include <cctype>
#include "vcf_parser.hpp"

std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    size_t start = 0, end;
    while ((end = str.find(delimiter, start)) != std::string::npos) {
        tokens.push_back(str.substr(start, end - start));
        start = end + 1;
    }
    tokens.push_back(str.substr(start));
    return tokens;
}

// Function to determine and extract an integer from the string
std::pair<int, int> determine_str(const std::string& s, int length_s, int i) {
    int current_int = 0;
    while (i < length_s && std::isdigit(s[i])) {
        current_int = current_int * 10 + (s[i] - '0');
        i++;
    }
    return {i, current_int};
}

// Function to decompose a string with snarl information
std::vector<std::string> decompose_string(const std::string& s) {
    std::vector<std::string> result;
    int i = 0;
    int length_s = s.size();
    int prev_int = -1; // Placeholder for "None"
    char prev_sym = '\0'; // Placeholder for "None"

    while (i < length_s) {
        char start_sym = s[i];
        i++;
        auto [next_i, current_int] = determine_str(s, length_s, i);
        i = next_i;

        if (prev_int != -1 && prev_sym != '\0') {
            result.push_back(std::string(1, prev_sym) + std::to_string(prev_int) + start_sym + std::to_string(current_int));
        }

        prev_int = current_int;
        prev_sym = start_sym;
    }
    return result;
}

// Function to decompose a list of snarl strings
std::vector<std::vector<std::string>> decompose_snarl(const std::vector<std::string>& lst) {
    std::vector<std::vector<std::string>> decomposed_list;
    for (const auto& s : lst) {
        decomposed_list.push_back(decompose_string(s));
    }
    return decomposed_list;
}

// Main function that parses the VCF file and fills the matrix
void fill_matrix(const std::string& vcf_path) {
    VCFParser vcfParser(vcf_path);  // Create an instance of VCFParser
    const auto& sampleNames = vcfParser.getSampleNames(); // Accessing sample names

    // Read and process variants
    while (vcfParser.hasNext()) {
        vcfParser.nextVariant();

        // Extract genotypes
        std::vector<std::pair<int, int>> genotypes;
        const auto& genotypes = vcfParser.getGenotypes();

        // Extract and split `AT` field from the INFO field
        std::vector<std::string> snarl_list = vcfParser.getATInfo();
        std::vector<std::vector<std::string>> list_list_decomposed_snarl = decompose_snarl(snarl_list);

        // Process genotypes and fill matrix
        for (size_t index_column = 0; index_column < genotypes.size(); ++index_column) {
            int allele_1 = genotypes[index_column].first;
            int allele_2 = genotypes[index_column].second;
            size_t col_idx = index_column * 2;

            if (allele_1 == -1 || allele_2 == -1) { // Handle missing genotypes (./.)
                continue;
            }

            std::cout << allele_1 << allele_2 << std::endl;
        }

        for (const auto& snarl : list_list_decomposed_snarl) {
            for (const auto& component : snarl) {
                std::cout << component << " ";
            }
            std::cout << std::endl;
        }
    }
}
