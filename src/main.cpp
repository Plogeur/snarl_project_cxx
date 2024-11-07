#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
#include "cxxopts.hpp"          
#include "snarl_parser.hpp"     
#include "vcf_parser.hpp"
#include "matrix.hpp"
#include "arg_parser.hpp"

// Function prototypes for file checking and parsing (implement these separately)
bool check_format_vcf_file(const std::string& path);
bool check_format_group_snarl(const std::string& path);
bool check_format_pheno_b(const std::string& path);
bool check_format_pheno_q(const std::string& path);
std::unordered_map<std::string, size_t> parse_snarl_path_file(const std::string& path);
std::vector<std::string> parse_group_file(const std::string& path);
std::vector<double> parse_pheno_file(const std::string& path);

int main(int argc, char* argv[]) {
    try {
        // Parse command-line arguments
        cxxopts::Options options("SnarlParser", "Parse and analyze snarl from VCF file");
        options.add_options()
            ("vcf_path", "Path to the VCF file (.vcf or .vcf.gz)", cxxopts::value<std::string>())
            ("snarl", "Path to the snarl file containing snarls and aT (.txt or .tsv)", cxxopts::value<std::string>())
            ("b,binary", "Path to the binary group file (.txt or .tsv)", cxxopts::value<std::string>())
            ("q,quantitative", "Path to the quantitative phenotype file (.txt or .tsv)", cxxopts::value<std::string>())
            ("o,output", "Path to the output file", cxxopts::value<std::string>()->default_value(""))
            ("help", "Print help");

        auto result = options.parse(argc, argv);

        // Display help if requested
        if (result.count("help") || argc == 1) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        // Extract mandatory arguments
        std::string vcf_path = result["vcf_path"].as<std::string>();
        std::string snarl_path = result["snarl"].as<std::string>();

        // Validate file formats
        if (!check_format_vcf_file(vcf_path) || !check_format_group_snarl(snarl_path)) {
            std::cerr << "Error: Invalid file format for VCF or Snarl file." << std::endl;
            return 1;
        }

        // Initialize the SnarlProcessor with the VCF path
        auto start = std::chrono::high_resolution_clock::now();
        SnarlProcessor vcf_object(vcf_path);
        vcf_object.fill_matrix();
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Time Matrix: " << std::chrono::duration<double>(end - start).count() << " s" << std::endl;

        // Parse the snarl file
        start = std::chrono::high_resolution_clock::now();
        auto snarl = parse_snarl_path_file(snarl_path);

        // Process binary group file if provided
        if (result.count("binary")) {
            std::string binary_path = result["binary"].as<std::string>();
            auto binary_group = parse_group_file(binary_path);

            if (result.count("output")) {
                std::string output_path = result["output"].as<std::string>();
                vcf_object.binary_table(snarl, binary_group, output_path);
            } else {
                vcf_object.binary_table(snarl, binary_group);
            }
        }

        // Process quantitative phenotype file if provided
        if (result.count("quantitative")) {
            std::string quantitative_path = result["quantitative"].as<std::string>();
            auto quantitative = parse_pheno_file(quantitative_path);

            if (result.count("output")) {
                std::string output_path = result["output"].as<std::string>();
                vcf_object.quantitative_table(snarl, quantitative, output_path);
            } else {
                vcf_object.quantitative_table(snarl, quantitative);
            }
        }

        end = std::chrono::high_resolution_clock::now();
        std::cout << "Time P-value: " << std::chrono::duration<double>(end - start).count() << " s" << std::endl;

    } catch (const cxxopts::OptionException& e) {
        std::cerr << "Error parsing options: " << e.what() << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
