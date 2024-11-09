#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
#include "snarl_parser.hpp"     
#include "vcf_parser.hpp"
#include "matrix.hpp"
#include "arg_parser.hpp"

void print_help() {
    std::cout << "Usage: SnarlParser [options]\n\n"
              << "Options:\n"
              << "  --vcf_path <path>           Path to the VCF file (.vcf or .vcf.gz)\n"
              << "  --snarl <path>              Path to the snarl file (.txt or .tsv)\n"
              << "  -b, --binary <path>         Path to the binary group file (.txt or .tsv)\n"
              << "  -q, --quantitative <path>   Path to the quantitative phenotype file (.txt or .tsv)\n"
              << "  -o, --output <path>         Path to the output file\n"
              << "  -h, --help                  Print this help message\n";
}

int main(int argc, char* argv[]) {
    // Declare variables to hold argument values
    std::string vcf_path, snarl_path, binary_path, quantitative_path, output_path;
    bool show_help = false;

    // Parse arguments manually
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--vcf_path" && i + 1 < argc) {
            vcf_path = argv[++i];
        } else if (arg == "--snarl" && i + 1 < argc) {
            snarl_path = argv[++i];
        } else if ((arg == "-b" || arg == "--binary") && i + 1 < argc) {
            binary_path = argv[++i];
        } else if ((arg == "-q" || arg == "--quantitative") && i + 1 < argc) {
            quantitative_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            show_help = true;
        }
    }

    // Display help if requested or if required arguments are missing
    if (show_help || vcf_path.empty() || snarl_path.empty()) {
        print_help();
        return 0;
    }

    try {
        // Check format of the VCF file
        check_format_vcf_file(vcf_path);
        std::cout << "VCF file format is correct." << std::endl;

        // Check format of the group/snarl file
        check_format_group_snarl(snarl_path);
        std::cout << "Group/Snarl file format is correct." << std::endl;

        // Initialize the SnarlProcessor with the VCF path
        SnarlParser vcf_object(vcf_path);
        auto start_1 = std::chrono::high_resolution_clock::now();
        vcf_object.fill_matrix();
        auto end_1 = std::chrono::high_resolution_clock::now();
        std::cout << "Time Matrix: " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;

        // Parse the snarl file
        start_1 = std::chrono::high_resolution_clock::now();
        auto snarl = parse_snarl_path_file(snarl_path);

        // Process binary group file if provided
        if (!binary_path.empty()) {
            auto binary_group = parse_group_file(binary_path);

            if (!output_path.empty()) {
                vcf_object.binary_table(snarl, binary_group, output_path);
            } else {
                vcf_object.binary_table(snarl, binary_group);
            }
        }

        // Process quantitative phenotype file if provided
        if (!quantitative_path.empty()) {
            auto quantitative = parse_pheno_file(quantitative_path);

            if (!output_path.empty()) {
                vcf_object.quantitative_table(snarl, quantitative, output_path);
            } else {
                vcf_object.quantitative_table(snarl, quantitative);
            }
        }

        end_1 = std::chrono::high_resolution_clock::now();
        std::cout << "Time P-value: " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}