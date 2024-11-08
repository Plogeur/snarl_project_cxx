#include "snarl_parser.hpp"
#include "vcf_parser.hpp"
#include "matrix.hpp"
#include "binary_analysis.hpp"
#include "quantitative_analysis.hpp"

SnarlParser::SnarlParser(const std::string& vcf_path) : vcf_path(vcf_path) {
    VCFParser vcfParser(vcf_path);
    list_samples = vcfParser.getSampleNames();
    matrix = Matrix(100000000, list_samples.size()*2);
}

// Function to determine and extract an integer from the string
std::pair<int, int> determine_str(const std::string& s, int length_s, int i) {
    int current_int = 0;
    while (i < length_s && std::isdigit(s[i])) {
        current_int = current_int * 10 + (s[i] - '0');
        i++;
    }
    return std::pair<int, int>(i, current_int);
}

// Function to decompose a string with snarl information
std::vector<std::string> decompose_string(const std::string& s) {
    std::vector<std::string> result;
    size_t i = 0;
    int length_s = s.size();
    int prev_int = -1; // Placeholder for "None"
    char prev_sym = '\0'; // Placeholder for "None"

    while (i < length_s) {
        char start_sym = s[i];
        i++;
        const auto& [next_i, current_int] = determine_str(s, length_s, i);
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
std::vector<std::vector<std::string> > decompose_snarl(const std::vector<std::string>& lst) {
    std::vector<std::vector<std::string> > decomposed_list;
    for (const auto& s : lst) {
        decomposed_list.push_back(decompose_string(s));
    }
    return decomposed_list;
}

// Retrieve the index of `key` if it exists in `ordered_map`. Otherwise, add it and return the new index.
size_t getOrAddIndex(std::unordered_map<std::string, size_t>& orderedMap, const std::string& key, int lengthOrderedMap) {
    auto it = orderedMap.find(key);
    if (it != orderedMap.end()) {
        return it->second;
    } else {
        int newIndex = lengthOrderedMap;
        orderedMap[key] = newIndex;
        return newIndex;
    }
}

// Add True to the matrix if snarl is found
void pushMatrix(const std::string& decomposedSnarl, std::unordered_map<std::string, size_t>& rowHeaderDict, size_t indexColumn) {
    // Retrieve or add the index in one step and calculate length once
    int lengthOrderedMap = rowHeaderDict.size();
    size_t idxSnarl = getOrAddIndex(rowHeaderDict, decomposedSnarl, lengthOrderedMap);

    // Check if a new matrix chunk is needed
    size_t currentRowsNumber = matrix.rowCount();
    if (lengthOrderedMap > currentRowsNumber - 1) {
        matrix.expandMatrix();
    }

    // Add data to the matrix
    matrix.add_data(idxSnarl, indexColumn);
}

// Main function that parses the VCF file and fills the matrix
void fill_matrix() {
    VCFParser vcfParser(vcf_path);  // Create an instance of VCFParser
    const std::vector<std::string>& sampleNames = vcfParser.getSampleNames();
    std::unordered_map<std::string, size_t> row_header_dict;

    // Read and process variants
    while (vcfParser.hasNext()) {
        vcfParser.nextVariant();

        // Extract genotypes
        const std::vector<std::vector<int> >& genotypes = vcfParser.getGenotypes();

        // Extract and split `AT` field from the INFO field
        std::vector<std::string> snarl_list = vcfParser.getATInfo();
        std::vector<std::vector<std::string> > list_list_decomposed_snarl = decompose_snarl(snarl_list);

        // Process genotypes and fill matrix
        for (size_t index_column = 0; index_column < genotypes.size(); ++index_column) {
            int allele_1 = genotypes[index_column][0];
            int allele_2 = genotypes[index_column][1];
            size_t col_idx = index_column * 2;

            if (allele_1 == -1 || allele_2 == -1) { // Handle missing genotypes (./.)
                continue;
            }

            for (auto& decompose_allele_1 : list_list_decomposed_snarl[allele_1]) {
                pushMatrix(decompose_allele_1, row_header_dict, col_idx);
            }

            for (auto& decompose_allele_1 : list_list_decomposed_snarl[allele_1]) {
                pushMatrix(decompose_allele_1, row_header_dict, col_idx);
            }
        }
        matrix.set_row_header(row_header_dict);
    }
}

// Binary Table Generation
void binary_table(const std::unordered_map<std::string, std::vector<std::string>>& snarls,
                                  const std::unordered_map<std::string, bool>& binary_groups,
                                  const std::string& output) 
{
    std::ofstream outf(output, std::ios::binary);
    if (!outf.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    // Write headers
    std::string headers = "CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP_Fisher\tP_Chi2\tGROUP_1_PATH_1\tGROUP_1_PATH_2\tGROUP_2_PATH_1\tGROUP_2_PATH_2\n";
    outf.write(headers.c_str(), headers.size());

    // Iterate over each snarl
    for (const auto& [snarl, list_snarl] : snarls) {
        auto df = create_binary_table(binary_groups, list_snarl);
        auto stats = binary_stat_test(df);

        std::string chrom = "NA", pos = "NA", type_var = "NA", ref = "NA", alt = "NA";
        // Assuming stats is a vector containing the values in the order: fisher_p_value, chi2_p_value, GIPI, GIPII, GIIPI, GIIPII
        std::stringstream data;
        data << chrom << "\t" << pos << "\t" << snarl << "\t" << type_var << "\t" << ref << "\t" << alt
             << "\t" << stats[0] << "\t" << stats[1] << "\t" << stats[2] << "\t" << stats[3]
             << "\t" << stats[4] << "\t" << stats[5] << "\n";
        
        outf.write(data.str().c_str(), data.str().size());
    }
}

// Quantitative Table Generation
void quantitative_table(const std::unordered_map<std::string, std::vector<std::string>>& snarls,
                                        const std::unordered_map<std::string, float>& quantitative,
                                        const std::string& output) 
{
    std::ofstream outf(output, std::ios::binary);
    if (!outf.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    // Write headers
    std::string headers = "CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP\n";
    outf.write(headers.c_str(), headers.size());

    // Iterate over each snarl
    for (const auto& [snarl, list_snarl] : snarls) {
        auto df = create_quantitative_table(list_snarl);
        double pvalue = linear_regression(df, quantitative);

        std::string chrom = "NA", pos = "NA", type_var = "NA", ref = "NA", alt = "NA";
        std::stringstream data;
        data << chrom << "\t" << pos << "\t" << snarl << "\t" << type_var << "\t" << ref << "\t" << alt << "\t" << pvalue << "\n";

        outf.write(data.str().c_str(), data.str().size());
    }
}
