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
std::pair<size_t, size_t> determine_str(const std::string& s, size_t length_s, size_t i) {
    size_t current_int = 0;
    while (i < length_s && std::isdigit(s[i])) {
        current_int = current_int * 10 + (s[i] - '0');
        i++;
    }
    return std::pair<size_t, size_t>(i, current_int);
}

// Function to decompose a string with snarl information
std::vector<std::string> decompose_string(const std::string& s) {
    std::vector<std::string> result;
    size_t i = 0;
    size_t length_s = s.size();
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
std::vector<std::vector<std::string>> decompose_snarl(const std::vector<std::string>& lst) {
    std::vector<std::vector<std::string>> decomposed_list;
    for (const auto& s : lst) {
        decomposed_list.push_back(decompose_string(s));
    }
    return decomposed_list;
}

// Retrieve the index of `key` if it exists in `ordered_map`. Otherwise, add it and return the new index.
size_t getOrAddIndex(std::unordered_map<std::string, size_t>& orderedMap, const std::string& key, size_t lengthOrderedMap) {
    auto it = orderedMap.find(key);
    if (it != orderedMap.end()) {
        return it->second;
    } else {
        size_t newIndex = lengthOrderedMap;
        orderedMap[key] = newIndex;
        return newIndex;
    }
}

// Add True to the matrix if snarl is found
void SnarlParser::pushMatrix(const std::string& decomposedSnarl, std::unordered_map<std::string, size_t>& rowHeaderDict, size_t indexColumn) {
    // Retrieve or add the index in one step and calculate length once
    size_t lengthOrderedMap = rowHeaderDict.size();
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
void SnarlParser::fill_matrix() {
    VCFParser vcfParser(vcf_path);
    std::unordered_map<std::string, size_t> row_header_dict;

    // Read and process variants
    while (vcfParser.hasNext()) {
        vcfParser.nextVariant();

        // Extract genotypes
        const std::vector<std::vector<int> >& genotypes = vcfParser.getGenotypes();

        // Extract and split `AT` field from the INFO field
        std::vector<std::string> snarl_list = vcfParser.getATInfo();
        std::vector<std::vector<std::string>> list_list_decomposed_snarl = decompose_snarl(snarl_list);

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

std::vector<int> identify_correct_path(
    const std::vector<std::string>& decomposed_snarl, 
    const std::unordered_map<std::string, size_t>& row_headers_dict, 
    std::vector<int>& srr_save, const Matrix& matrix) {

    std::vector<int> rows_to_check;

    // Loop through decomposed snarl
    for (const auto& snarl : decomposed_snarl) {
        if (snarl == "*") {
            continue;  // Skip this snarl if it's "*"
        }

        auto it = row_headers_dict.find(snarl);
        if (it != row_headers_dict.end()) {
            rows_to_check.push_back(it->second);  // Add row index if found
        } else {
            return {};  // Return empty vector if snarl is not found in the dictionary
        }
    }

    // Extract the rows to check
    std::vector<std::vector<bool>> extracted_rows;
    for (int row : rows_to_check) {
        extracted_rows.push_back(matrix.get_matrix()[row]);
    }

    // Find columns where all values are 1 (or true if row is std::vector<bool>)
    std::vector<int> idx_srr_save;
    if (!extracted_rows.empty()) {
        int num_cols = extracted_rows[0].size();
        
        // Ensure all rows have the same number of columns
        for (const auto& row : extracted_rows) {
            if (row.size() != num_cols) {
                // Handle case where row sizes are inconsistent (optional)
                // You might want to throw an error or handle it in a different way
                throw std::runtime_error("Inconsistent row size in extracted_rows");
            }
        }

        // Check each column
        for (int col = 0; col < num_cols; ++col) {
            bool all_ones = true;
            
            for (const auto& row : extracted_rows) {
                // Compare to true (if rows are std::vector<bool>)
                if (row[col] != true) {  // Use 'true' for boolean comparison
                    all_ones = false;
                    break;
                }
            }

            if (all_ones) {
                idx_srr_save.push_back(col);  // Store column index where all values are true (or 1)
            }
        }
    }

    // Assign to srr_save and return the result
    srr_save = idx_srr_save;
    return srr_save;
}

// Binary Table Generation
void SnarlParser::binary_table(const std::unordered_map<std::string, std::vector<std::string>>& snarls,
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
        std::vector<std::vector<int>> df = create_binary_table(binary_groups, list_snarl, list_samples, matrix);
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
void SnarlParser::quantitative_table(const std::unordered_map<std::string, std::vector<std::string>>& snarls,
                                        const std::unordered_map<std::string, float>& quantitative_phenotype,
                                        const std::string& output) 
{
    std::ofstream outf(output, std::ios::binary);
    if (!outf.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    // Write headers
    std::string headers = "CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tSE\tBETA\tP\n";
    outf.write(headers.c_str(), headers.size());

    // Iterate over each snarl
    for (const auto& [snarl, list_snarl] : snarls) {
        std::unordered_map<std::string, std::vector<int>> df = create_quantitative_table(list_snarl, list_samples, matrix);
        
        // std::make_tuple(se, beta, p_value)
        std::tuple<double, double, double> tuple_info = linear_regression(df, quantitative_phenotype);

        std::string chrom = "NA", pos = "NA", type_var = "NA", ref = "NA", alt = "NA";
        std::stringstream data;
        data << chrom << "\t" << pos << "\t" << snarl << "\t" << type_var << "\t" << ref << "\t" << alt
            << "\t" << std::get<0>(tuple_info) << "\t" << std::get<1>(tuple_info) 
            << "\t" << std::get<2>(tuple_info) << "\n";

        outf.write(data.str().c_str(), data.str().size());
    }
}
