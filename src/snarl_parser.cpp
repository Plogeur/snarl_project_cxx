#include "snarl_parser.hpp"
#include "matrix.hpp"
#include "binary_analysis.hpp"
#include "quantitative_analysis.hpp"

SnarlParser::SnarlParser(const std::string& vcf_path) : filename(vcf_path), file(vcf_path), matrix(1000000, parseHeader().size() * 2) {
    sampleNames = parseHeader();
}

std::vector<std::string> SnarlParser::parseHeader() {
    file.clear();                     // Clear any flags in the file stream
    file.seekg(0, std::ios::beg);     // Move the file stream to the beginning

    std::vector<std::string> sampleNames;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;                 // Skip any empty lines
        }

        // Stop reading the header once we encounter a non-comment line
        if (line[0] != '#') {
            file.seekg(-static_cast<int>(line.size()) - 1, std::ios::cur);  // Go back to the start of this line
            break;
        }

        // Process the header line that starts with "#CHROM"
        if (line.substr(0, 6) == "#CHROM") {
            std::istringstream headerStream(line);
            std::string sampleName;

            // Skip the mandatory VCF columns
            for (int i = 0; i < 9; ++i) {
                headerStream >> sampleName;
            }

            // Read remaining entries as sample names
            while (headerStream >> sampleName) {
                sampleNames.push_back(sampleName);
            }
        }
    }
    return sampleNames;
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
const std::vector<std::vector<std::string>> decompose_snarl(const std::vector<std::string>& lst) {
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
    size_t currentRowsNumber = matrix.getRows();
    if (lengthOrderedMap > currentRowsNumber - 1) {
        matrix.expandMatrix();
    }

    // Add data to the matrix
    matrix.set(idxSnarl, indexColumn);
}

// Main function that parses the VCF file and fills the matrix
void SnarlParser::fill_matrix() {

    // Track overall start time
    auto start_total = std::chrono::high_resolution_clock::now();
    std::unordered_map<std::string, size_t> row_header_dict;

    for (std::string line; std::getline(file, line);) {

        Variant variant(line, sampleNames);
        const std::vector<std::vector<int>>& genotypes = variant.genotypes;
        const std::vector<std::string>& snarl_list = variant.atInfo;
        const std::vector<std::vector<std::string>> list_list_decomposed_snarl = decompose_snarl(snarl_list);

        // Process genotypes and fill matrix
        auto start_processing = std::chrono::high_resolution_clock::now();
        for (size_t index_column = 0; index_column < genotypes.size(); ++index_column) {
            int allele_1 = genotypes[index_column][0];
            int allele_2 = genotypes[index_column][1];
            size_t col_idx = index_column * 2;

            if (allele_1 == -1 || allele_2 == -1) { // Handle missing genotypes (./.)
                continue;
            }

            // Measure the time taken to push each allele into the matrix
            auto start_allele = std::chrono::high_resolution_clock::now();
            for (auto& decompose_allele_1 : list_list_decomposed_snarl[allele_1]) {
                pushMatrix(decompose_allele_1, row_header_dict, col_idx);
            }

            for (auto& decompose_allele_1 : list_list_decomposed_snarl[allele_2]) {
                pushMatrix(decompose_allele_1, row_header_dict, col_idx);
            }

            auto end_allele = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_allele = end_allele - start_allele;
            std::cout << "Allele processing time for column " << index_column << ": " << elapsed_allele.count() << " seconds\n";
        }

        auto end_processing = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_processing = end_processing - start_processing;
        std::cout << "Genotype processing time: " << elapsed_processing.count() << " seconds\n";

        // Time setting row headers
        auto start_row_header = std::chrono::high_resolution_clock::now();
        matrix.set_row_header(row_header_dict);
        auto end_row_header = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_row_header = end_row_header - start_row_header;
        std::cout << "Row header setting time: " << elapsed_row_header.count() << " seconds\n";
        std::cout << "---------------------------------------------------------------------------" << std::endl;
        std::cout << "" << std::endl;
    }

    // Track overall end time
    auto end_total = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_total = end_total - start_total;
    std::cout << "Total fill_matrix execution time: " << elapsed_total.count() << " seconds\n";

    std::cout << "Matrix : " << std::endl;
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 8; ++j) {
            std::cout << matrix(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

std::vector<int> identify_correct_path(
    const std::vector<std::string>& decomposed_snarl, 
    const std::unordered_map<std::string, size_t>& row_headers_dict, 
    std::vector<int>& srr_save, const Matrix& matrix, const size_t num_cols) {
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

    std::cout << "identify_correct_path 1" << std::endl;
    std::vector<std::vector<bool>> extracted_rows(rows_to_check.size());

    for (size_t i = 0; i < rows_to_check.size(); i++) {
        int row = rows_to_check[i];
        
        // Loop over the columns in the row and add each column's data to extracted_rows[i]
        for (size_t j = row; j < row + num_cols; ++j) {
            extracted_rows[i].push_back(matrix.get_matrix()[j]);
        }
    }
    std::cout << "identify_correct_path 2" << std::endl;

    // Find columns where all values are 1 (or true if row is std::vector<bool>)
    std::vector<int> idx_srr_save;
    if (!extracted_rows.empty()) {
        
        // Ensure all rows have the same number of columns
        for (const auto& row : extracted_rows) {
            if (row.size() != num_cols) {
                // Handle case where row sizes are inconsistent (optional)
                // You might want to throw an error or handle it in a different way
                throw std::runtime_error("Inconsistent row size in extracted_rows");
            }
        }
        std::cout << "identify_correct_path 3" << std::endl;

        // Check each column
        for (size_t col = 0; col < num_cols; ++col) {
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
    
    std::cout << "identify_correct_path 4" << std::endl;

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
        std::vector<std::vector<int>> df = create_binary_table(binary_groups, list_snarl, sampleNames, matrix);
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
        std::cout << "machin" << std::endl;

        std::unordered_map<std::string, std::vector<int>> df = create_quantitative_table(list_snarl, sampleNames, matrix);
        std::cout << "chose" << std::endl;

        // std::make_tuple(se, beta, p_value)
        std::tuple<double, double, double> tuple_info = linear_regression(df, quantitative_phenotype);

        std::cout << "bidule" << std::endl;

        std::string chrom = "NA", pos = "NA", type_var = "NA", ref = "NA", alt = "NA";
        std::stringstream data;
        data << chrom << "\t" << pos << "\t" << snarl << "\t" << type_var << "\t" << ref << "\t" << alt
            << "\t" << std::get<0>(tuple_info) << "\t" << std::get<1>(tuple_info) 
            << "\t" << std::get<2>(tuple_info) << "\n";

        outf.write(data.str().c_str(), data.str().size());
    }
}

// Parses a variant line from the VCF file and extracts genotype and AT field
Variant::Variant(std::string& line, std::vector<std::string> sampleNames) {
    std::istringstream variantStream(line);
    std::string columnData;
    std::string infoField;
    std::string genotypeStr;
    // Skip the first 9 columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
    for (int i = 0; i < 9; ++i) {
        variantStream >> columnData;
        if (i == 7) {
            infoField = columnData; // Capture the INFO column for parsing
        }
    }
    // Read the genotype data for each sample
    for (const auto& sample : sampleNames) {
        if (!(variantStream >> genotypeStr)) {
            throw std::runtime_error("Error reading genotype for sample: " + sample);
        }
        genotypes.push_back(extractGenotype(genotypeStr));
    }
    // Extract AT field from the INFO column
    atInfo = extractATField(infoField);
}

// Extracts the genotype from the genotype string
std::vector<int> extractGenotype(const std::string& genotypeStr) {
    std::vector<int> alleles;
    std::vector<std::string> alleleStrs = split(genotypeStr, '/');

    for (const auto& alleleStr : alleleStrs) {
        if (alleleStr == ".") {
            alleles.push_back(-1); // Use -1 for missing data
        } else {
            alleles.push_back(std::stoi(alleleStr)); // Convert allele string to integer
        }
    }
    return alleles;
}

// Extracts the AT field from the INFO field and returns a vector of strings
std::vector<std::string> extractATField(const std::string& infoField) {
    std::vector<std::string> atValues; // Change to vector of strings
    std::vector<std::string> infoParts = split(infoField, ';');
    
    for (const std::string& part : infoParts) {
        if (part.rfind("AT=", 0) == 0) {
            std::string atData = part.substr(3);  // Remove "AT="
            // Split the AT data using ',' only
            atValues = split(atData, ','); // Directly assign the split result to atValues
            return atValues; // Return the result
        }
    }
    return atValues;  // Return empty if AT field is not found
}

// Utility function to split strings
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}
