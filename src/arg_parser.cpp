#include <arg_parser.hpp>

// Function to parse the binary phenotype file
std::unordered_map<std::string, bool> parse_group_file(const std::string& group_file) {
    std::unordered_map<std::string, bool> group;

    std::ifstream file(group_file);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << group_file << std::endl;
    }

    std::string line;
    std::getline(file, line);  // Skip header line

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string sample;
        int group_line;

        std::getline(iss, sample, '\t');
        iss >> group_line;

        if (group_line == 0) {
            group[sample] = 0;
        } else if (group_line == 1) {
            group[sample] = 1;
        }
    }

    file.close();
    return group;
}
// Function to parse the phenotype file
std::unordered_map<std::string, float> parse_pheno_file(const std::string& file_path) {
    std::unordered_map<std::string, float> parsed_pheno;

    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << file_path << std::endl;
        return parsed_pheno;  // Return empty map if file can't be opened
    }

    std::string line;
    std::getline(file, line);  // Skip header line

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string col1, iid;
        float pheno;

        // Read the first column, IID (second column), and PHENO (third column)
        std::getline(iss, col1, '\t');  // Skip the first column
        std::getline(iss, iid, '\t');   // Get the IID
        iss >> pheno;                   // Read the PHENO as a float

        parsed_pheno[iid] = pheno;  // Add to map
    }

    file.close();
    return parsed_pheno;
}

// Function to parse the snarl path file
std::unordered_map<std::string, std::vector<std::string>> parse_snarl_path_file(const std::string& path_file) {
    std::unordered_map<std::string, std::vector<std::string>> snarl_paths;
    std::ifstream file(path_file);

    if (!file.is_open()) {
        std::cerr << "Could not open the file: " << path_file << std::endl;
        return snarl_paths;
    }

    std::string line, snarl, path_list;
    
    // Read header
    std::getline(file, line);

    // Process each line
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::getline(ss, snarl, '\t');   // snarl column
        std::getline(ss, path_list, '\t'); // paths column
        
        if (!path_list.empty()) {
            std::istringstream path_stream(path_list);
            std::string path;
            while (std::getline(path_stream, path, ',')) {
                snarl_paths[snarl].push_back(path);
            }
        }
    }
    file.close();

    return snarl_paths;
}