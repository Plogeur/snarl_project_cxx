#include "vcf_parser.hpp"

// Constructor: Opens the VCF file and parses the header
VCFParser::VCFParser(const std::string& filename) : file(filename) {
    file.open(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open VCF file.");
    }
    parseHeader();
}

// Destructor: Closes the VCF file if it is open
VCFParser::~VCFParser() {
    if (file.is_open()) {
        file.close();
    }
}

// Checks if there are more lines to read in the VCF file
bool VCFParser::hasNext() {
    return !file.eof();
}

// Reads the next variant from the VCF file
void VCFParser::nextVariant() {
    std::string line;
    if (std::getline(file, line)) {
        parseVariant(line);
    } else {
        throw std::runtime_error("Error reading next variant.");
    }
}

// Gets the AT information for the current variant
std::vector<std::string> VCFParser::getATInfo() const {
    return atInfo;
}

// Gets the genotypes for the current variant
const std::vector<std::vector<int> >& VCFParser::getGenotypes() const {
    return genotypes;
}

// Gets the sample names from the VCF header
const std::vector<std::string>& VCFParser::getSampleNames() const {
    return sampleNames;
}

// Parses the VCF header to extract sample names
void VCFParser::parseHeader() {

    file.clear();  // Ensure the file stream is ready
    file.seekg(0, std::ios::beg);  // Reset to the beginning of the file

    std::string line;
    while (std::getline(file, line)) {

        if (line.empty()) {
            continue;  // Skip any empty lines
        }

        if (line[0] != '#') {
            break;  // Stop when we hit the data lines
        }

        if (line.substr(0, 6) == "#CHROM") {
            std::istringstream headerStream(line);
            std::string sampleName;
            while (headerStream >> sampleName) {
                if (sampleName != "#CHROM" && sampleName != "POS" &&
                    sampleName != "ID" && sampleName != "REF" &&
                    sampleName != "ALT" && sampleName != "QUAL" &&
                    sampleName != "FILTER" && sampleName != "INFO" &&
                    sampleName != "FORMAT") {
                    sampleNames.push_back(sampleName);
                }
            }
        }
    }
}

// Parses a variant line from the VCF file
void VCFParser::parseVariant(const std::string& line) {
    std::istringstream variantStream(line);
    std::string infoField, formatField;
    std::string genotypeStr;

    // Skip the first 8 columns
    for (int i = 0; i < 8; ++i) {
        variantStream >> genotypeStr; // Just reading the columns, we can ignore them
    }

    // Read the FORMAT field
    variantStream >> formatField;

    // Clear previous genotype data
    genotypes.clear();

    // Now read the genotype data for each sample
    for (const auto& sample : sampleNames) {
        if (!(variantStream >> genotypeStr)) {
            throw std::runtime_error("Error reading genotype for sample: " + sample);
        }
        genotypes.push_back(extractGenotype(genotypeStr));
    }

    // Extract AT information from the INFO field
    if (variantStream >> infoField) {
        atInfo = extractATField(infoField);
    }
}

// Extracts the genotype from the genotype string
std::vector<int> VCFParser::extractGenotype(const std::string& genotypeStr) {
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
std::vector<std::string> VCFParser::extractATField(const std::string& infoField) {
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
std::vector<std::string> VCFParser::split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}
