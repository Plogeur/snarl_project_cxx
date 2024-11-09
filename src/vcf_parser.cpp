#include "vcf_parser.hpp"

// Constructor that opens the VCF file
VCFParser::VCFParser(const std::string& filename) : file(filename) {
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open VCF file: " + filename);
    }
}

// Destructor that closes the file stream
VCFParser::~VCFParser() {
    if (file.is_open()) {
        file.close();
    }
}

// Constructor for the begin iterator
VCFParser::Iterator::Iterator(std::ifstream* file) : file(file) {
    ++(*this);  // Load the first line
}

// Default constructor for the end iterator
VCFParser::Iterator::Iterator() : file(nullptr) {}


// Prefix increment to load the next variant
VCFParser::Iterator& VCFParser::Iterator::operator++() {
    loadNextVariant();
    return *this;
}

// Dereference operator to access the current variant
const Variant& VCFParser::Iterator::operator*() const {
    return currentVariant;
}

// Comparison operator for end detection
bool VCFParser::Iterator::operator!=(const Iterator& other) const {
    return file != other.file;
}

// Load the next variant from the file
void VCFParser::Iterator::loadNextVariant() {
    std::string line;
    if (file && std::getline(*file, line)) {
        currentVariant = VCFParser::parseVariant(line);  // Parse the line into a Variant object
    } else {
        file = nullptr;  // End of file reached
    }
}

// Iterator access methods for range-based for loop
VCFParser::Iterator VCFParser::begin() {
    return Iterator(&file);
}

VCFParser::Iterator VCFParser::end() {
    return Iterator();
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
    file.clear();                      // Clear any flags in the file stream
    file.seekg(0, std::ios::beg);      // Move the file stream to the beginning

    std::string line;
    while (std::getline(file, line)) {

        if (line.empty()) {
            continue;                  // Skip any empty lines
        }

        // Stop reading the header once we encounter a non-comment line
        if (line[0] != '#') {
            file.seekg(-line.size() - 1, std::ios::cur);  // Go back to the start of this line
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
}

// Parses a variant line from the VCF file and extracts genotype and AT field
Variant VCFParser::parseVariant(const std::string& line) {
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

    // Clear previous genotype data
    genotypes.clear();

    // Read the genotype data for each sample
    for (const auto& sample : sampleNames) {
        if (!(variantStream >> genotypeStr)) {
            throw std::runtime_error("Error reading genotype for sample: " + sample);
        }
        genotypes.push_back(extractGenotype(genotypeStr));
    }

    // Extract AT field from the INFO column
    atInfo = extractATField(infoField);
    
    return Variant(atInfo, genotypes)
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
