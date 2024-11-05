#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>

class VCFParser {
public:
    VCFParser(const std::string& filename);
    ~VCFParser();

    bool hasNext();
    void nextVariant();
    std::vector<std::string> getATInfo() const;
    const std::vector<std::vector<int> >& getGenotypes() const;
    const std::vector<std::string>& getSampleNames() const;

private:
    std::ifstream file;
    std::vector<std::string> sampleNames;
    std::vector<std::string> atInfo;
    std::vector<std::vector<int> > genotypes;

    void parseHeader();
    void parseVariant(const std::string& line);
    std::vector<int> extractGenotype(const std::string& genotypeStr);
    std::vector<std::string> extractATField(const std::string& infoField);
    std::vector<std::string> split(const std::string& s, char delimiter);
};