#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>

class Variant {
public:
    std::vector<std::string> atInfo;
    std::vector<std::vector<int>> genotypes;

    Variant(const std::vector<std::string>& atInfo, const std::vector<std::vector<int>>& genotypes)
        : atInfo(atInfo), genotypes(genotypes) {}
};

class VCFParser {
public:
    VCFParser(const std::string& filename);
    ~VCFParser();

   // Getter of list_samples
    const std::vector<std::string>& getSampleNames();

    static Variant parseVariant(const std::string& line);

    std::vector<std::string> extractATField(const std::string& infoField);

    std::vector<std::string> split(const std::string& s, char delimiter);

    std::vector<std::string> getATInfo() const;

    const std::vector<std::vector<int> >& getGenotypes() const;

    std::vector<int> extractGenotype(const std::string& genotypeStr);

    void parseHeader();

    // Iterator class for range-based for loop
    class Iterator {
    public:
        Iterator(std::ifstream* file);
        Iterator();  // End iterator

        // Prefix increment to load the next Variant
        Iterator& operator++();
        const Variant& operator*() const;
        bool operator!=(const Iterator& other) const;

    private:
        std::ifstream* file;
        Variant currentVariant;
        void loadNextVariant();
    };

    Iterator begin();
    Iterator end();

private:
    std::ifstream file;
    std::vector<std::string> sampleNames;
};
