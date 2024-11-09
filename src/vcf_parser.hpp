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

    Variant parseVariant(const std::string& line);
};
