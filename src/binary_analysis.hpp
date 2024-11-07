#ifndef BINARY_ANALYSIS_HPP
#define BINARY_ANALYSIS_HPP

#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <numeric>

// ------------------------ Chi2 exact test ------------------------

// Function to calculate the Chi-square test statistic
double chiSquareStatistic(const std::vector<std::vector<int>>& observed);

// Function to calculate the degrees of freedom for a 2D table
int calculateDegreesOfFreedom(int rows, int cols);

// Function to calculate the p-value based on the Chi-square statistic using an approximation
double chiSquarePValue(double chiSquare, int degreesOfFreedom);

// Function to perform the Chi-square test
std::string chi2Test(const std::vector<std::vector<int>>& observed);

// ------------------------ Fisher exact test ------------------------

// Function to initialize the log factorials array
void initLogFacs(double* logFacs, int n);

// Function to calculate the log probability of the hypergeometric distribution
double logHypergeometricProb(double* logFacs, int a, int b, int c, int d);

// Function to perform Fisher's exact test
double fastFishersExactTest(const std::vector<std::vector<int>>& table);

#endif // CHI2FISHERTEST_HPP
