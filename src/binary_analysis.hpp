#ifndef BINARY_ANALYSIS_HPP
#define BINARY_ANALYSIS_HPP

#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <iomanip>
#include <sstream>
#include "matrix.hpp"
#include "snarl_parser.hpp"

// ------------------------ Chi2 test ------------------------

// Function to calculate the Chi-square test statistic
double chiSquareStatistic(const std::vector<std::vector<int>>& observed);

// Function to calculate the degrees of freedom for a 2D table
int calculateDegreesOfFreedom(int rows, int cols);

// Function to compute the regularized incomplete gamma function (for Chi-square CDF approximation)
double gammaIncomplete(double s, double x);

// Function to calculate the p-value from Chi-square statistic using incomplete gamma function
double chiSquarePValue(double chiSquare, int degreesOfFreedom);

// Function to perform the Chi-square test
std::string chi2Test(const std::vector<std::vector<int>>& observed);

// ------------------------ Fisher exact test ------------------------

// Function to initialize the log factorials array
void initLogFacs(long double* logFacs, int n);

// Function to calculate the log probability of the hypergeometric distribution
long double logHypergeometricProb(long double* logFacs , int a, int b, int c, int d);

// Function to perform Fisher's exact test
long double fastFishersExactTest(const std::vector<std::vector<int>>& table);

// ------------------------ Binary table & stats ------------------------

std::vector<std::string> binary_stat_test(const std::vector<std::vector<int>>& df);

std::vector<std::vector<int>> create_binary_table(
    const std::unordered_map<std::string, bool>& groups, 
    const std::vector<std::string>& list_path_snarl, 
    const std::vector<std::string>& list_samples, Matrix& matrix);

#endif // CHI2FISHERTEST_HPP
