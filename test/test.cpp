#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <tuple>
#include <numeric>

double mean(const std::vector<double>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double variance(const std::vector<double>& v) {
    double m = mean(v);
    double var = 0.0;
    for (double val : v) {
        var += (val - m) * (val - m);
    }
    return var / (v.size() - 1);
}

double covariance(const std::vector<double>& x, const std::vector<double>& y) {
    double mean_x = mean(x);
    double mean_y = mean(y);
    double cov = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        cov += (x[i] - mean_x) * (y[i] - mean_y);
    }
    return cov / (x.size() - 1);
}

std::tuple<double, double, double> linear_regression(const std::vector<std::vector<int>>& df, const std::unordered_map<std::string, float>& args) {
    if (df.size() < 2 || df[0].size() != 2) {
        throw std::invalid_argument("Input data should be a non-empty 2D vector with 2 columns (X, Y).");
    }

    // Separate X and Y from the dataset
    std::vector<double> X, Y;
    for (const auto& row : df) {
        X.push_back(static_cast<double>(row[0]));
        Y.push_back(static_cast<double>(row[1]));
    }

    // Calculate beta (slope) and intercept (alpha) for Y = alpha + beta * X
    double beta = covariance(X, Y) / variance(X);
    double alpha = mean(Y) - beta * mean(X);

    // Calculate standard error (SE) of beta
    double sum_errors_squared = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        double predicted_y = alpha + beta * X[i];
        double error = Y[i] - predicted_y;
        sum_errors_squared += error * error;
    }
    double se = std::sqrt(sum_errors_squared / (X.size() - 2)) / std::sqrt(variance(X) * (X.size() - 1));

    // Calculate t-statistic for beta and two-tailed p-value
    double t_stat = beta / se;
    double p_value = 2 * (1.0 - std::erf(std::abs(t_stat) / std::sqrt(2))); // Approximation using erf

    return std::make_tuple(p_value, se, beta);
}

int main() {
    // Sample data
    std::vector<std::vector<int>> df = {
        {1, 2},
        {2, 4},
        {3, 5},
        {4, 4},
        {5, 5}
    };

    std::unordered_map<std::string, float> args; // Not used in this example

    // Perform linear regression
    auto [p_value, se, beta] = linear_regression(df, args);

    // Output results
    std::cout << "P-value: " << p_value << std::endl;
    std::cout << "Standard Error: " << se << std::endl;
    std::cout << "Beta: " << beta << std::endl;

    return 0;
}
