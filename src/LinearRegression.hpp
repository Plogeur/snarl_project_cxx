#ifndef LINEAR_REGRESSION_HPP
#define LINEAR_REGRESSION_HPP

#include <Eigen/Dense>
#include <vector>

// Function to perform linear regression and return beta coefficients, standard errors, and p-values
void linearRegression(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Eigen::VectorXd& beta, Eigen::VectorXd& se, Eigen::VectorXd& p_values);

#endif // LINEAR_REGRESSION_HPP
