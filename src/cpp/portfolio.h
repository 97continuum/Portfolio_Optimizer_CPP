//
// Created by Talha Jamal on 20/05/2024.
//

#ifndef COURSEWORK_PORTFOLIO_H
#define COURSEWORK_PORTFOLIO_H

#include <iostream>

class Portfolio {
public:
    Portfolio(const std::vector< std::vector<double> >& returns, double targetReturn, int numAssets, int numPeriods);

    std::vector<double> calculateMeanReturn(); // Function to calculate
    std::vector< std::vector<double> > calculateCovarianceMatrix(); // function to calculate Covariance Matrix
    std::vector<double> solveOptimization(); // Function to create Qx = b systems of linear equations

private:
    std::vector< std::vector<double> > returns; // Matrix to Store Returns
    double targetReturn; // Variable to hold Target Return for Portfolio
    int numAssets; // Total Number of Assets
    int numPeriods; // Total Number of Periods
    std::vector<double> meanReturns; // The vector of Mean Returns for N Assets
    std::vector< std::vector<double> > covarianceMatrix; // Variance Covariance Matrix for N Assets
    std::vector<double> conjugateGradient(const std::vector<std::vector<double>> &Q, const std::vector<double> &b,
                                          const std::vector<double> &x0) const; // Function to solve Qx = b equation
};


#endif //COURSEWORK_PORTFOLIO_H
