//
// Created by Talha Jamal on 20/05/2024.
//

#include "portfolio.h"
#include <cmath>
#include <numeric>
#include <iostream>
#include "linearAlgebra.h"

// Constructor for Portfolio Class
Portfolio::Portfolio(const std::vector<std::vector<double>> &returns, double targetReturn, int numAssets, int numPeriods)
    : returns(returns), targetReturn(targetReturn), numAssets(numAssets), numPeriods(numPeriods){}

// Calculate Mean Returns
std::vector<double> Portfolio::calculateMeanReturn()
{
    //size_t numAssets = returns.size();
    //size_t numPeriods = returns[0].size();
    meanReturns.resize(numAssets);

    for (int i = 0; i < numAssets; ++i)
    {
        meanReturns[i] = std::accumulate(returns[i].begin(), returns[i].end(), 0.0) / numPeriods;
    }
    return meanReturns;
}

// Calculate Covariance Matrix
std::vector< std::vector<double> > Portfolio::calculateCovarianceMatrix() {
    //size_t numAssets = returns.size();
    //size_t numPeriods = returns[0].size();
    covarianceMatrix.resize(numAssets, std::vector<double>(numAssets, 0.0));

    // Calculate Mean Returns if they have not been calculated
    if (meanReturns.empty()) {
        calculateMeanReturn();
    }

    for (size_t i = 0; i < numAssets; ++i) {
        for (size_t j = 0; j < numAssets; ++j) {
            double cov = 0.0;
            for (size_t k = 0; k < numPeriods; ++k) {
                cov += (returns[i][k] - meanReturns[i]) * (returns[j][k] - meanReturns[j]);
            }
            covarianceMatrix[i][j] = cov / (numPeriods - 1);
        }
    }
    return covarianceMatrix;
}

// solve Optimization Problem by creating the System of Linear Equations needed for Conjugate Gradient Method
std::vector<double> Portfolio::solveOptimization()
{
    //size_t n = meanReturns.size(); // number of Assets
    std::vector< std::vector<double> > Q(numAssets + 2, std::vector<double>(numAssets + 2, 0.0)); // Dimensions of Matrix Q
    std::vector<double> b(numAssets + 2, 0.0); // Dimensions of Vector b
    std::vector<double> x0(numAssets + 2, 0.0); // Dimensions of Vector X

    // Fill Q Matrix
    for (size_t i=0; i < numAssets; ++i)
    {
        for (size_t j = 0; j < numAssets; ++j)
        {
            Q[i][j] = covarianceMatrix[i][j];
        }
        Q[i][numAssets] = -meanReturns[i];
        Q[i][numAssets+1] = -1.0;
        Q[numAssets][i] = -meanReturns[i];
        Q[numAssets+1][i] = -1.0;
    }
    printMatrix(Q, "Q");
    // Fill b vector
    b[numAssets] = -targetReturn;
    b[numAssets+1] = -1.0;
    printVector(b, "b");

    // Solve for Qx = b by Conjugate Method
    return conjugateGradient(Q, b, x0);
}

std::vector<double> Portfolio::conjugateGradient(const std::vector<std::vector<double>> &Q,
                                                 const std::vector<double> &b,
                                                 const std::vector<double> &x0) const
{
    //size_t n = b.size();
    std::vector<double> x = x0;
    printVector(x, "x");
    std::vector<double> s = vectorSubtraction(b, matrixVectorMultiplication(Q, x0)); // s: b - Qx
    printVector(s, "s");
    std::vector<double> p = s; // set Initial Direction
    printVector(p, "p");
    double sTs = vectorDotProduct(s, s);
    std::cout << "sTs : " << sTs << std::endl;

    for (size_t i=0; i < numAssets; ++i)
    {
        double alpha = sTs / (vectorDotProduct(p, matrixVectorMultiplication(Q, p))); // step size
        std::cout << "alpha" << alpha << std::endl;
        x = vectorAddition(x, scalarMultiplication(p, alpha));
        printVector(x, "x");
        s = vectorSubtraction(s, scalarMultiplication(p, alpha));
        printVector(s, "s");
        double sTsNew = vectorDotProduct(s, s);
        std::cout << "sTsNew " << sTsNew << std::endl;
        if (sTsNew < 1.0E-6)
        {
            std::cout << "sTsNew less than epsilon threshold - breaking out from code" << std::endl;
            break;
        }
        p = vectorAddition(s, scalarMultiplication(p, sTsNew/sTs));
        printVector(p, "p");
        sTs = sTsNew;
        std::cout << "sTs : " << sTs << std::endl;
    }
    std::cout << "Finished Optimization" << std::endl;
    return x;
}
