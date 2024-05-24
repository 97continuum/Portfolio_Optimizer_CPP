//
// Created by Talha Jamal on 20/05/2024.
//

#include "portfolio.h"
#include <cmath>
#include <numeric>
#include <iostream>
#include "linearAlgebra.h"
#include "math.h"

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
    std::vector<double> x0(numAssets, 1/numAssets); // Dimensions of Vector X
    x0.push_back(0.005);
    x0.push_back(0.005);

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
    //printMatrix(Q, "Q");
    // Fill b vector
    b[numAssets] = -targetReturn;
    b[numAssets+1] = -1.0;
    //printVector(b, "b");

    // Solve for Qx = b by Conjugate Method
    return conjugateGradient(Q, b, x0);
}

std::vector<double> Portfolio::conjugateGradient(const std::vector<std::vector<double>> &Q,
                                                 const std::vector<double> &b,
                                                 const std::vector<double> &x0)
{
    std::vector<double>  x, s, p, Qp, alphaQp;
    double alpha, beta, sTs, pTQp, sTsNew;
    size_t n = b.size();
    x = x0;
    s = vectorSubtraction(b, matrixVectorMultiplication(Q, x0)); // s: b - Qx
    p = s; // set Initial Direction
    sTs = vectorDotProduct(s, s);
    //printVector(x, "x");
    //printVector(s, "s");
    //printVector(p, "p");
    //std::cout << "sTs : " << sTs << std::endl;

    for (size_t k = 0; k < n; ++k)
    {
        Qp = matrixVectorMultiplication(Q, p);
        //printVector(Qp, "Q*p");
        pTQp = vectorDotProduct(p, Qp);
        //std::cout << "Dot Product " << pTQp << std::endl;
        if (std::abs(pTQp) < 1e-10) {
            std::cerr << "Division by nearly zero in alpha computation, terminating algorithm." << std::endl;
            break;
        }
        alpha = sTs / pTQp; // step size
        //std::cout << "alpha " << alpha << std::endl;
        x = vectorAddition(x, scalarMultiplication(p, alpha));
        //printVector(x, "x");
        alphaQp = scalarMultiplication(Qp, alpha);
        s = vectorSubtraction(s, alphaQp);
        //printVector(s, "s");
        sTsNew = vectorDotProduct(s, s);
        //std::cout << "sTsNew " << sTsNew << std::endl;

        if (sTsNew < 1.0E-6) {
            //std::cout << "sTsNew less than epsilon threshold - breaking from code" << std::endl;
            break;
        }
        beta = sTsNew / sTs;
        p = vectorAddition(s, scalarMultiplication(p, beta));
        //printVector(p, "p");
        sTs = sTsNew;
        //std::cout << "sTs : " << sTs << std::endl;
    }

    //std::cout << "Finished Optimization" << std::endl;
    return x;
}
