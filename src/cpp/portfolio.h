//
// Created by Talha Jamal on 20/05/2024.
//

#ifndef COURSEWORK_PORTFOLIO_H
#define COURSEWORK_PORTFOLIO_H

#include <iostream>
#include "linearAlgebra.h"

using namespace std;

Matrix sliceMatrixByRows(const Matrix& matrix, int row_start, int row_end);
Vector calculateAverage(const Matrix& m);
Matrix calculateCovMatrix(const Matrix& m, Vector meanReturns);

class BacktestingEngine {
public:
    BacktestingEngine(int isWindow_,int oosWindow_,int slidingWindow_,int numOfSlidingWindows_,
                      const Matrix& AssetReturns_, double targetReturn_);
    // IS Functions
    void calculateIsMean();
    void calculateIsCovMat();
    void calculateIsQ();
    void optimizer(double epsilon);

    void setb(double targetReturn_)
    {
        targetReturn = targetReturn_;
        Vector b(numAssets, 0.0);
        b.push_back(-targetReturn);
        b.push_back(-1);
        isb = b;
    }

    // OOS Functions
    void calculateOOSMean();
    void calculateOOSCovMatrix();

    ~BacktestingEngine();

protected:
    int isWindow, oosWindow, slidingWindow, numOfSlidingWindows;;
    Matrix AssetReturns;
    size_t numAssets, numReturns;

    // OOS Variables
    Matrix oosMean;
    vector<Matrix> oosCovMatrix;

    // IS Variables - including variables needed for Portfolio Optimization
    double targetReturn;
    Vector isb, isLambda, isMu;
    Matrix isMean, isWeights, X;
    vector<Matrix> isCovMatrix, isQ;
};

class Portfolios : public BacktestingEngine {
public:
    Portfolios(int isWindow_, int oosWindow_, int slidingWindow_, int numOfSlidingWindows_,
               const Matrix& AssetReturns_, double targetReturn_);

    void runBacktest();
    Vector getActualAverageReturn();
    Vector getActualCovMatrix();
    double getAvgAbnormalReturn();
    double getCumulativeAvgAbnormalReturn();
    double getStd();
    double getPortfolioSharpeRatio() const;
protected:
    Vector actualAverageReturn;
    Vector actualCovMatrix;
    double avgAbnormalReturn;
    double cumulativeAvgAbnormalReturn;
    double standardDeviation;
    double portfolioSharpeRatio;
};

#endif //COURSEWORK_PORTFOLIO_H
