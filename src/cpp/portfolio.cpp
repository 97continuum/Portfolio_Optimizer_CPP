#include <numeric>
#include "portfolio.h"
#include "linearAlgebra.h"
#include <stdexcept>
#include <cmath>

using namespace std;

// Function to calculate the mean of a vector
double calculateAverage(const std::vector<double>& avgReturns) {
    double sum = accumulate(avgReturns.begin(), avgReturns.end(), 0.0);
    return sum / avgReturns.size();
}

// Function to calculate the variance of a vector
double calculateVariance(const std::vector<double>& avgReturns) {
    double mean = calculateAverage(avgReturns);
    double variance = 0.0;

    for (const double& return_value : avgReturns) {
        variance += pow(return_value - mean, 2);
    }

    return variance / avgReturns.size();
}

Vector calculateAverage(const Matrix& m)
{
    if (m.empty()){return {};}
    size_t numRows = m.size();
    size_t numCols = m[0].size();
    Vector result(numCols, 0.0);
    // Sum each column's elements
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result[j] += m[i][j];
        }
    }
    // Divide by number of Rows to get the mean
    for (int i = 0; i < numCols; ++i)
    {
        result[i] /= numRows;
    }
    return result;
}

Matrix calculateCovMatrix(const Matrix& m, Vector meanReturns)
{
    size_t numCols = m[0].size();
    size_t numRows = m.size();
    Matrix covMatrix(numCols, Vector(numCols, 0.0));
    //Vector meanReturns = calculateAverage(m);

    for (size_t i = 0; i < numCols; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            double cov = 0.0;
            for (size_t k = 0; k < numRows; ++k) {
                cov += (m[k][i] - meanReturns[i]) * (m[k][j] - meanReturns[j]);
            }
            covMatrix[i][j] = cov / (numRows - 1);
        }
    }
    return covMatrix;
}

// Backtesting Class Constructor
BacktestingEngine::BacktestingEngine(int isWindow_,int oosWindow_,int slidingWindow_,int numOfSlidingWindows_,
                                     const Matrix& AssetReturns_, double targetReturn_):
                                     isWindow(isWindow_),oosWindow(oosWindow_),slidingWindow(slidingWindow_),
                                     numOfSlidingWindows(numOfSlidingWindows_),AssetReturns(AssetReturns_),
                                     targetReturn(targetReturn_)
{
    numReturns = AssetReturns.size();
    numAssets = AssetReturns[0].size();
    targetReturn = targetReturn_;
    Vector temp(numAssets, 0.0);
    temp.push_back(-targetReturn_);
    temp.push_back(-1);
    isb = temp;
}

Matrix BacktestingEngine::getISMean(){ return isMean;};
vector<Matrix> BacktestingEngine::getISCovMat() {return isCovMatrix;};
Matrix BacktestingEngine::getOSMean() { return oosMean;};
vector<Matrix> BacktestingEngine::getOSCovMat() { return oosCovMatrix;};
vector<Matrix> BacktestingEngine::getQ() { return Q;};
Matrix BacktestingEngine::getWeights() { return isWeights;};

// Destructor for Backtesting Engine Class
BacktestingEngine::~BacktestingEngine() = default;

void BacktestingEngine::calculateIsMean()
{
    Matrix M; // Matrix to Store Returns for each IS Period
    Vector V; // Store Mean Return for each IS Period -> Mean Return Vectors for each window pushed to isMean Matrix
    for (int i = 0; i < numOfSlidingWindows; i++)
    {
        for (int j = slidingWindow * i; j < (slidingWindow * i) + isWindow; j++)
        {
            M.push_back(AssetReturns[j]);
        }
        V = calculateAverage(M);
        isMean.push_back(V);
        M.clear();
    }
}

void BacktestingEngine::calculateIsCovMat()
{
    Matrix M1; // Matrix to Store Returns for each OOS Period
    Matrix M2; // Matrix to store Cov Matrix for each IS Period -> Each IS Window's CovMat pushed to oosCovMatrix
    for (int i = 0; i < numOfSlidingWindows; i++)
    {
        for (int j = slidingWindow * i; j < (slidingWindow * i) + isWindow; j++)
        {
            M1.push_back(AssetReturns[j]);
        }
        M2 = calculateCovMatrix(M1, isMean[i]);
        isCovMatrix.push_back(M2);
        M1.clear();
    }
}

void BacktestingEngine::calculateOOSMean()
{
    Matrix M; // Matrix to Store Returns for each OOS Period
    Vector V; // Store Mean Return for each OOS Period -> Mean Return Vectors for each window pushed to oosMean Matrix
    for (int i = 0; i < numOfSlidingWindows; i++)
    {
        for (int j = isWindow + (slidingWindow * i); j < isWindow + (slidingWindow * i) + oosWindow; j++)
        {
            M.push_back(AssetReturns[j]);
        }
        V = calculateAverage(M);
        oosMean.push_back(V);
        M.clear();
    }
}

void BacktestingEngine::calculateOOSCovMatrix()
{
    Matrix M1; // Matrix to Store Returns for each OOS Period
    Matrix M2; // Matrix to store Cov Matrix for each IS Period -> Each IS Window's CovMat pushed to oosCovMatrix
    for (int i = 0; i < numOfSlidingWindows; i++)
    {
        for (int j = isWindow + (slidingWindow * i); j < isWindow + (slidingWindow * i) + oosWindow; j++)
        {
            M1.push_back(AssetReturns[j]);
        }
        M2 = calculateCovMatrix(M1, oosMean[i]);
        oosCovMatrix.push_back((M2));
        M1.clear();
    }
}

// solve Optimization Problem by creating the System of Linear Equations needed for Conjugate Gradient Method
void BacktestingEngine::calculateQ() {
    Matrix tempQ;
    Matrix tempCov;
    Vector tempMean;
    Vector tempOnes(isCovMatrix[0].size(), -1.0);

    tempOnes.push_back(0);
    tempOnes.push_back(0);

    for (int i = 0; i < numOfSlidingWindows; i++) {
        tempCov = isCovMatrix[i];
        tempMean = isMean[i] * (-1);
        for (int j = 0; j < isCovMatrix[0].size(); j++) {
            tempCov[j].push_back(tempMean[j]);
            tempCov[j].push_back(-1.0);
            tempQ.push_back(tempCov[j]);
        }

        tempMean.push_back(0);
        tempMean.push_back(0);

        tempQ.push_back(tempMean);
        tempQ.push_back(tempOnes);
        Q.push_back(tempQ);
        tempQ.clear();
    }
}

void BacktestingEngine::optimizer(double epsilon)
{
    Vector  s_k, s_k1, p_k, p_k1, Qp, alphaQp, weights;
    double alpha, beta, pTQp;

    for (int i = 0; i < numOfSlidingWindows; i++)
    {
        Vector x0(numAssets, 1/numAssets); // initializing vector x in Qx = b equation
        x0.push_back(0.1); // initializing value for lambda
        x0.push_back(0.1); // initializing value for mu

        //printVector(x0, "x0");

        s_k = isb - (Q[i]*x0);
        //printVector(s_k, "s_k");
        p_k = s_k;

        //double sTs = s_k * s_k;
        //cout << "sTs " << sTs << endl;

        while ((s_k * s_k) > epsilon)
        {
            alpha = (s_k * s_k)/(p_k * (Q[i] * p_k));
            //cout << "Alpha " << alpha << endl;

            x0 = x0 + (p_k * alpha);
            //printVector(x0, "x0");

            s_k1 = s_k - ((Q[i] * p_k) * alpha);
            //printVector(s_k1, "s_k1");

            beta = (s_k1 * s_k1)/(s_k * s_k);
            //cout << "Beta " << beta << endl;

            p_k1 = s_k1 + (p_k * beta);
            //printVector(p_k1, "p_k1");

            s_k = s_k1;
            p_k = p_k1;
        }
        //printVector(x0, "x0");
        X.push_back(x0);
        for (int j = 0; j < x0.size()-2; j++)
        {
            weights.push_back(x0[j]); // Add weights except Lagrange Multipliers from x
        }
        //printVector(weights, "weights");
        isWeights.push_back(weights);
        weights.clear();
        isLambda.push_back(x0[x0.size()-2]);
        isMu.push_back(x0[x0.size()-1]);
    }
}

// Constructor for Portfolios Class
Portfolios::Portfolios(int isWindow_, int oosWindow_, int slidingWindow_, int numOfSlidingWindows_,
                       const Matrix& AssetReturns_, double targetReturn_)
        : BacktestingEngine(isWindow_, oosWindow_, slidingWindow_, numOfSlidingWindows_, AssetReturns_, targetReturn_) {}

void Portfolios::runBacktest()
{
    avgReturn = 0.0;
    avgCovariance = 0.0;
    double variance = 0.0;
    for (int i = 0; i < numOfSlidingWindows; i++)
    {
        avgReturn = isWeights[i] * oosMean[i];
        vectorOfAverageReturn.push_back(avgReturn); // Actual Average Return for Each Backtest Period
        avgCovariance = isWeights[i] * (isCovMatrix[i] * isWeights[i]);
        vectorOfAverageCov.push_back(avgCovariance);// Actual Cov Matrix for each Backtest
    }
    // Calculate Average Return & Cov per backtest period
    avgReturnPerBacktest = calculateAverage(vectorOfAverageReturn);
    avgCovPerBacktest = calculateAverage(vectorOfAverageCov);
    standardDeviation = sqrt(calculateVariance(vectorOfAverageReturn)/numOfSlidingWindows); // calculate standard deviation
    portfolioSharpeRatio = (calculateAverage(vectorOfAverageReturn) - 0) / standardDeviation; // risk free rate = 0
}

// Get Methods
double Portfolios::getAvgReturnPerBacktest(){ return avgReturnPerBacktest;};
double Portfolios::getAvgCovPerBacktest(){ return avgCovPerBacktest;};
Vector Portfolios::getVectorOfActualAverageReturn(){ return vectorOfAverageReturn;}
Vector Portfolios::getVectorOfActualCovMatrix(){ return vectorOfAverageCov;};
double Portfolios::getStd(){ return standardDeviation;};
double Portfolios::getPortfolioSharpeRatio() const { return portfolioSharpeRatio;}
