//
// Created by Talha Jamal on 19/05/2024.
//
#include <iostream>
#include <fstream>
#include "read_data.h"
#include "portfolio.h"

using namespace std;

int  main (int  argc, char  *argv[])
{
    int numberAssets = 83; // Initialize Number of Assets
    int numberReturns = 700; // Max Length of Returns Data

    // Set list of increasing target Portfolio Returns
    const int steps = 20;
    double temp[steps];
    double start = 0.005;
    double end = 0.1;
    double increment = (end - start) / (steps - 1);

    // Backtesting Parameters
    double epsilon = 1e-6;
    int isWindow = 100;
    int oosWindow = 12;
    int slidingWindow = 12;
    int numOfSlidingWindows = 1 + ((numberReturns - isWindow - oosWindow) / (slidingWindow));

    vector<Portfolios> TargetReturnsPortfolios;
    string desiredDirectory = "/Users/talhajamal/Desktop/Code/Portfolio_Optimizer_CPP"; // Home Directory
    changeWorkingDirectory(desiredDirectory); // Change to Correct Working Directory
    string fileName = "data/asset_returns.csv"; // Return Data FileName
    checkFileInCurrentDirectory(fileName); // Check if File Exists and File Path is correct

    Matrix returnMatrix;
    returnMatrix.resize(numberReturns);
    for(int i=0;i<numberReturns;i++)
        returnMatrix[i].resize(numberAssets);
    readData(returnMatrix,fileName); // Read return data from the file and store in 2D returnMatrix

    for (int i = 0; i < steps; ++i) {temp[i] = start + i * increment;}
    Vector tReturns(temp, temp + steps);

    ofstream returnsFile; // Create CSV File to send returns to
    ofstream covFile; // Create CSV File to send covariances to
    returnsFile.open("data/results/PortfolioReturns.csv");
    covFile.open("data/results/Covariances.csv");
    returnsFile << "Target Returns,";
    covFile << "Target Returns,";

    for (int i = 0; i < numOfSlidingWindows - 1; ++i)
    {
        returnsFile << "Backtest Period: " << i+1 << ",";
        covFile << "Backtest Period: " << i+1 << ",";
    }
    returnsFile << "Backtest Period: " << numOfSlidingWindows << endl;
    covFile << "Backtest Period: " << numOfSlidingWindows << endl;

    for (int i = 0; i < 20; i++)
    {
        // Run backtest for current target return
        cout << "=================================================================" << endl;
        cout << "Running a Backtest for Portfolio with Target Returns of : " << tReturns[i] << endl;
        Portfolios portfolio(isWindow, oosWindow, slidingWindow, numOfSlidingWindows, returnMatrix, tReturns[i]);
        portfolio.calculateIsMean();
        portfolio.calculateIsCovMat();
        portfolio.calculateOOSMean();
        portfolio.calculateOOSCovMatrix();
        portfolio.calculateQ();
        portfolio.optimizer(epsilon);
        portfolio.runBacktest();
        TargetReturnsPortfolios.push_back(portfolio);

        // Results
        Vector actualAvgReturn = portfolio.getVectorOfActualAverageReturn();
        Vector actualCovMat = portfolio.getVectorOfActualCovMatrix();
        cout << "Portfolio's Actual Average Return is: " << portfolio.getAvgReturnPerBacktest() << endl;
        cout << "Portfolio's Actual Average Covariance is: " << portfolio.getAvgCovPerBacktest()<< endl;
        cout << "Portfolio's Sharpe Ratio is: " << portfolio.getPortfolioSharpeRatio() << endl;

        // Push Results to CSV
        returnsFile << tReturns[i] << ",";
        covFile << tReturns[i] << ",";
        for (int j = 0; j < actualAvgReturn.size()-1; j++)
        {
            returnsFile << actualAvgReturn[j] << ",";
        }
        returnsFile << actualAvgReturn[actualAvgReturn.size()-1] << endl;
        for (int k = 0; k < actualCovMat.size()-1; k++)
        {
            covFile << actualCovMat[k] << ",";
        }
        covFile << actualCovMat[actualCovMat.size()-1] << endl;
    }
    return 0;
}
