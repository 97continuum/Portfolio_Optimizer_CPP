//
// Created by Talha Jamal on 19/05/2024.
//
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include "csv.h"
#include "read_data.h"
#include "linearAlgebra.h"
#include "unitTests.h"
#include "portfolio.h"

int  main (int  argc, char  *argv[])
{
    std::string desiredDirectory = "/Users/talhajamal/Desktop/Code/Portfolio_Optimizer_CPP"; // Home Directory
    changeWorkingDirectory(desiredDirectory); // Change to Correct Working Directory

    // Initialize Variables
    int numberAssets = 5; // Initialize Number of Assets
    int numberReturns = 10; // Max Length of Returns Data
    auto **returnMatrix = new double*[numberAssets]; // a matrix to store the return data
    //allocate memory for return data
    for(int i=0;i<numberAssets;i++)
        returnMatrix[i]=new double[numberReturns];

    //read the data from the file and store it into the return matrix
    std::string fileName = "data/asset_returns_small.csv";
    checkFileInCurrentDirectory(fileName); // Check if File Exists and File Path is correct
    readData(returnMatrix,fileName); // Read return data from the file and store in 2D returnMatrix
    // Convert to vector of vectors
    std::vector<std::vector<double>> returns = convertToVectorMatrix(returnMatrix, numberAssets, numberReturns);
    //testAllFunctions(); // Test Linear Algebra Functions

    // Instantiate Portfolio Object
    double targetReturns = 0.001;
    Portfolio portfolio(returns, targetReturns, numberAssets, numberReturns);

    std::vector<double> portfolioMeanReturns = portfolio.calculateMeanReturn(); // Calculate mean returns
    printVector(portfolioMeanReturns, "Mean Returns");
    std::vector< std::vector<double> > portfolioCovMatrix = portfolio.calculateCovarianceMatrix(); // Calculate Cov Matrix
    printMatrix(portfolioCovMatrix, "Cov Matrix"); // Print Covariance Matrix

    std::vector<double> weights = portfolio.solveOptimization();
    printVector(weights, "Weights");

    // Delete Memory from Double Pointer
    deleteDoublePointer(returnMatrix, numberAssets);
    return 0;
}
