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
#include "parameter_estimation.h"
#include "linearAlgebra.h"
#include "unitTests.h"
#include "portfolio.h"

int  main (int  argc, char  *argv[])
{
    std::string desiredDirectory = "/Users/talhajamal/Desktop/Code/Portfolio_Optimizer_CPP"; // Home Directory
    changeWorkingDirectory(desiredDirectory); // Change to Correct Working Directory
    // Initialize Variables
    int numberAssets = 83; // Initialize Number of Assets
    int numberReturns = 700; // Max Length of Returns Data
    double **returnMatrix = new double*[numberAssets]; // a matrix to store the return data
    //allocate memory for return data
    for(int i=0;i<numberAssets;i++)
        returnMatrix[i]=new double[numberReturns];

    cout << "Reading Data" << std::endl;

    //read the data from the file and store it into the return matrix
    std::string fileName = "data/asset_returns.csv";
    checkFileInCurrentDirectory(fileName); // Check if File Exists and File Path is correct

    try {
        // Attempt to read the file
        std::cout << "Attempting to Open the File to read through it" << std::endl;
        readData(returnMatrix,fileName); // Read return data from the file and store in 2D returnMatrix

    } catch (const std::exception& e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }

    std::cout << "Converting to returnMatrix to vector" << std::endl;
    // Convert to vector of vectors
    std::vector<std::vector<double>> returns = convertToVectorMatrix(returnMatrix, numberAssets, numberReturns);
    std::cout << "Converted to vector successfully" << std::endl;
    // Test Linear Algebra Functions
    //testAllFunctions();

    // Instantiate Portfolio Object
    double targetReturns = 0.10;
    Portfolio portfolio(returns, targetReturns, 10, numberReturns);

    std::cout << "Checking if both meanReturn and covariance matrix function are calculated correctly" << std::endl;

    std::cout << "Parameter Estimation Script:" << std::endl;
    std::cout << " =========================== " << std::endl;
    std::cout << "Mean Return Calculation for First 10 Assets" << std::endl;
    std::vector<double> meanReturns = calculateMean(returnMatrix, 10, numberReturns);
    std::cout << "Mean Returns: " << std::endl;
    for (double mean: meanReturns)
        {
            std::cout << mean << " ";
        }
    std::cout << std::endl;

    std::cout << "Portfolio Class Method:" << std::endl;
    std::cout << " =========================== " << std::endl;
    std::cout << "Mean Return Calculation for First 10 Assets" << std::endl;

    // Calculate mean returns
    std::vector<double> portfolioMeanReturns = portfolio.calculateMeanReturn();
    std::cout << "Portfolio Mean Returns:\n";
    for (double r : portfolioMeanReturns) {
        std::cout << r << " ";
    }
    std::cout << std::endl;


    std::cout << "Parameter Estimation Script:" << std::endl;
    std::cout << " =========================== " << std::endl;
    std::cout << "Covariance Matrix Calculation for First 10 Assets" << std::endl;
    std::vector< std::vector<double> > covarianceMatrix = calculateCovarianceMatrix(returnMatrix, 10, numberReturns);
    std::cout << "CovarianceMatrix" << std::endl;
    for (const auto& row: covarianceMatrix)
    {
        for(const auto& value: row)
        {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Portfolio Class Method:" << std::endl;
    std::cout << " =========================== " << std::endl;
    std::cout << "Covariance Matrix Calculation for First 10 Assets" << std::endl;
    std::vector< std::vector<double> > portfolioCovMatrix = portfolio.calculateCovarianceMatrix();
    std::cout << "Covariance Matrix:" << std::endl;
    for (const auto& row: portfolioCovMatrix)
    {
        for (const auto& value: row)
        {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    // Delete Memory from Double Pointer
    for(int i=0;i<numberAssets;i++)
        delete[] returnMatrix[i];
    delete[] returnMatrix;
    cout << "Deleted Memory" << std::endl;
    return 0;
}
