#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "csv.h"  // Assuming you have a Csv class defined somewhere

double string_to_double(const std::string& s) {
    try {
        return std::stod(s);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << std::endl;
        return 0.0;
    } catch (const std::out_of_range& e) {
        std::cerr << "Out of range: " << e.what() << std::endl;
        return 0.0;
    }
}

void readData(double **data, const std::string& fileName) {
    //std::cout << "Attempting to open file: " << fileName << std::endl;

    std::ifstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << fileName << std::endl;
        return;
    }

    Csv csv(file);
    std::string line;
    int i = 0;
    while (csv.getline(line) != 0) {
        for (int j = 0; j < csv.getnfield(); j++) {
            double temp = string_to_double(csv.getfield(j));
            data[j][i] = temp;
        }
        i++;
    }
    file.close();
    //std::cout << "File read successfully: " << fileName << std::endl;
}

bool changeWorkingDirectory(const std::string& newDir) {
    if (chdir(newDir.c_str()) != 0)
    {
        perror("chdir() error");
        std::cerr << "Failed to change working directory to: " << newDir << std::endl;
        return false;
    }
    else
    {
        //std::cout << "Successfully changed working directory to: " << newDir << std::endl;
        return true;
    }
}

bool checkFileInCurrentDirectory(const std::string& fileName) {
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        std::cout << "Current working dir: " << cwd << std::endl;
    } else {
        perror("getcwd() error");
        return false;
    }

    std::string fullPath = std::string(cwd) + "/" + fileName;
    std::ifstream file(fullPath);
    if (!file) {
        std::cerr << "File " << fullPath << " does not exist." << std::endl;
        return false;
    } else {
        //std::cout << "File " << fullPath << " exists" << std::endl;
    }
    file.close();
    return true;
}

std::vector<std::vector<double>> convertToVectorMatrix(double **returnMatrix, int numberAssets, int numberReturns) {
    std::vector<std::vector<double>> vectorMatrix(numberAssets, std::vector<double>(numberReturns));

    for (int i = 0; i < numberAssets; ++i) {
        for (int j = 0; j < numberReturns; ++j) {
            vectorMatrix[i][j] = returnMatrix[i][j];
        }
    }

    return vectorMatrix;
}

void deleteDoublePointer(double** matrix, int numberRows) {
    for (int i = 0; i < numberRows; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
    std::cout << "Deleted Memory" << std::endl;
}