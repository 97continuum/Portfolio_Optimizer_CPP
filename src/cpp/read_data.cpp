#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "csv.h"  // Assuming you have a Csv class defined somewhere

using namespace std;

double string_to_double(const string& s) {
    try {
        return stod(s);
    } catch (const invalid_argument& e) {
        cerr << "Invalid argument: " << e.what() << endl;
        return 0.0;
    } catch (const out_of_range& e) {
        cerr << "Out of range: " << e.what() << endl;
        return 0.0;
    }
}

void readData(double **data, const string& fileName) {
    //cout << "Attempting to open file: " << fileName << endl;

    ifstream file(fileName);
    if (!file.is_open()) {
        cerr << "Failed to open file: " << fileName << endl;
        return;
    }

    Csv csv(file);
    string line;
    int i = 0;
    while (csv.getline(line) != 0) {
        for (int j = 0; j < csv.getnfield(); j++) {
            double temp = string_to_double(csv.getfield(j));
            data[j][i] = temp;
        }
        i++;
    }
    file.close();
    //cout << "File read successfully: " << fileName << endl;
}

bool changeWorkingDirectory(const string& newDir) {
    if (chdir(newDir.c_str()) != 0)
    {
        perror("chdir() error");
        cerr << "Failed to change working directory to: " << newDir << endl;
        return false;
    }
    else
    {
        //cout << "Successfully changed working directory to: " << newDir << endl;
        return true;
    }
}

bool checkFileInCurrentDirectory(const string& fileName) {
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        cout << "Current working dir: " << cwd << endl;
    } else {
        perror("getcwd() error");
        return false;
    }

    string fullPath = string(cwd) + "/" + fileName;
    ifstream file(fullPath);
    if (!file) {
        cerr << "File " << fullPath << " does not exist." << endl;
        return false;
    } else {
        //cout << "File " << fullPath << " exists" << endl;
    }
    file.close();
    return true;
}

Matrix convertToVectorMatrix(double **returnMatrix, int numberAssets, int numberReturns) {
    Matrix vectorMatrix(numberAssets, Vector(numberReturns));

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
    cout << "Deleted Memory" << endl;
}