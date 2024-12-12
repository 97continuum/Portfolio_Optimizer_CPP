#ifndef COURSEWORK_READ_DATA_H
#define COURSEWORK_READ_DATA_H

#include <string>
#include "linearAlgebra.h"

double string_to_double( const string& s );
void readData(Matrix& data, const string& fileName);
bool changeWorkingDirectory(const string& newDir);
bool checkFileInCurrentDirectory(const string& fileName);
void deleteDoublePointer(double** matrix, int numberRows);
Matrix convertToVectorMatrix(double **returnMatrix, int numberAssets, int numberReturns);
#endif //COURSEWORK_READ_DATA_H
