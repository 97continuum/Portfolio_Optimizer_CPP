//
// Created by Talha Jamal on 19/05/2024.
//

#ifndef COURSEWORK_READ_DATA_H
#define COURSEWORK_READ_DATA_H

#include <string>

double string_to_double( const string& s );
void readData(double **data,const string& fileName);
bool changeWorkingDirectory(const string& newDir);
bool checkFileInCurrentDirectory(const string& fileName);
void deleteDoublePointer(double** matrix, int numberRows);
MatrixconvertToVectorMatrix(double **returnMatrix, int numberAssets, int numberReturns);
#endif //COURSEWORK_READ_DATA_H
