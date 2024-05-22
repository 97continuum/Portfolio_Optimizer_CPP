//
// Created by Talha Jamal on 19/05/2024.
//

#ifndef COURSEWORK_READ_DATA_H
#define COURSEWORK_READ_DATA_H

#include <string>

double string_to_double( const std::string& s );
void readData(double **data,const std::string& fileName);
bool checkFileInCurrentDirectory(const std::string& fileName);
std::vector< std::vector<double> > convertToVectorMatrix(double **returnMatrix, int numberAssets, int numberReturns);
#endif //COURSEWORK_READ_DATA_H
