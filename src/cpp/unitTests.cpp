//
// Created by Talha Jamal on 20/05/2024.
//

#include "unitTests.h"
#include "iostream"
#include <vector>
#include "linearAlgebra.h"



void testVectorAdd() {
    std::vector<double> a = {1.0, 2.0, 3.0};
    std::vector<double> b = {4.0, 5.0, 6.0};
    std::vector<double> result = vectorAddition(a, b);

    std::cout << "Vector Addition Result: ";
    for (double val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void testVectorSub() {
    std::vector<double> a = {4.0, 5.0, 6.0};
    std::vector<double> b = {1.0, 2.0, 3.0};
    std::vector<double> result = vectorSubtraction(a, b);

    std::cout << "Vector Subtraction Result: ";
    for (double val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void testDotProduct() {
    std::vector<double> a = {1.0, 2.0, 3.0};
    std::vector<double> b = {4.0, 5.0, 6.0};
    double result = vectorDotProduct(a, b);

    std::cout << "Dot Product Result: " << result << std::endl;
}

void testScalarMult() {
    std::vector<double> a = {1.0, 2.0, 3.0};
    double scalar = 2.0;
    std::vector<double> result = scalarMultiplication(a, scalar);

    std::cout << "Scalar Multiplication Result: ";
    for (double val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void testMatVecMult() {
    std::vector<std::vector<double>> mat = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0}
    };
    std::vector<double> vec = {1.0, 2.0, 3.0};
    std::vector<double> result = matrixVectorMultiplication(mat, vec);

    std::cout << "Matrix-Vector Multiplication Result: ";
    for (double val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void testMatMatMult() {
    std::vector<std::vector<double>> A = {
            {1.0, 2.0},
            {3.0, 4.0},
            {5.0, 6.0}
    };
    std::vector<std::vector<double>> B = {
            {7.0, 8.0, 9.0},
            {10.0, 11.0, 12.0}
    };
    std::vector<std::vector<double>> result = matrixMultiplication(A, B);

    std::cout << "Matrix-Matrix Multiplication Result: " << std::endl;
    for (const auto& row : result) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

void testAllFunctions(){
    std::cout << "Testing Vector Addition..." << std::endl;
    testVectorAdd();
    std::cout << std::endl;

    std::cout << "Testing Vector Subtraction..." << std::endl;
    testVectorSub();
    std::cout << std::endl;

    std::cout << "Testing Dot Product..." << std::endl;
    testDotProduct();
    std::cout << std::endl;

    std::cout << "Testing Scalar Multiplication..." << std::endl;
    testScalarMult();
    std::cout << std::endl;

    std::cout << "Testing Matrix-Vector Multiplication..." << std::endl;
    testMatVecMult();
    std::cout << std::endl;

    std::cout << "Testing Matrix-Matrix Multiplication..." << std::endl;
    testMatMatMult();
    std::cout << std::endl;
}