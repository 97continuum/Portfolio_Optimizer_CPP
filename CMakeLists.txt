cmake_minimum_required(VERSION 3.27)
project(Portfolio_Optimizer_CPP)

set(CMAKE_CXX_STANDARD 14)

include_directories(src/cpp)

add_executable(Portfolio_Optimizer_CPP
        src/cpp/csv.cpp
        src/cpp/linearAlgebra.cpp
        src/cpp/linearAlgebra.h
        src/cpp/main.cpp
        src/cpp/portfolio.cpp
        src/cpp/portfolio.h
        src/cpp/read_data.cpp
        src/cpp/read_data.h
)
