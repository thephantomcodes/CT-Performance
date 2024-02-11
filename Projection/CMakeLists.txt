cmake_minimum_required(VERSION 3.5)

project(Projection LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SOURCES
    ./src/*.cpp
)

file(GLOB HEADERS
    ./include/*.h
)

include_directories("${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}/include")

add_executable(paramtest ${SOURCES} ${HEADERS})
