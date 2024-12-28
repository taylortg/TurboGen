cmake_minimum_required(VERSION 3.10)

# Project name
project(TurboGen)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(COOLPROP_LIB_DIR "${CMAKE_SOURCE_DIR}/externals/CoolProp/lib/CoolPropd.lib")

# Set the CoolProp include directory
set(COOLPROP_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/externals/CoolProp/include")

# Add the executable
add_executable(TurboGen 
    src/main.cpp
    src/common.cpp
    src/tgParser.cpp
    src/thermo.cpp
    src/impeller.cpp
)

# Include the CoolProp headers
target_include_directories(TurboGen PRIVATE ${COOLPROP_INCLUDE_DIR})

# Link the CoolProp library
target_link_libraries(TurboGen PRIVATE ${COOLPROP_LIB_DIR})

# Add definitions for using C++11
add_definitions(-DCOOLPROP_USE_CXX11)
