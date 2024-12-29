cmake_minimum_required(VERSION 3.20)

# Project name
project(TurboGen)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable building shared libraries (DLLs)
set(BUILD_SHARED_LIBS ON)

# Set fmt and CoolProp directories
set(fmt_SOURCE_DIR "${CMAKE_SOURCE_DIR}/externals/fmt")
set(COOLPROP_DIR "${CMAKE_SOURCE_DIR}/externals/CoolProp")

# Add fmt as a subdirectory if it’s part of your project
add_subdirectory(${fmt_SOURCE_DIR} ${CMAKE_BINARY_DIR}/fmt)

# Include directories for fmt and CoolProp
set(fmt_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/externals/fmt/include")
set(COOLPROP_INCLUDE_DIR "${COOLPROP_DIR}/include")

# Add the executable
add_executable(TurboGen 
    src/main.cpp
    src/common.cpp
    src/tgParser.cpp
    src/thermo.cpp
    src/impeller.cpp
)

# Include CoolProp and fmt headers
target_include_directories(TurboGen PRIVATE ${COOLPROP_INCLUDE_DIR} ${fmt_INCLUDE_DIR})

# Link to the fmt shared library and CoolProp static library
target_link_libraries(TurboGen PRIVATE 
    ${COOLPROP_DIR}/lib/CoolPropd.lib   # Static link to CoolPropd
    fmt::fmt                            # Link to the fmt shared library (DLL)
)

# Add definitions for using C++11
add_definitions(-DCOOLPROP_USE_CXX11)

# Ensure fmt.dll is available at runtime
add_custom_command(TARGET TurboGen POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_BINARY_DIR}/fmt/bin/Debug/fmtd.dll ${CMAKE_BINARY_DIR}/Debug
)
