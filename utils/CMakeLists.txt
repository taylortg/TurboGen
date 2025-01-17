cmake_minimum_required(VERSION 3.20)

# Project name
project(Utils)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable building shared libraries (DLLs)
set(BUILD_SHARED_LIBS ON)

# Set fmt and CoolProp directories
set(fmt_SOURCE_DIR "${CMAKE_SOURCE_DIR}/../externals/fmt")

# Add fmt as a subdirectory if it’s part of your project
add_subdirectory(${fmt_SOURCE_DIR} ${CMAKE_BINARY_DIR}/fmt)

# Include directories for fmt and CoolProp
set(fmt_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../externals/fmt/include")

# Add the executable
add_executable(Utils 
    cmds/makeBlade.cpp
    common/config.cpp
)

# Include CoolProp and fmt headers
target_include_directories(Utils PRIVATE ${fmt_INCLUDE_DIR})

# Link to the fmt shared library and CoolProp static library
target_link_libraries(Utils PRIVATE 
    fmt::fmt                            # Link to the fmt shared library (DLL)
)

# Ensure fmt.dll is available at runtime
add_custom_command(TARGET Utils POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_BINARY_DIR}/fmt/bin/Debug/fmtd.dll ${CMAKE_BINARY_DIR}/Debug
)
