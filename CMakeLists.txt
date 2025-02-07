cmake_minimum_required(VERSION 3.20)

# Project name
project(TurboGen)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable building shared libraries (DLLs)
if(WIN32)
    set(BUILD_SHARED_LIBS ON)
endif()

# Set fmt and CoolProp directories
set(fmt_SOURCE_DIR "${CMAKE_SOURCE_DIR}/externals/fmt")
set(COOLPROP_DIR "${CMAKE_SOURCE_DIR}/externals/CoolProp")

# Platform-specific configuration
if(APPLE)
    find_package(fmt REQUIRED)
elseif(WIN32)
    # Windows: Use the existing local fmt subdirectory
    add_subdirectory(${fmt_SOURCE_DIR} ${CMAKE_BINARY_DIR}/fmt)
endif()

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
    src/Aungier.cpp
    src/plotter.cpp
    src/cli.cpp
    correlations/prelimCalcs.cpp
)

# Include CoolProp and fmt headers
target_include_directories(TurboGen PRIVATE 
    ${COOLPROP_INCLUDE_DIR} 
    ${fmt_INCLUDE_DIR}
)

# Platform-specific linking
if(APPLE)
    # macOS: Link fmt using find_package
    target_link_libraries(TurboGen PRIVATE 
        ${COOLPROP_DIR}/lib/libCoolProp.a
        fmt::fmt)
elseif(WIN32)
    # Windows: Link fmt as a subdirectory (static or shared, depending on your setup)
    target_link_libraries(TurboGen PRIVATE 
        ${COOLPROP_DIR}/lib/CoolPropd.lib
        fmt::fmt                            
    )
endif()

# Add definitions for using C++11
add_definitions(-DCOOLPROP_USE_CXX11)

# Ensure fmt.dll is available at runtime
if(WIN32)
    add_custom_command(TARGET TurboGen POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_BINARY_DIR}/fmt/bin/Debug/fmtd.dll ${CMAKE_BINARY_DIR}/Debug
    )
endif()
