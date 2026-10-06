cmake_minimum_required(VERSION 3.19)

if(APPLE)
    message(STATUS "Detected Apple platform (macOS)")

    message(STATUS "Attempting to detect macOS SDK path using xcrun...")
    execute_process(
        COMMAND xcrun --sdk macosx --show-sdk-path
        OUTPUT_VARIABLE MACOSX_SDK_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    message(STATUS "macOS SDK path detected: ${MACOSX_SDK_PATH}")

    set(CMAKE_OSX_SYSROOT "${MACOSX_SDK_PATH}" CACHE STRING "macOS SDK path")
    message(STATUS "Set CMAKE_OSX_SYSROOT to: ${CMAKE_OSX_SYSROOT}")

    # Optionally also set deployment target
    set(CMAKE_OSX_DEPLOYMENT_TARGET "12.0" CACHE STRING "Minimum macOS version")
    message(STATUS "Set CMAKE_OSX_DEPLOYMENT_TARGET to: ${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif()

project(fronts)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Armadillo and HDF5 packages
find_package(Armadillo REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)

# Find SuperLU
find_library(SUPERLU_LIB superlu HINTS /opt/homebrew/lib /usr/local/lib)
if(NOT SUPERLU_LIB)
    message(FATAL_ERROR "SuperLU library not found! Pass -DSUPERLU_LIB=/path/to/libsuperlu")
endif()
find_path(SUPERLU_INCLUDE_DIR "slu_ddefs.h" HINTS /opt/homebrew/include /usr/local/include PATH_SUFFIXES superlu)

# Create the executable
add_executable(fronts main.cxx InitialCondition.h Coefficient.h FisherBasicSolver.h FisherBasicSolver.cxx Fisher23Solver.h Fisher23Solver.cxx options.h options.cxx H5DatasetBatchLoader.h GreensIntegrals.h)
target_compile_options(fronts PRIVATE -Wall -Wextra)

# Include directories and link libraries
target_link_libraries(fronts PRIVATE ${ARMADILLO_LIBRARIES} ${HDF5_LIBRARIES} ${SUPERLU_LIB})
target_include_directories(fronts PRIVATE ${ARMADILLO_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${SUPERLU_INCLUDE_DIR})

# Required definitions for HDF5
add_definitions(-DARMA_USE_SUPERLU)
target_compile_definitions(fronts PRIVATE ${HDF5_DEFINITIONS})

# Output CMake configuration debug information
message(STATUS "Armadillo include dirs: ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "SuperLU include dir: ${SUPERLU_INCLUDE_DIR}")
