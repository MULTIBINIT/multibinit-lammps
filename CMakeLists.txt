# -*- CMake -*- file for Multibinit LAMMPS plugin
# Build instructions:
#   mkdir build && cd build
#   cmake .. -DLAMMPS_SOURCE_DIR=/path/to/lammps/src -DABINIT_ROOT=/path/to/abinit
#   make
#
# Usage in LAMMPS:
#   plugin load /path/to/multibinitplugin.so
#   pair_style multibinit 2 2 2 ngqpt 2 2 2 dipdip 1
#   pair_coeff * * BaTiO3_DDB
##########################################

cmake_minimum_required(VERSION 3.20)

# Enforce out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. Create a build directory.")
endif()

project(multibinitplugin VERSION 1.0 LANGUAGES CXX)

# LAMMPS source directory (required)
if(NOT LAMMPS_SOURCE_DIR)
  get_filename_component(LAMMPS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/../lammps/src ABSOLUTE)
endif()
set(LAMMPS_HEADER_DIR ${LAMMPS_SOURCE_DIR} CACHE PATH "Location of LAMMPS headers")
if(NOT EXISTS ${LAMMPS_HEADER_DIR}/lammps.h)
  message(FATAL_ERROR "Must set LAMMPS_SOURCE_DIR to LAMMPS src directory")
endif()

# Default install to $HOME/.local
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Default install path" FORCE)
endif()

# C++17 required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

# Debug flags for bounds checking and uninitialized variables
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wuninitialized -fsanitize=address")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wuninitialized -fcheck=bounds -finit-integer=-999 -finit-real=snan -ffpe-trap=invalid,zero,overflow")

# Enable debug mode
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type")

# MPI is required
find_package(MPI REQUIRED COMPONENTS CXX)
set(MPI_CXX_SKIP_MPICXX TRUE)

# Find ABINIT
set(ABINIT_ROOT "" CACHE PATH "Path to ABINIT installation")
find_path(ABINIT_INCLUDE_DIR
  NAMES m_mb_cbinding.mod m_mb_cbinding.h
  PATHS
    ${ABINIT_ROOT}/include
    ${ABINIT_ROOT}/modules
    ${ABINIT_ROOT}/../modules
    $ENV{ABINIT_ROOT}/include
    $ENV{ABINIT_ROOT}/modules
    $ENV{ABINIT_DIR}/include
    $ENV{ABINIT_DIR}/modules
)

# ABINIT libraries - may be a single libabinit or multiple module libraries
# Search in both lib/ and install/lib/ (ABINIT CMake builds put libs in install/lib)
find_library(ABINIT_LIBRARY
  NAMES abinit libabinit
  PATHS
    ${ABINIT_ROOT}/install/lib
    ${ABINIT_ROOT}/lib
    $ENV{ABINIT_ROOT}/install/lib
    $ENV{ABINIT_ROOT}/lib
    $ENV{ABINIT_DIR}/install/lib
    $ENV{ABINIT_DIR}/lib
)

# Look for effpot module (contains Multibinit)
find_library(EFFPOT_LIBRARY
  NAMES effpot lib78_effpot
  PATHS
    ${ABINIT_ROOT}/install/lib
    ${ABINIT_ROOT}/lib
    $ENV{ABINIT_ROOT}/install/lib
    $ENV{ABINIT_ROOT}/lib
)

# Look for common ABINIT modules
find_library(ABINIT_COMMON_LIBRARY
  NAMES abi_common lib01_common
  PATHS
    ${ABINIT_ROOT}/install/lib
    ${ABINIT_ROOT}/lib
    $ENV{ABINIT_ROOT}/install/lib
    $ENV{ABINIT_ROOT}/lib
)

find_library(XMPI_LIBRARY
  NAMES xmpi lib12_hide_mpi
  PATHS
    ${ABINIT_ROOT}/install/lib
    ${ABINIT_ROOT}/lib
    $ENV{ABINIT_ROOT}/install/lib
    $ENV{ABINIT_ROOT}/lib
)

# Fortran runtime (required for linking Fortran code)
find_library(GFORTRAN_LIBRARY gfortran)
find_library(M_LIBRARY m)

# NetCDF (required by ABINIT)
find_library(NETCDF_LIBRARY netcdf)
find_library(NETCDFF_LIBRARY netcdff)

# LibXC (required by ABINIT)
find_library(XC_LIBRARY xc)

if(ABINIT_INCLUDE_DIR)
  message(STATUS "Found ABINIT includes: ${ABINIT_INCLUDE_DIR}")
else()
  message(WARNING "ABINIT include directory not found. Set ABINIT_ROOT.")
endif()

# Collect ABINIT libraries
set(ABINIT_LIBS "")
if(EFFPOT_LIBRARY)
  list(APPEND ABINIT_LIBS ${EFFPOT_LIBRARY})
  message(STATUS "Found effpot: ${EFFPOT_LIBRARY}")
endif()
if(XMPI_LIBRARY)
  list(APPEND ABINIT_LIBS ${XMPI_LIBRARY})
  message(STATUS "Found xmpi: ${XMPI_LIBRARY}")
endif()
if(ABINIT_COMMON_LIBRARY)
  list(APPEND ABINIT_LIBS ${ABINIT_COMMON_LIBRARY})
  message(STATUS "Found abi_common: ${ABINIT_COMMON_LIBRARY}")
endif()
if(ABINIT_LIBRARY)
  list(APPEND ABINIT_LIBS ${ABINIT_LIBRARY})
  message(STATUS "Found abinit: ${ABINIT_LIBRARY}")
endif()
if(NETCDFF_LIBRARY)
  list(APPEND ABINIT_LIBS ${NETCDFF_LIBRARY})
  message(STATUS "Found netcdff: ${NETCDFF_LIBRARY}")
endif()
if(NETCDF_LIBRARY)
  list(APPEND ABINIT_LIBS ${NETCDF_LIBRARY})
  message(STATUS "Found netcdf: ${NETCDF_LIBRARY}")
endif()
if(XC_LIBRARY)
  list(APPEND ABINIT_LIBS ${XC_LIBRARY})
  message(STATUS "Found xc: ${XC_LIBRARY}")
endif()
if(GFORTRAN_LIBRARY)
  list(APPEND ABINIT_LIBS ${GFORTRAN_LIBRARY})
  message(STATUS "Found gfortran: ${GFORTRAN_LIBRARY}")
endif()
if(M_LIBRARY)
  list(APPEND ABINIT_LIBS ${M_LIBRARY})
endif()

# Create LAMMPS interface library (for plugin builds)
add_library(lammps INTERFACE)
target_include_directories(lammps INTERFACE
  ${LAMMPS_HEADER_DIR}
  ${LAMMPS_HEADER_DIR}/OPENMP  # For OMP compatibility if needed
)
target_link_libraries(lammps INTERFACE MPI::MPI_CXX)

# Build the plugin as a shared library
add_library(multibinitplugin MODULE
  src/multibinitplugin.cpp
  src/pair_multibinit.cpp
)

target_include_directories(multibinitplugin PRIVATE
  ${LAMMPS_HEADER_DIR}
  ${ABINIT_INCLUDE_DIR}
  ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(multibinitplugin PRIVATE
  lammps
  ${ABINIT_LIBS}
)

# Set plugin properties
set_target_properties(multibinitplugin PROPERTIES
  PREFIX ""
  OUTPUT_NAME "multibinitplugin"
)

# Platform-specific settings
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set_target_properties(multibinitplugin PROPERTIES
    SUFFIX ".so"
    # dynamic_lookup allows undefined symbols (for LAMMPS symbols loaded at runtime)
    # But we still link ABINIT to get those symbols
    LINK_FLAGS "-Wl,-undefined,dynamic_lookup"
  )
  # Set rpath to find ABINIT library at runtime
  set_target_properties(multibinitplugin PROPERTIES
    INSTALL_RPATH "${ABINIT_ROOT}/install/lib"
    BUILD_WITH_INSTALL_RPATH TRUE
  )
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set_target_properties(multibinitplugin PROPERTIES
    WINDOWS_EXPORT_ALL_SYMBOLS TRUE
  )
else()
  set_target_properties(multibinitplugin PROPERTIES
    SUFFIX ".so"
    LINK_FLAGS "-rdynamic"
  )
endif()

# Install
install(TARGETS multibinitplugin
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

# Print summary
message(STATUS "")
message(STATUS "Multibinit LAMMPS Plugin Configuration Summary:")
message(STATUS "  LAMMPS source: ${LAMMPS_HEADER_DIR}")
message(STATUS "  ABINIT root: ${ABINIT_ROOT}")
message(STATUS "  ABINIT includes: ${ABINIT_INCLUDE_DIR}")
message(STATUS "  ABINIT libs: ${ABINIT_LIBS}")
message(STATUS "  Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
