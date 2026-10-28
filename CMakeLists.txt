cmake_minimum_required(VERSION 3.21)

# Project name and version
project(TerraNeo VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable verbose makefile output
set(CMAKE_VERBOSE_MAKEFILE OFF)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "16.0" )
  message( FATAL_ERROR "TERRA-NG requires at least version 16 of Clang. Found: ${CMAKE_CXX_COMPILER_VERSION}" )
endif ()

# Make sure the serial backend is always build
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Always enable Serial backend" FORCE)

# Optionally build against ADIOS2 for I/O
option( TNG_BUILD_WITH_ADIOS2 "Build with ADIOS2" OFF )

if( TNG_BUILD_WITH_ADIOS2 )

  # We need an ADIOS2 installation built against a matching Kokkos installation, so we
  # enforce using an external Kokkos here
  set( TNG_USE_EXTERNAL_KOKKOS ON CACHE BOOL "Use an external installation of Kokkos" FORCE)
  find_package( Kokkos REQUIRED )

  # Find ADIOS2
  find_package( ADIOS2 REQUIRED COMPONENTS CXX MPI )

endif()

# Always make compiler output coloured by default
set( TNG_COLOUR_COMPILER_OUTPUT ON CACHE BOOL "Make compiler messages use colour" )
mark_as_advanced( TNG_COLOUR_COMPILER_OUTPUT )
if( TNG_COLOUR_COMPILER_OUTPUT )
  add_compile_options($<$<CXX_COMPILER_ID:Clang>:-fdiagnostics-color=always>)
  add_compile_options($<$<CXX_COMPILER_ID:GNU>:-fdiagnostics-color=always>)
  if( CMAKE_CXX_COMPILER_ID MATCHES "NVHPC" )
    add_compile_options( -Xcompiler -fdiagnostics-color=always )
  endif()
endif()

# Suppress selected compiler warnings
set( TNG_SUPPRESS_WARNINGS "-Wno-unknown-pragmas -Wno-sign-compare" CACHE STRING "Suppress certain compiler warnings" )
mark_as_advanced( TNG_SUPPRESS_WARNINGS )
if( NOT "${TNG_SUPPRESS_WARNINGS}" STREQUAL "" )
  separate_arguments( TNG_SUPPRESS_LIST NATIVE_COMMAND "${TNG_SUPPRESS_WARNINGS}" )
  # why does the generator expression not work?
  # add_compile_options( $<$<CXX_COMPILER_ID:Clang>:$<JOIN:${TNG_SUPPRESS_LIST},$<SEMICOLON>>> )
  if( CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
    add_compile_options( ${TNG_SUPPRESS_LIST} )
  endif()
  message( STATUS "Compiler uses the following suppressions: ${TNG_SUPPRESS_WARNINGS}" )
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND Kokkos_ENABLE_HIP )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -munsafe-fp-atomics -ffast-math -Rpass-analysis=kernel-resource-usage")
endif ()

find_package(MPI REQUIRED)
if (MPI_FOUND)
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH})
    foreach (LIB ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
        if (LIB)
            list(APPEND SERVICE_LIBS ${LIB})
        endif ()
    endforeach ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")

    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
endif ()

set(EIGEN_DIR extern/eigen-3.4.0 CACHE PATH "eigen path")
include_directories(SYSTEM ${EIGEN_DIR})
link_directories(${EIGEN_DIR})

# Check if both CUDA/NVHPC and OpenMP are being used and enforce correct openmp flag for GCC
if( CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC" AND Kokkos_ENABLE_OPENMP )
  if( TARGET OpenMP::OpenMP_CXX )
    set_target_properties( OpenMP::OpenMP_CXX PROPERTIES
      INTERFACE_COMPILE_OPTIONS "-Xcompiler;-fopenmp"
    )
    message( STATUS "Overriding properties of target OpenMP::OpenMP_CXX with '-Xcompiler;-fopenmp'" )
  endif()
  # It seems that OpenMP_CXX_FLAGS doesn't get used anymore (but better safe than sorry)
  set( OpenMP_CXX_FLAGS "-Xcompiler -fopenmp" CACHE STRING "Flags used by the CXX compiler during OpenMP builds" FORCE )
endif()

# Enable testing
enable_testing()

# Add the src and tests subdirectories
add_subdirectory(apps)
add_subdirectory(extern)
add_subdirectory(src)
add_subdirectory(tests)

# A custom target to be built by the CI.
# Mostly to limit compile time and still cover a bit of everything.
add_custom_target(ci-subset
        DEPENDS
        mantlecirculation
        benchmark_baseline
        benchmark_operators
        read_checkpoint
        visualize_radial_profile
        visualize_spherical_harmonics
        visualize_thick_spherical_shell_mesh
        test_laplace_multigrid
        test_epsilon_divdiv_stokes
        test_unsteady_advection_diffusion_supg)
