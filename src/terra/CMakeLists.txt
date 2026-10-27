add_library(terra)

target_link_libraries(terra PUBLIC Kokkos::kokkos MPI::MPI_CXX)

if( TNG_BUILD_WITH_ADIOS2 )
  target_link_libraries(terra PUBLIC adios2::cxx_mpi)
endif()

target_include_directories(terra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(terra
        PRIVATE
        library_dummy.cpp
	terra_definitions.in.hpp
        types.hpp)

if(NESMIK_HOME)
	target_include_directories(terra PUBLIC ${NESMIK_HOME}/include)
	target_link_directories(terra PUBLIC ${NESMIK_HOME}/lib ${NESMIK_HOME}/lib64)
	target_link_libraries(terra PUBLIC nesmik)
        target_compile_definitions(terra PUBLIC TERRANEO_USE_NESMIK)
endif()

# Generate preprocessor macros from TNG options
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/terra_definitions.in.hpp
                ${CMAKE_CURRENT_BINARY_DIR}/terra_definitions.hpp )
target_include_directories( terra PUBLIC ${CMAKE_CURRENT_BINARY_DIR} )

add_subdirectory(communication)
add_subdirectory(dense)
add_subdirectory(eigen)
add_subdirectory(fe)
add_subdirectory(geophysics)
add_subdirectory(grid)
add_subdirectory(io)
add_subdirectory(kokkos)
add_subdirectory(linalg)
add_subdirectory(mpi)
add_subdirectory(shell)
add_subdirectory(util)
add_subdirectory(plates)
