cmake_minimum_required(VERSION 3.1.0)
project (SpheroidalScattering)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O0")
find_package(MPI REQUIRED)

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)

include_directories(include test ${MPI_INCLUDE_PATH} "/usr/include/python3.7m")

set(PROJECT_LINK_LIBS ${PROJECT_LINK_LIBS} "libspecialfun.a" "libsplinter-static-3-0.a")
link_directories("lib" "lib/suitesparse")

file(GLOB SOURCES "src/*.cpp" "test/*.cpp")

#${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} hello

add_executable(main ${SOURCES})

## threading: -lpthread
##for Boost.Python link to: -lpython3.7m -lboost_python37 -lboost_numpy37 -lboost_system -lboost_filesystem
target_link_libraries(main ${PROJECT_LINK_LIBS} ${MPI_CXX_LIBRARIES} -lgfortran -llapack -llapacke -lopenblas -lumfpack -lpthread)

