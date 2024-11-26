cmake_minimum_required(VERSION 3.10)

# export CMAKE_PREFIX_PATH="/home/vol07/scarf237/src/specle:${CMAKE_PREFIX_PATH}"

# set the project name
project(specle LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(FFTW3_DIR /work4/scd/scarf1354/work/proj2/specle/)
# set(FFTW3_DIR /work/c01/c01/terence411/work/specle/)

find_package(MPI REQUIRED)

# find_package(Kokkos REQUIRED PATHS ${KOKKOS_DIR})
find_package(HDF5 REQUIRED)
find_package(FFTW3 REQUIRED PATHS ${FFTW3_DIR})

# find_package(FFTW3 REQUIRED)

# set(THREADS_PREFER_PTHREAD_FLAG ON)
# find_package(Threads REQUIRED)

add_executable(specle
    src/commands/basicio.cpp
    src/commands/hdfio.cpp
    src/commands/mutate.cpp
    src/commands/setters.cpp
    src/commands/setup.cpp
    src/data/data-util.cpp
    src/specle.cpp
    src/globals.cpp
    src/io.cpp
    src/parse.cpp
    src/transfer.cpp
)

target_include_directories(specle PUBLIC ${FFTW3_INCLUDE_DIRS} ${MPI_INCLUDE_PATH} ${HDF5_INCLUDE_DIRS})
target_link_libraries(specle PUBLIC fftw3_mpi fftw3 ${MPI_LIBRARIES} ${HDF5_LIBRARIES})

# target_include_directories(specle PUBLIC ${FFTW3_INCLUDE_DIRS})
# target_link_libraries(specle PUBLIC fftw3 fftw3_threads)
# target_link_libraries(specle PRIVATE Threads::Threads)

# target_compile_definitions(specle PRIVATE KOKKOS_DEPENDENCE)

target_link_directories(specle PUBLIC ${FFTW3_LIBRARY_DIRS})

add_executable(dataToImg EXCLUDE_FROM_ALL src/dataToImg.cpp)
add_executable(dataToAscii EXCLUDE_FROM_ALL src/dataToAscii.cpp)

message("MPI ${MPI_INCLUDE_PATH} | ${MPI_LIBRARIES}")
