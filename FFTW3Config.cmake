# defined since 2.8.3
if (CMAKE_VERSION VERSION_LESS 2.8.3)
  get_filename_component (CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
endif ()


# Allows loading FFTW3 settings from another project
set (FFTW3_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")

set (FFTW3_LIBRARIES fftw3)
set (FFTW3_LIBRARY_DIRS /work4/scd/scarf1354/work/proj2/resources/lib)
set (FFTW3_INCLUDE_DIRS /work4/scd/scarf1354/work/proj2/resources/include)

# set (FFTW3_LIBRARY_DIRS /work4/scd/scarf1354/work/proj2/resources/lib)
# set (FFTW3_INCLUDE_DIRS /work4/scd/scarf1354/work/proj2/resources/include)

# set (FFTW3_LIBRARY_DIRS /home/vol07/scarf237/arch/amd/lib)
# set (FFTW3_INCLUDE_DIRS /home/vol07/scarf237/arch/amd/include)

message("-- FFTW3Config.cmake File --")
message("CMAKE_CURRENT_LIST_DIR  : ${CMAKE_CURRENT_LIST_DIR}")
message("FFTW3_CONFIG_FILE       : ${CMAKE_CURRENT_LIST_FILE}")
message("FFTW3_LIBRARY_DIRS      : ${FFTW3_LIBRARY_DIRS}")
message("FFTW3_INCLUDE_DIRS      : ${FFTW3_INCLUDE_DIRS}")

include ("${CMAKE_CURRENT_LIST_DIR}/FFTW3LibraryDepends.cmake")

if (CMAKE_VERSION VERSION_LESS 2.8.3)
  set (CMAKE_CURRENT_LIST_DIR)
endif ()