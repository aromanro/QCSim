cmake_minimum_required(VERSION 3.10)
project(QCSim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(default_build_type "Release")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(SOURCES QCSim/QCSim.cpp 
            QCSim/Tests.cpp 
            QCSim/AdderTests.cpp
            QCSim/ErrorCorrectionTests.cpp
            QCSim/FFT.cpp
            QCSim/GatesTests.cpp
            QCSim/MeasurementsTests.cpp
            QCSim/ParadoxesTests.cpp
            QCSim/SimulationTests.cpp
            QCSim/GamesTests.cpp
            QCSim/DistributedTests.cpp
            QCSim/PhaseEstimationTests.cpp
            QCSim/TeleportationTests.cpp
            QCSim/DeutschJozsaTests.cpp
            QCSim/QMLTests.cpp
            QCSim/IsingTests.cpp
            QCSim/VQETests.cpp)
 
SET( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please set EIGEN3_INCLUDE_DIR env variable to the include directory of Eigen3")
ENDIF()

SET( FFTW3_DIR "$ENV{FFTW3_DIR}" )
IF( NOT FFTW3_DIR )
    MESSAGE( FATAL_ERROR "Please set FFTW3_DIR env variable to the include directory of fftw3")
ENDIF()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

INCLUDE_DIRECTORIES ( QCSim "${EIGEN3_INCLUDE_DIR}" "${FFTW3_DIR}" )
LINK_DIRECTORIES (${FFTW3_DIR})
IF(WIN32)
LINK_LIBRARIES(libfftw3-3)
ELSE()
LINK_LIBRARIES(libfftw3.a libfftw3_threads.a)
ENDIF()

add_executable (QCSim ${SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(QCSim PUBLIC OpenMP::OpenMP_CXX)
endif()
