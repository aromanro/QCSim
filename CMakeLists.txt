cmake_minimum_required(VERSION 3.10)
project(QCSim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()

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
            QCSim/VQETests.cpp
			QCSim/MPSSimulatorTests.cpp
			QCSim/CliffordTests.cpp)
 
SET( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please set EIGEN3_INCLUDE_DIR env variable to the include directory of Eigen3")
ENDIF()

SET( FFTW3_DIR "$ENV{FFTW3_DIR}" )
IF( NOT FFTW3_DIR )
    MESSAGE( FATAL_ERROR "Please set FFTW3_DIR env variable to the include directory of fftw3")
ENDIF()

if(APPLE)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
      set(OpenMP_C "${CMAKE_C_COMPILER}")
      set(OpenMP_C_FLAGS "-fopenmp")
      set(OpenMP_CXX_LIB_NAMES "omp")
      set(OpenMP_omp_LIBRARY omp)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      set(OpenMP_CXX_FLAGS "-fopenmp")
      set(OpenMP_CXX_LIB_NAMES "omp")
      set(OpenMP_omp_LIBRARY omp)
    endif()
endif()

find_package(OpenMP REQUIRED)

include(CheckCXXCompilerFlag)
function(enable_cxx_flag flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_set)
    if(flag_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()

if(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64" OR CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "amd64" OR CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "AMD64")
    if(UNIX OR APPLE)
        if (NOT CMAKE_OSX_ARCHITECTURES STREQUAL "arm64")
            set(SIMD_FLAGS "-mfma;-mavx2")
            enable_cxx_flag("-mpopcnt")
        endif()
    elseif(MSVC)
	    set(SIMD_FLAGS "/arch:AVX2")
    endif()
endif()

INCLUDE_DIRECTORIES ( QCSim "${EIGEN3_INCLUDE_DIR}" "${FFTW3_DIR}" )
LINK_DIRECTORIES (${FFTW3_DIR})
IF(WIN32)
    LINK_LIBRARIES(libfftw3-3)
ELSE()
    LINK_LIBRARIES(${FFTW3_DIR}/../lib/libfftw3.a ${FFTW3_DIR}/../lib/libfftw3_threads.a)
ENDIF()

add_executable (QCSim ${SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(QCSim PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(QCSim PRIVATE ${OpenMP_CXX_FLAGS})
    target_compile_options(QCSim PRIVATE ${OpenMP_CXX_FLAGS})
endif()

target_compile_options(QCSim PUBLIC ${SIMD_FLAGS})