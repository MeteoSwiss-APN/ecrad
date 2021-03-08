# This CMakeLists.txt file is the entry point for the configuration and the
# compilation of the ECRAD.

FUNCTION(PREPEND var prefix)
   SET(listVar "")
   FOREACH(f ${ARGN})
      LIST(APPEND listVar "${prefix}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION()

cmake_minimum_required(VERSION 3.0)

set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

IF(CMAKE_BUILD_TYPE MATCHES Debug)
  message("Debug build")
  set(BUILD_TYPE "Debug")
ELSEIF(CMAKE_BUILD_TYPE MATCHES Release)
  message("Release build")
  set(BUILD_TYPE "Release")
ELSE()
  message( FATAL_ERROR "Unsupported build type ${CMAKE_BUILD_TYPE}")
ENDIF()

project("ECRAD-MCH ${BUILD_TYPE}" VERSION 1.4.0)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/module")

enable_language(Fortran)
enable_language(CXX)

set(NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)

option(GPU_ACC "Compile for GPU target using OpenACC" OFF)

set(GPU_ACC_TARGET_CC "" CACHE STRING "Compute capability of the GPU target")
set(GPU_ACC_TARGET_MAXREGCOUNT "32" CACHE STRING "Maximum number of registers per thread for GPU target")

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
  include(cmake/compiler/pgi.cmake)
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  include(cmake/compiler/intel.cmake)
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Cray")
  include(cmake/compiler/cray.cmake)
else()
  include(cmake/compiler/gnu.cmake)
endif()

option(SINGLE_PRECISION "Use single precision" OFF)
option(USE_OPENMP "Use OpenMP" ON)
option(USE_DR_HOOK "Use Dr Hook profiling system" OFF)

option(USE_REPORT "Generate XML report" ON)

#if("${NetCDF_VERSION}" VERSION_GREATER_EQUAL "4")
  set(FC_FPPFLAGS "${FC_FPPFLAGS} -DNC_NETCDF")
#endif()

if(SINGLE_PRECISION)
  set(FC_FPPFLAGS "${FC_FPPFLAGS} -DSINGLE_PRECISION")
endif()

if(USE_DR_HOOK)
  message( FATAL_ERROR "Dr hook not supported in this build")
endif()

set(Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

if(USE_REPORT)
set(FC_FPPFLAGS "${FC_FPPFLAGS} -DXML_REPORT_SUPPORTED")
endif()

set(CMAKE_Fortran_FLAGS_RELEASE "${FC_RELEASE_WARNFLAGS} ${FC_BASIC_FLAGS} ${FC_FPPFLAGS} ${FC_RELEASE_OPTFLAGS} ${FC_RELEASE_DEBUGFLAGS} -I ${Fortran_MODULE_DIRECTORY}")
set(CMAKE_Fortran_FLAGS_RELEASE "${FC_RELEASE_WARNFLAGS} ${FC_BASIC_FLAGS} ${FC_FPPFLAGS} ${FC_RELEASE_OPTFLAGS} ${FC_RELEASE_DEBUGFLAGS} -I ${Fortran_MODULE_DIRECTORY}")
set(CMAKE_Fortran_FLAGS_DEBUG "${FC_DEBUG_WARNFLAGS} ${FC_BASIC_FLAGS} ${FC_FPPFLAGS} ${FC_DEBUG_OPTFLAGS} ${FC_DEBUG_DEBUGFLAGS} -I ${Fortran_MODULE_DIRECTORY}")

set(CMAKE_CXX_FLAGS_RELEASE "-fpic")

if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_Fortran_FLAGS}")
endif()

include_directories(${CMAKE_SOURCE_DIR}/include ${NETCDF_INCLUDES})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(ifsaux)
set_target_properties( ifsaux PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(ifsaux)

add_library(drhook)
set_target_properties( drhook PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(drhook)
add_dependencies(drhook ifsaux)

add_library(ifsrrtm)
set_target_properties( ifsrrtm PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(ifsrrtm)
add_dependencies(ifsrrtm drhook)

add_library(utilities)
set_target_properties( utilities PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(utilities)
add_dependencies(utilities drhook ifsaux)

add_library(radiation)
set_target_properties( radiation PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(radiation)
add_dependencies(radiation utilities)

add_library(radsurf STATIC)
set_target_properties( radsurf PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(radsurf)
add_dependencies(radsurf radiation)

add_subdirectory(driver)
message(
  STATUS
  "NETCDF_LIBRARIES_F90: ${NETCDF_LIBRARIES_F90}"
)

if(USE_REPORT)
add_library(driver_report_c)

add_library(driver_report_fortran_wrapper SHARED)
set_target_properties( driver_report_fortran_wrapper PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_dependencies(driver_report_fortran_wrapper driver_report_c)
target_link_libraries(driver_report_fortran_wrapper PRIVATE driver_report_c)
endif()

add_executable(ecrad ${ecrad_SRC})
set_target_properties( ecrad PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
set(STATIC_LIBS ifsaux drhook ifsrrtm utilities radiation radsurf)
if(USE_REPORT)
set(STATIC_LIBS ${STATIC_LIBS})
endif()
set(ECRAD_LIBS ${STATIC_LIBS} ${NETCDF_LIBRARIES_F90})
if(USE_OPENMP)
  set(ECRAD_LIBS ${ECRAD_LIBS} ${OpenMP_Fortran_LIBRARIES})
endif()

target_link_libraries(ecrad PRIVATE ${ECRAD_LIBS} ${ECRAD_LIBS} ${ECRAD_LIBS} driver_report_fortran_wrapper)

find_package(Python 3.6 REQUIRED)

set(PY_VENV_DIR "${CMAKE_BINARY_DIR}/pyenv")

set(PY_PACKAGES xarray numpy netCDF4 generateDS py-cpuinfo matplotlib)

add_custom_target(create-python-virtual-environment)
add_custom_command(
 TARGET create-python-virtual-environment
 COMMAND echo "${Python_EXECUTABLE} -m venv ${PY_VENV_DIR}"
 COMMAND ${Python_EXECUTABLE} -m venv ${PY_VENV_DIR}
 COMMAND ${PY_VENV_DIR}/bin/pip3 install ${PY_PACKAGES}
 COMMENT "Creating python virtual environment"
)

# Regenerate serialization code for driver report
add_custom_target(update-driver-report-format)
add_subdirectory(driver_report)

# Regenerate serialization code for test reports
add_custom_target(update-test-report-formats)

add_custom_target(reset-ci-tests)
add_dependencies(reset-ci-tests ecrad create-python-virtual-environment)
add_subdirectory(ci_tests)
add_custom_command(
   TARGET reset-ci-tests
   COMMAND ${PY_VENV_DIR}/bin/python3 ${CI_TESTS_RUNTIME_DIR}/tests_runner.py --reset-ci-tests
   COMMENT "Generate reference results for ci testcases"
   )

add_custom_target(run-ci-tests)
add_dependencies(run-ci-tests ecrad create-python-virtual-environment)
add_custom_command(
 TARGET run-ci-tests
 COMMAND ${PY_VENV_DIR}/bin/python3 ${CI_TESTS_RUNTIME_DIR}/tests_runner.py --run-ci-tests
 COMMENT "Running continious integration tests"
)

add_custom_target(run-cpu-regression-tests)
add_dependencies(run-cpu-regression-tests ecrad create-python-virtual-environment)
add_custom_command(
 TARGET run-cpu-regression-tests
 COMMAND ${PY_VENV_DIR}/bin/python3 ${CI_TESTS_RUNTIME_DIR}/tests_runner.py --run-cpu-reg-tests
 COMMENT "Running continious integration CPU regression tests"
)




