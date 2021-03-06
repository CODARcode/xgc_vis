cmake_minimum_required (VERSION 3.1.3)
project (xgc_vis)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

configure_file (def.h.cmake def.h)

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

option (WITH_ADIOS "Build with ADIOS" ON)
option (WITH_CUDA "Build with CUDA" OFF)
option (WITH_QT "Build with Qt" OFF)
option (WITH_VTK "Build with VTK" OFF)
option (WITH_PNG "Build with PNG" OFF)
option (WITH_NC "Build with NetCDF" OFF)
option (WITH_H5 "Build with HDF5" OFF)
option (WITH_ROCKSDB "Build with RocksDB" OFF)
option (WITH_OPENMP "Build with OpenMP" OFF)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_RUNTIME OFF)

find_package (FTK REQUIRED)
find_package (ZLIB REQUIRED)
find_package (Boost REQUIRED COMPONENTS program_options system)

set (CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories (
  ${PROJECT_BINARY_DIR}
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/third_party/libtourtre/include
  ${PROJECT_SOURCE_DIR}/third_party/ftk/include
  ${PROJECT_SOURCE_DIR}/third_party/json
  ${PROJECT_SOURCE_DIR}/third_party/websocketpp
  ${PROJECT_SOURCE_DIR}/third_party/concurrentqueue
  ${PROJECT_SOURCE_DIR}/third_party/cutil
  ${Boost_INCLUDE_DIR}
  ${FTK_INCLUDE_DIR}
)

if (WITH_OPENMP)
  find_package (OpenMP REQUIRED)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

if (WITH_CUDA)
  find_package (CUDA REQUIRED)
  set (CUDA_PROPAGATE_HOST_FLAGS OFF)
endif ()

if (WITH_PNG)
  find_package (PNG REQUIRED)
  include_directories (${PNG_INCLUDE_DIRS})
endif ()

if (WITH_NC)
  find_package (NetCDF REQUIRED)
  include_directories (${NETCDF_INCLUDES})
endif ()

if (WITH_ADIOS)
  find_package (ADIOS REQUIRED)
  include_directories (${ADIOS_INCLUDE_DIRS})
endif ()

if (WITH_H5)
  find_package (HDF5 REQUIRED)
  include_directories (${HDF5_INCLUDE_DIRS})
endif ()

if (WITH_ROCKSDB)
  find_package (RocksDB REQUIRED)
  include_directories (${RocksDB_INCLUDE_DIR})
endif ()

if (WITH_QT)
  find_package (OpenGL REQUIRED)
  find_package (GLUT REQUIRED)
  find_package (GLEW REQUIRED)
  find_package (Qt4 REQUIRED COMPONENTS QtCore QtGui QtOpenGL)
  include (${QT_USE_FILE})
  include_directories (${GLEW_INCLUDE_DIR})
  add_subdirectory (gui)
  add_subdirectory (gui2D)
endif ()

if (WITH_VTK)
  find_package (VTK REQUIRED)
  include (${VTK_USE_FILE})
  add_subdirectory (converter)
endif ()

add_subdirectory (third_party/libtourtre)
add_subdirectory (core)
add_subdirectory (common)
add_subdirectory (volren)
add_subdirectory (io)
