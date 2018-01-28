cmake_minimum_required (VERSION 3.1.3)
project (xgc_vis)

option (XGCVIS_BUILD_CORE "Build core analysis code" ON)
option (XGCVIS_BUILD_GUI "Build GUI" OFF)
option (XGCVIS_BUILD_CONVERTER "Build VTK converter" OFF)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

find_package (ZLIB REQUIRED)
set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_RUNTIME OFF)
find_package (Boost REQUIRED COMPONENTS program_options system)

include_directories (
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/third_party/libtourtre/include
  ${PROJECT_SOURCE_DIR}/third_party/json
  ${PROJECT_SOURCE_DIR}/third_party/websocketpp
  ${Boost_INCLUDE_DIR}
)

if (XGCVIS_BUILD_CORE)
  find_package (ADIOS REQUIRED)
  include_directories (
    ${ADIOS_INCLUDE_DIRS} 
  )
  add_subdirectory (core)
endif ()

if (XGCVIS_BUILD_GUI)
  find_package (OpenGL REQUIRED)
  # find_package (GLEW REQUIRED)
  # find_package (GLUT REQUIRED)
  find_package (Qt4 REQUIRED COMPONENTS QtCore QtGui QtOpenGL)
  include (${QT_USE_FILE})
  # include_directories (${GLEW_INCLUDE_DIR})
  add_subdirectory (gui)
endif ()

if (XGCVIS_BUILD_CONVERTER)
  find_package (VTK REQUIRED)
  include (${VTK_USE_FILE})
  add_subdirectory (converter)
endif ()

add_subdirectory (third_party/libtourtre)
