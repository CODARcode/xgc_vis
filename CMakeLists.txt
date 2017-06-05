cmake_minimum_required (VERSION 3.1.3)
project (xgc_vis)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

find_package (OpenGL REQUIRED)
find_package (GLEW REQUIRED)
find_package (GLUT REQUIRED)
find_package (ZLIB REQUIRED)
find_package (Qt4 REQUIRED COMPONENTS QtCore QtGui QtOpenGL)
find_package (ADIOS REQUIRED)

include_directories (${ADIOS_INCLUDE_DIRS})

add_executable (main main.cpp)
target_link_libraries (main ${ADIOS_LIBRARIES} ${ZLIB_LIBRARIES})
