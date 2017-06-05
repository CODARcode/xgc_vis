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

include (${QT_USE_FILE})
include_directories (${ADIOS_INCLUDE_DIRS} ${GLEW_INCLUDE_DIR})

add_executable (main main.cpp widget.cpp trackball.cpp)
target_link_libraries (main ${ADIOS_LIBRARIES} ${ZLIB_LIBRARIES} ${QT_LIBRARIES} ${GLUT_LIBRARY} ${GLEW_LIBRARY} ${OPENGL_LIBRARIES})
set_target_properties (main PROPERTIES AUTOMOC ON)
