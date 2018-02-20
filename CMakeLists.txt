cmake_minimum_required (VERSION 3.1.3)
project (xgc_vis)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

#set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_MULTITHREADED ON)
#set (Boost_USE_STATIC_RUNTIME OFF)

find_package (Boost REQUIRED COMPONENTS program_options system)
find_package (OpenGL REQUIRED)
find_package (GLEW REQUIRED)
find_package (GLUT REQUIRED)
find_package (ZLIB REQUIRED)
find_package (VTK REQUIRED)
find_package (Qt4 REQUIRED COMPONENTS QtCore QtGui QtOpenGL)
find_package (ADIOS REQUIRED)

find_package (CUDA REQUIRED)

# set (CUDA_NVCC_FLAGS "-arch=sm_60;-std=c++11")
set (CUDA_PROPAGATE_HOST_FLAGS OFF)

add_subdirectory (libtourtre)

include (${VTK_USE_FILE})
include (${QT_USE_FILE})
include_directories (
  ${ADIOS_INCLUDE_DIRS} 
  ${GLEW_INCLUDE_DIR} 
  ${PROJECT_SOURCE_DIR}/libtourtre/include
  ${PROJECT_SOURCE_DIR}/websocketpp
)

cuda_add_executable (main main.cpp widget.cpp trackball.cpp bvh.cpp volren.cu)
target_link_libraries (main tourtre ${Boost_LIBRARIES} ${ADIOS_LIBRARIES} ${ZLIB_LIBRARIES} ${QT_LIBRARIES} ${GLUT_LIBRARY} ${GLEW_LIBRARY} ${OPENGL_LIBRARIES})
set_target_properties (main PROPERTIES AUTOMOC ON)

add_executable (converter converter.cpp)
target_link_libraries (converter vtkIOCore vtkIOLegacy ${ADIOS_LIBRARIES} ${ZLIB_LIBRARIES})
