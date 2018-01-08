add_executable (converter converter.cpp)
target_link_libraries (converter vtkIOCore vtkIOLegacy ${ADIOS_LIBRARIES} ${ZLIB_LIBRARIES})
