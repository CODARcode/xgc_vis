add_executable (xgcvis_core main.cpp)
target_link_libraries (xgcvis_core tourtre ${ADIOS_LIBRARIES} ${ZLIB_LIBRARIES})
