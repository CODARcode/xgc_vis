add_executable (xgcvis_core main.cpp xgcBlobExtractor.cpp)
target_link_libraries (xgcvis_core tourtre ${ADIOS_LIBRARIES} ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})

#add_executable (wsproxy wsproxy.cpp)
#target_link_libraries (wsproxy ${Boost_LIBRARIES})
