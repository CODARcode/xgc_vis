add_executable (xgcvis_core main.cpp xgcBlobExtractor.cpp xgcMultiscaleAnalysis.cpp xgcLevelSetAnalysis.cpp)
target_link_libraries (xgcvis_core tourtre xgcvis_common xgcvis_volren ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})

#add_executable (test1 test.cpp xgcBlobExtractor.cpp)
#target_link_libraries (test1 tourtre common ${ADIOS_LIBRARIES} ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})

#add_executable (wsproxy wsproxy.cpp)
#target_link_libraries (wsproxy ${Boost_LIBRARIES})
