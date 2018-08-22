add_library (xgcvis_core xgcMultiscaleAnalysis.cpp xgcLevelSetAnalysis.cpp)
target_link_libraries (xgcvis_core tourtre xgcvis_common xgcvis_volren ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable (xgcvis main.cpp)
target_link_libraries (xgcvis xgcvis_core)

if (WITH_ROCKSDB)
  target_link_libraries (xgcvis ${RocksDB_LIBRARY})
endif ()

if (WITH_QT)
  target_link_libraries (xgcvis xgcvis_gui)
endif ()

# add_executable (test1 test.cpp xgcBlobExtractor.cpp xgcLevelSetAnalysis.cpp)
# target_link_libraries (test1 tourtre xgcvis_common ${ADIOS_LIBRARIES} ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})

#add_executable (wsproxy wsproxy.cpp)
#target_link_libraries (wsproxy ${Boost_LIBRARIES})
