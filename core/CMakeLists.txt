add_executable (xgcvis_core main.cpp xgcBlobExtractor.cpp)
if (XGCVIS_BUILD_VOLREN)
  target_link_libraries (xgcvis_core tourtre common volren ${ADIOS_LIBRARIES} ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${PNG_LIBRARIES})
else ()
  target_link_libraries (xgcvis_core tourtre common ${ADIOS_LIBRARIES} ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})
endif()

#add_executable (test1 test.cpp xgcBlobExtractor.cpp)
#target_link_libraries (test1 tourtre common ${ADIOS_LIBRARIES} ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})

#add_executable (wsproxy wsproxy.cpp)
#target_link_libraries (wsproxy ${Boost_LIBRARIES})
