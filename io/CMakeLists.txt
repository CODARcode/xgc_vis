add_library (xgcvis_io xgcMesh.cpp xgcData.cpp xgcEq.cpp)
target_link_libraries (xgcvis_io xgcvis_common ${ADIOS_LIBRARIES} ${ZLIB_LIBRARY})
  
add_executable (testio main.cpp)
  target_link_libraries (testio xgcvis_io)

if (WITH_VTK)
  target_link_libraries (xgcvis_io ${VTK_LIBRARIES})
endif ()

if (WITH_NC)
  target_link_libraries (xgcvis_io ${NETCDF_LIBRARIES})
endif ()

if (WITH_H5)
  target_link_libraries (xgcvis_io ${HDF5_C_LIBRARIES})
endif ()
