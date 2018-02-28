add_library (xgcvis_io xgcMesh.cpp xgcData.cpp)
target_link_libraries (xgcvis_io ${ADIOS_LIBRARIES})

if (WITH_VTK)
  target_link_libraries (xgcvis_io ${VTK_LIBRARIES})

  add_executable (testio main.cpp)
  target_link_libraries (testio xgcvis_io)
endif ()
