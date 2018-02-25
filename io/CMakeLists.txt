add_library (xgcvis_io xgcMesh.cpp xgcData.cpp)

if (WITH_VTK)
  target_link_libraries (xgcvis_io ${VTK_LIBRARIES})
endif ()
