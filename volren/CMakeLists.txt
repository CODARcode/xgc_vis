set (volren_sources bvh.cpp volrenEngine.cpp)
set (volren_cuda_sources volren.cu)
set (volren_cuda_cxx_sources volren.cu.cpp)

if (WITH_CUDA)
  cuda_add_library (xgcvis_volren SHARED ${volren_sources} ${volren_cuda_sources})
  # cuda_add_executable (test_volren main.cpp)
else ()
  add_library (xgcvis_volren STATIC ${volren_sources} ${volren_cuda_cxx_sources})
endif ()

target_link_libraries (xgcvis_volren xgcvis_io ${PNG_LIBRARIES})
# target_link_libraries (test_volren volren)
