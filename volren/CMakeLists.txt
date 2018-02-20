set (volren_sources bvh.cpp volren.cu)
cuda_add_library (volren SHARED ${volren_sources})

cuda_add_executable (test_volren main.cpp)
target_link_libraries (test_volren volren)
