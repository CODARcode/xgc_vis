set (volren_sources bvh.cpp volren.cu)
cuda_add_library (volren STATIC ${volren_sources})
