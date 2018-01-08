add_executable (xgcvis_gui main.cpp widget.cpp trackball.cpp)
target_link_libraries (xgcvis_gui tourtre ${ADIOS_LIBRARIES} ${ZLIB_LIBRARIES} ${QT_LIBRARIES} ${GLUT_LIBRARY} ${GLEW_LIBRARY} ${OPENGL_LIBRARIES})
set_target_properties (xgcvis_gui PROPERTIES AUTOMOC ON)
