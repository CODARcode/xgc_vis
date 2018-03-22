add_library (xgcvis_gui widget.cpp trackball.cpp)
set_target_properties (xgcvis_gui PROPERTIES AUTOMOC ON)
target_link_libraries (xgcvis_gui tourtre xgcvis_io ${GLUT_LIBRARY} ${QT_LIBRARIES} ${OPENGL_LIBRARIES} ${Boost_LIBRARIES})

add_executable (test_gui main.cpp)
target_link_libraries (test_gui xgcvis_gui)

# add_executable (xgcvis_gui main.cpp widget.cpp trackball.cpp)
# target_link_libraries (xgcvis_gui tourtre ${ADIOS_LIBRARIES} ${ZLIB_LIBRARIES} ${QT_LIBRARIES} ${GLUT_LIBRARY} ${GLEW_LIBRARY} ${OPENGL_LIBRARIES})
# target_link_libraries (xgcvis_gui tourtre xgcvis_io ${GLUT_LIBRARY} ${QT_LIBRARIES} ${OPENGL_LIBRARIES} ${Boost_LIBRARIES})
# set_target_properties (xgcvis_gui PROPERTIES AUTOMOC ON)
