add_executable (xgcvis_gui2D main.cpp widget.cpp trackball.cpp)
set_target_properties (xgcvis_gui2D PROPERTIES AUTOMOC ON)
target_link_libraries (xgcvis_gui2D ${GLUT_LIBRARY} ${QT_LIBRARIES} ${OPENGL_LIBRARIES} ${Boost_LIBRARIES} ${GLEW_LIBRARY})
