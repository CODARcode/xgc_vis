#include "def.h"
#include <QMouseEvent>
#include <QFileDialog>
#include <QInputDialog>
#include <QDebug>
#include <fstream>
#include <iostream>
#include <queue>
#include <json.hpp>
#include "widget.h"
#include "io/xgcMesh.h"
#include "io/xgcData.h"
#include "core/xgcLevelSetAnalysis.h"

#ifdef __APPLE__
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glu.h>
// #include <GL/glut.h>
#endif

using nlohmann::json;

typedef websocketpp::client<websocketpp::config::asio_client> client;
using websocketpp::lib::placeholders::_1;
using websocketpp::lib::placeholders::_2;
using websocketpp::lib::bind;

CGLWidget::CGLWidget(XGCMesh &m_, XGCData &d_, XGCDataReader &r_, ftk::tracking_graph<> &g_, std::vector<std::vector<std::set<size_t> > > &cc_,
      const QGLFormat& fmt, QWidget *parent, QGLWidget *sharedWidget)
  : m(m_), d(d_), r(r_), g(g_), cc(cc_),
    QGLWidget(fmt, parent, sharedWidget), 
    _fovy(30.f), _znear(0.1f), _zfar(10.f), 
    _eye(0, 0, 2.5), _center(0, 0, 0), _up(0, 1, 0), 
    toggle_mesh(false), toggle_wireframe(false), toggle_extrema(false), toggle_labels(true), 
    current_slice(0)
{
  // r.read(m, d);

  updateMesh();
  updateData();
  // thread_ws = new std::thread(&CGLWidget::connectToWebSocketServer, this, "ws://red:9002");

  // contour = d.sampleAlongPsiContour(m, 0.2); 
  // m.sampleScalarsAlongPsiContour(d.dpot, 10, 0.2); // TODO FIXME
  // contour = m.sampleScalarsAlongPsiContour(d.dpot, 10, 0.2); // TODO FIXME
  // contour = m.testMarchingTriangles(m.psi, 0.2);
}

CGLWidget::~CGLWidget()
{
  if (thread_ws != NULL) {
    thread_ws->join();
    delete thread_ws;
  }
}

void CGLWidget::onMessage(client *c, websocketpp::connection_hdl, message_ptr msg)
{

}

void CGLWidget::connectToWebSocketServer(const std::string& uri)
{
  try {
    // Set logging to be pretty verbose (everything except message payloads)
    c.set_access_channels(websocketpp::log::alevel::all);
    c.clear_access_channels(websocketpp::log::alevel::frame_payload);

    // Initialize ASIO
    c.init_asio();

    // Register our message handler
    c.set_message_handler(bind(&CGLWidget::onMessage, this, &c, ::_1, ::_2));

    websocketpp::lib::error_code ec;
    client::connection_ptr con = c.get_connection(uri, ec);
    if (ec) {
      std::cout << "could not create connection because: " << ec.message() << std::endl;
      // return 0;
    }

    // Note that connect here only requests a connection. No network messages are
    // exchanged until the event loop starts running in the next line.
    c.connect(con);

    // Start the ASIO io_service run loop
    // this will cause a single connection to be made to the server. c.run()
    // will exit when this connection is closed.
    c.run();
  } catch (websocketpp::exception const & e) {
    std::cout << e.what() << std::endl;
  }
}

void CGLWidget::loadBranchesFromJsonFile(const std::string& filename)
{
  json j;
  std::ifstream ifs(filename, std::ifstream::in);
  ifs >> j;
  ifs.close();

  for (int i=0; i<j.size(); i++) {
    branches.push_back(new ctBranch);
    label_colors[i] = QColor(rand()%256, rand()%256, rand()%256);
  }

  for (json::iterator it = j.begin(); it != j.end(); it ++) {
    json jbranch = *it;
    int id = jbranch["id"];
    ctBranch *b = branches[id];

    b->extremum = jbranch["extremum"];
    b->saddle = jbranch["saddle"];

    if (!jbranch["children"].is_null())
      b->children.head = branches[jbranch["children"][0]];

    if (!jbranch["nextChild"].is_null())
      b->nextChild = branches[jbranch["nextChild"]];
    
    if (!jbranch["prevChild"].is_null())
      b->prevChild = branches[jbranch["prevChild"]];
  }

  fprintf(stderr, "Branch decompositions loaded: %lu\n", branches.size());
}

void CGLWidget::mousePressEvent(QMouseEvent* e)
{
  _trackball.mouse_rotate(e->x(), e->y()); 
}

void CGLWidget::mouseMoveEvent(QMouseEvent* e)
{
  _trackball.motion_rotate(e->x(), e->y()); 
  updateGL(); 
}

void CGLWidget::keyPressEvent(QKeyEvent* e)
{
  switch (e->key()) {
  case Qt::Key_M: 
    toggle_mesh = !toggle_mesh;
    updateGL();
    break;

  case Qt::Key_W:
    toggle_wireframe = !toggle_wireframe; 
    updateGL();
    break;

#if 0
  case Qt::Key_E:
    toggle_extrema = !toggle_extrema;
    updateGL();
    break;
#endif
  case Qt::Key_L:
    toggle_labels = !toggle_labels;
    updateGL();
    break;

  case Qt::Key_Up:
    if (r.advanceTimestep() >= 0) {
      fprintf(stderr, "advancing timestep... %zu\n", r.getCurrentTimestep());
      r.read(m, d);
      updateData();
      updateGL();
    }
    break;

  case Qt::Key_Down:
    if (r.recedeTimestep() >= 0) {
      fprintf(stderr, "receding timestep... %zu\n", r.getCurrentTimestep());
      r.read(m, d);
      updateData();
      updateGL();
    }
    break;

  case Qt::Key_Right:
    current_slice = (current_slice + 1) % m.nPhi;
    updateGL();
    break;

  case Qt::Key_Left:
    current_slice = (current_slice - 1 + m.nPhi) % m.nPhi;
    updateGL();
    break;

  default: break;
  }
}

void CGLWidget::wheelEvent(QWheelEvent* e)
{
  _trackball.wheel(e->delta());
  updateGL(); 
}

void CGLWidget::initializeGL()
{
  // glewInit();
  _trackball.init();

  // opengl smooth rendering
  {
    glEnable(GL_MULTISAMPLE);

    GLint bufs, samples; 
    glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs); 
    glGetIntegerv(GL_SAMPLES, &samples); 

    glEnable(GL_LINE_SMOOTH); 
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST); 
    
    glEnable(GL_POLYGON_SMOOTH); 
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST); 
    
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1, 1);

    glEnable(GL_BLEND); 
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
  
  // initialze light for tubes
  {
    GLfloat ambient[]  = {0.1, 0.1, 0.1}, 
            diffuse[]  = {0.5, 0.5, 0.5}, 
            specular[] = {0.8, 0.8, 0.8}; 
    GLfloat dir[] = {0, 0, -1}; 
    GLfloat pos[] = {1, 1, 4, 1};
    GLfloat shiness = 100; 

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient); 
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse); 
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular); 
    glLightfv(GL_LIGHT0, GL_POSITION, pos); 
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir); 
    
    GLfloat light1_position[] = {-4.0, 4.0, 0.0, 1.0};
    GLfloat light1_spot_direction[] = {1.0, -1.0, 0.0};

    glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, light1_spot_direction);

    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE); 
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE); 

    glEnable(GL_NORMALIZE); 
    glEnable(GL_COLOR_MATERIAL); 
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular); 
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shiness); 
  }

  CHECK_GLERROR();
}

void CGLWidget::resizeGL(int w, int h)
{
  _trackball.reshape(w, h); 
  glViewport(0, 0, w, h);

  CHECK_GLERROR(); 
}

void CGLWidget::paintGL()
{
  glClearColor(1, 1, 1, 0); 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

  _projmatrix.setToIdentity(); 
  _projmatrix.perspective(_fovy, (float)width()/height(), _znear, _zfar); 
  _mvmatrix.setToIdentity();
  _mvmatrix.lookAt(_eye, _center, _up);
  _mvmatrix.rotate(_trackball.getRotation());
  _mvmatrix.scale(_trackball.getScale());

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glLoadMatrixd(_projmatrix.data()); 
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  glLoadMatrixd(_mvmatrix.data()); 

  // glutWireTeapot(1.0);

  renderSinglePlane();
  // renderMultiplePlanes();

#if 0
  glTranslatef(-m.coords_centroid_x, -m.coords_centroid_y, 0.f);
  glColor3f(0, 0, 0);
  glBegin(GL_LINE_LOOP);
  // glBegin(GL_POINTS);
  for (int i=0; i<contour.size()/3; i++) 
    glVertex2f(contour[i*3], contour[i*3+1]);
  glEnd();
#endif

  CHECK_GLERROR();
}

void CGLWidget::renderMultiplePlanes()
{
  if (f_vertices.size() == 0) return;

  glColor3f(0, 0, 0);
 
  if (toggle_wireframe) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // render wireframe
  }
  else {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

  glEnable(GL_DEPTH_TEST);
  // glDisable(GL_DEPTH_TEST);
  // glBlendFunc(GL_SRC_ALPHA, GL_ONE);

  for (int i=0; i<m.nPhi; i++) {
    glPushMatrix();
    glRotatef(360.f/m.nPhi*i, 0, 1, 0);

#if 0
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    glVertexPointer(2, GL_FLOAT, 0, f_vertices.data());
    glColorPointer(3, GL_FLOAT, 0, f_colors.data());
    glDrawArrays(GL_TRIANGLES, 0, f_vertices.size()/2);

    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
#endif

    glColor4f(1, 0, 0, 0.8);
    glBegin(GL_POINTS);
    for (int j=0; j<maximum[i].size(); j++) {
      int k = maximum[i][j];
      glVertex2f(m.coords[k*2], m.coords[k*2+1]);
    }
    glEnd();
    
    glColor4f(0, 1, 0, 0.8);
    glBegin(GL_POINTS);
    for (int j=0; j<minimum[i].size(); j++) {
      int k = minimum[i][j];
      glVertex2f(m.coords[k*2], m.coords[k*2+1]);
    }
    glEnd();

    glPopMatrix();
  }
}

void CGLWidget::renderSinglePlane()
{
  if (f_vertices.size() == 0) return;

  glColor3f(0, 0, 0);

  glTranslatef(-1.7, 0, 0);
    
  if (toggle_labels) {
    glPointSize(5.f);
    glBegin(GL_POINTS);
    for (int i=0; i<m.nNodes; i++) {
      int label = labels[i];
      if (label != 0) {
        QColor c = label_colors[label];
        glColor3ub(c.red(), c.green(), c.blue());
        glVertex2f(m.coords[i*2], m.coords[i*2+1]);
      }
    }
    glEnd();
  } else if (toggle_mesh) {
    if (toggle_wireframe) {
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // render wireframe
    }
    else {
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    glVertexPointer(2, GL_FLOAT, 0, f_vertices.data());
    glColorPointer(3, GL_FLOAT, 0, f_colors.data() );
    glDrawArrays(GL_TRIANGLES, 0, f_vertices.size()/2);

    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
  }

  CHECK_GLERROR();
}

template <typename T>
static T clamp(T min, T max, T val)
{
  return std::max(min, std::min(max, val));
}

template <typename T>
static T clamp_normalize(T min, T max, T val)
{
  return (clamp(min, max, val) - min) / (max - min);
}
 
void CGLWidget::updateMesh()
{
  f_vertices.clear();

  for (int i=0; i<m.nTriangles; i++) {
    int i0 = m.conn[i*3], i1 = m.conn[i*3+1], i2 = m.conn[i*3+2];
    f_vertices.push_back(m.coords[i0*2]);
    f_vertices.push_back(m.coords[i0*2+1]);
    f_vertices.push_back(m.coords[i1*2]);
    f_vertices.push_back(m.coords[i1*2+1]);
    f_vertices.push_back(m.coords[i2*2]);
    f_vertices.push_back(m.coords[i2*2+1]);
  }
}

void CGLWidget::updateData()
{
  fprintf(stderr, "update data...\n");
#if 1
  const double threshold = 80;

  auto components = cc[r.getCurrentTimestep()];

  labels.resize(m.nNodes * m.nPhi);
  std::fill(labels.begin(), labels.end(), 0);
  for (size_t i = 0; i < components.size(); i++) {
    auto gLabel = g.get_global_label(r.getCurrentTimestep(), i);
    if (label_colors.find(gLabel+1) == label_colors.end())
      label_colors[gLabel+1] = QColor(rand()%256, rand()%256, rand()%256);
    // qDebug() << label_colors[i];
    for (const auto v : components[i]) {
      if (v <= m.nNodes)
        labels[v] = gLabel+1; 
    }
  }
#endif

#if 1
  // const float min = -100, max = 100;
  const float min = -1, max = 1;

  f_colors.clear();

  for (int plane = 0; plane < m.nPhi; plane ++) {
    for (int i=0; i<m.nTriangles; i++) {
      int v[3] = {m.conn[i*3], m.conn[i*3+1], m.conn[i*3+2]};

      for (int j=0; j<3; j++) {
        // float val = clamp_normalize(min, max, (float)d.dpot[plane*m.nNodes + v[j]]);
        float val = d.dneOverne0[plane*m.nNodes + v[j]];
        f_colors.push_back(val);
        f_colors.push_back(1-val);
        f_colors.push_back(0);
      }
    }
  }
#endif
}
