#ifndef _WIDGET_H
#define _WIDGET_H

#include "io/xgcMesh.h"
#include "io/xgcData.h"
#include "io/xgcDataReader.h"

// #include <GL/glew.h>
#include <QGLWidget>
#include <QList>
#include <QVector>
#include <QVector3D>
#include <QMatrix4x4>
#include <tourtre.h>
#include <cmath>
#include <set>
#include <ftk/tracking_graph.hh>
#include "trackball.h"

#ifndef Q_MOC_RUN
#include <websocketpp/config/asio_no_tls_client.hpp>
#include <websocketpp/client.hpp>
#endif

class QMouseEvent;
class QKeyEvent; 
class QWheelEvent; 

/* 
 * \class   CGLWidget
 * \author  Hanqi Guo
 * \brief   A light-weight Qt-based vortex viewer
*/
class CGLWidget : public QGLWidget
{
  Q_OBJECT

public:
  CGLWidget(XGCMesh &m, XGCData &d, XGCDataReader &r, ftk::tracking_graph<> &g, std::vector<std::vector<std::set<size_t> > > &cc,
      const QGLFormat& fmt=QGLFormat::defaultFormat(), QWidget *parent=NULL, QGLWidget *sharedWidget=NULL); 
  ~CGLWidget(); 

  void updateMesh();
  void updateData();

  void loadBranchesFromJsonFile(const std::string& filename);

protected:
  void initializeGL(); 
  void resizeGL(int w, int h); 
  void paintGL();

  void mousePressEvent(QMouseEvent*); 
  void mouseMoveEvent(QMouseEvent*);
  void keyPressEvent(QKeyEvent*); 
  void wheelEvent(QWheelEvent*); 

protected:
  void renderSinglePlane();
  void renderMultiplePlanes();
  void renderExtremum();
  void renderLabels();

public: // client
  std::thread *thread_ws;
  void connectToWebSocketServer(const std::string& uri);
  
  typedef websocketpp::client<websocketpp::config::asio_client> client;
  typedef websocketpp::config::asio_client::message_type::ptr message_ptr;
 
  client c;
  void onMessage(client* c, websocketpp::connection_hdl hdl, message_ptr msg);

private: 
  XGCMesh &m; 
  XGCData &d;
  XGCDataReader &r;

  ftk::tracking_graph<> &g;
  std::vector<std::vector<std::set<size_t> > > &cc;
  // int currentTimestep = 10;

  std::vector<double> contour;

private:
  CGLTrackball _trackball;
  QMatrix4x4 _projmatrix, _mvmatrix; 

private: // camera
  const float _fovy, _znear, _zfar; 
  const QVector3D _eye, _center, _up;

  bool toggle_mesh, toggle_wireframe, toggle_extrema, toggle_labels;
  int current_slice;

private: // mesh
  std::vector<float> f_vertices;
  std::vector<float> f_colors, f_label_colors;

private: // analysis
  std::vector<std::set<int> > nodeGraph; // node->{neighbor nodes}
  std::map<int, std::vector<int> > maximum, minimum;
  std::vector<size_t> labels;
  std::map<size_t, QColor> label_colors; 
  // std::map<int, std::vector<size_t> > all_labels;
  // std::vector<int> maximum, minimum;
  
  std::vector<ctBranch*> branches;
  std::map<ctBranch*, size_t> branchSet;
}; 

#endif
