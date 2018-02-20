#ifndef _WIDGET_H
#define _WIDGET_H

#include <GL/glew.h>
#include <QGLWidget>
#include <QList>
#include <QVector>
#include <QVector3D>
#include <QMatrix4x4>
#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>
#include <tourtre.h>
#include <cmath>
#include <set>
#include "trackball.h"

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
  CGLWidget(const QGLFormat& fmt=QGLFormat::defaultFormat(), QWidget *parent=NULL, QGLWidget *sharedWidget=NULL); 
  ~CGLWidget(); 

  void setTriangularMesh(int nNodes, int nTriangles, int nPhi, double *coords, int *conn);
  void setData(double *dpot);
  void updateDpot(double *dpot);

  void setNLimit(size_t n) {nLimit = n;}

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

private: 
  void buildContourTree(int plane, double *dpot);
  void addExtremumFromBranchDecomposition(int plane, ctBranch *root, ctBranch *b, void *);
  void simplifyBranchDecompositionByThreshold(ctBranch *b, double threshold, void *);
  std::set<ctBranch*> simplifyBranchDecompositionByNumbers(ctBranch* root, std::map<ctBranch*, size_t>&, int nLimit, void *);
  void buildSegmentation(ctBranch *b, std::vector<size_t> &labels, void*);

  void buildContourTree3D(double *dpot); 
  void buildSegmentation3D(ctBranch *b, std::vector<size_t> &labels, void*); 

  void extractExtremum(int plane, double *dpot);
  void constructDiscreteGradient(double *dpot);

  void extractStreamers(int plane, ctBranch *root, std::map<ctBranch*, size_t>& branchSet, int nStreamers, double percentage, void *);
  void extractStreamersFromExtremum(int plane, double *dpot, double percentage=0.1);
  int flood2D(size_t seed, size_t id, std::vector<size_t> &labels, double min, double max, void *d);

private: // server
  typedef websocketpp::server<websocketpp::config::asio> server;
  typedef server::message_ptr message_ptr;
  server wss;
  std::thread *thread_wss;

public:
  void startServer(int port);
  void onMessage(server* s, websocketpp::connection_hdl hdl, message_ptr msg);

private:
  CGLTrackball _trackball;
  QMatrix4x4 _projmatrix, _mvmatrix; 

private: // camera
  const float _fovy, _znear, _zfar; 
  const QVector3D _eye, _center, _up;

  bool toggle_mesh, toggle_wireframe, toggle_extrema, toggle_labels;
  int current_slice;

private: // volren
  struct ctx_rc *rc;
  float *framebuf;

private: // mesh
  double *coords; 
  int *conn;
  int nNodes, nTriangles, nPhi;

  std::vector<float> f_vertices;
  std::vector<float> f_colors;

  std::map<size_t, QColor> label_colors;

private: // analysis
  size_t nLimit;

  std::vector<std::set<int> > nodeGraph; // node->{neighbor nodes}
  std::map<int, std::vector<int> > maximum, minimum;
  std::map<int, std::vector<size_t> > all_labels;
  // std::vector<int> maximum, minimum;
}; 

#endif
