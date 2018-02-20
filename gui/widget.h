#ifndef _WIDGET_H
#define _WIDGET_H

// #include <GL/glew.h>
#include <QGLWidget>
#include <QList>
#include <QVector>
#include <QVector3D>
#include <QMatrix4x4>
#include <tourtre.h>
#include <cmath>
#include <set>
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
  CGLWidget(const QGLFormat& fmt=QGLFormat::defaultFormat(), QWidget *parent=NULL, QGLWidget *sharedWidget=NULL); 
  ~CGLWidget(); 

  void setTriangularMesh(int nNodes, int nTriangles, int nPhi, double *coords, int *conn);
  void loadMeshFromJsonFile(const std::string& filename);
  void loadBranchesFromJsonFile(const std::string& filename);
  void loadLabels(const std::string& filename);
  // void setData(double *dpot);

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
  void simplifyBranchDecompositionByThreshold(ctBranch *b, double threshold, void *);
  void simplifyBranchDecompositionByNumbers(ctBranch *b, int nLimit, void *);
  void buildSegmentation(ctBranch *b, std::vector<size_t> &labels, void*); 

  void buildContourTree3D(double *dpot); 
  void buildSegmentation3D(ctBranch *b, std::vector<size_t> &labels, void*); 

  void extractExtremum(int plane, double *dpot);
  void constructDiscreteGradient(double *dpot);

public: // client
  void connectToWebSocketServer(const std::string& uri);
  
  typedef websocketpp::client<websocketpp::config::asio_client> client;
  typedef websocketpp::config::asio_client::message_type::ptr message_ptr;
  void onMessage(client* c, websocketpp::connection_hdl hdl, message_ptr msg);

private:
  CGLTrackball _trackball;
  QMatrix4x4 _projmatrix, _mvmatrix; 

private: // camera
  const float _fovy, _znear, _zfar; 
  const QVector3D _eye, _center, _up;

  bool toggle_mesh, toggle_wireframe, toggle_extrema, toggle_labels;
  int current_slice;

private: // mesh
  std::vector<double> coords; 
  std::vector<int> conn;
  int nNodes, nTriangles, nPhi;

  std::vector<float> f_vertices;
  std::vector<float> f_colors;

  std::map<size_t, QColor> label_colors;

private: // analysis
  std::vector<std::set<int> > nodeGraph; // node->{neighbor nodes}
  std::map<int, std::vector<int> > maximum, minimum;
  std::map<int, std::vector<size_t> > all_labels;
  // std::vector<int> maximum, minimum;
  
  std::vector<ctBranch*> branches;
  std::map<ctBranch*, size_t> branchSet;
}; 

#endif
