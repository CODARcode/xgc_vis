#ifndef _WIDGET_H
#define _WIDGET_H

#include <GL/glew.h>
#include <QGLWidget>
#include <QList>
#include <QVector>
#include <QVector3D>
#include <QMatrix4x4>
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

private: 
  void extractExtremum(int plane, double *dpot);
  void constructDiscreteGradient(double *dpot);

private:
  CGLTrackball _trackball;
  QMatrix4x4 _projmatrix, _mvmatrix; 

private: // camera
  const float _fovy, _znear, _zfar; 
  const QVector3D _eye, _center, _up;

  bool toggle_wireframe, toggle_extrema;

private: // mesh
  double *coords; 
  int *conn;
  int nNodes, nTriangles, nPhi;

  std::vector<float> f_vertices;
  std::vector<float> f_colors;

private: // analysis
  std::vector<std::set<int> > nodeGraph; // node->{neighbor nodes}
  std::map<int, std::vector<int> > maximum, minimum;
  // std::vector<int> maximum, minimum;
}; 

#endif
