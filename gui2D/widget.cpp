#include "def.h"
#include <QMouseEvent>
#include <QFileDialog>
#include <QInputDialog>
#include <QDebug>
#include <fstream>
#include <iostream>
#include <queue>
#include <functional>
#include <hdf5.h>
// #include <json.hpp>
#include "widget.h"

#ifdef __APPLE__
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glu.h>
// #include <GL/glut.h>
#endif

CGLWidget::CGLWidget(const QGLFormat& fmt, QWidget *parent, QGLWidget *sharedWidget) :
  _fovy(30.f), _znear(0.1f), _zfar(10.f), 
  _eye(0, 0, 2.5), _center(0, 0, 0), _up(0, 1, 0)
{
}

CGLWidget::~CGLWidget()
{
}

int CGLWidget::advanceSlice() {
  currentSlice ++;
  if (currentSlice >= nphi) currentSlice -= nphi;
  return currentSlice;
}

int CGLWidget::recedeSlice() {
  currentSlice --;
  if (currentSlice < 0) currentSlice += nphi;
  return currentSlice;
}

int CGLWidget::advanceTimestep() {
  if (currentTimestep < nt-1) currentTimestep ++;
  return currentTimestep;
}

int CGLWidget::recedeTimestep() {
  if (currentTimestep > 0) currentTimestep --;
  return currentTimestep;
}
  
void CGLWidget::trackSuperLevelsetT()
{
  typedef std::vector<std::set<size_t> > Components;

  std::shared_ptr<Components> components = 
    std::make_shared<Components>(ftk::extractConnectedComponents<size_t, std::set<size_t> >(
          W*H*nt,
          std::bind(ftk::Get26Neighbors3DRegular<size_t>, W, H, nt, std::placeholders::_1), 
          [this](size_t i) {return data[i] >= 0.2;})); // FIXME

  fprintf(stderr, "#components=%zu\n", components->size());

  labelsRGB.resize(W*H*nt*3);
  for (auto c : *components) {
    if (c.size() < 100) continue;
    GLubyte r = rand()%256, g = rand()%256, b = rand()%256;
    for (auto i : c) {
      labelsRGB[i*3] = r;
      labelsRGB[i*3+1] = g;
      labelsRGB[i*3+2] = b;
    }
  }
}

void CGLWidget::trackSuperLevelset()
{
#if 0
  typedef std::vector<std::set<size_t> > Components;

  std::shared_ptr<Components> components = 
    std::make_shared<Components>(ftk::extractConnectedComponents<size_t, std::set<size_t> >(
          W*H,
          std::bind(ftk::Get8Neighbors2DRegular<size_t>, W, H, std::placeholders::_1), 
          [this](size_t i) {return data[currentTimestep][W*H*currentSlice+i] >= 0.2;}));

  fprintf(stderr, "#components=%zu\n", components->size());

  std::vector<float> d(W*H*4);
  for (auto c : *components) {
    // if (c.size() < 20) continue;
    float r = (float)rand() / RAND_MAX, 
          g = (float)rand() / RAND_MAX, 
          b = (float)rand() / RAND_MAX;
    for (auto i : c) {
      d[i*4] = 1; r; 
      d[i*4+1] = 1; g;
      d[i*4+2] = 1; b;
      d[i*4+3] = 1;
    }
  }
  
  glBindTexture(GL_TEXTURE_2D, texLabels);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, W, H, 0, GL_RGBA, GL_FLOAT, &d[0]);
#endif
}

void CGLWidget::loadData(const std::string& path)
{
  data.resize(W*H*nphi*nt);

  for (int i=0; i<nt; i++) { // time
    for (int j=0; j<nphi; j++) { // phi
      const std::string filename = path + "/norm-" + std::to_string(i+1) + "-" + std::to_string(j) + ".dat";
      FILE *fp = fopen(filename.c_str(), "rb");
      assert(fp);

      std::vector<double> dd(W*H);
      fread(&dd[0], sizeof(double), W*H, fp);

      for (int k=0; k<W*H; k++) {
        data[W*H*nphi*i + W*H*j + k] = dd[k];
      }
      fclose(fp);
    }
  }

  trackSuperLevelsetT();
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
  case Qt::Key_Up: 
    advanceTimestep();
    fprintf(stderr, "%d.%d\n", currentTimestep, currentSlice);
    updateDataGL();
    updateGL();
    break;

  case Qt::Key_Down:
    recedeTimestep();
    fprintf(stderr, "%d.%d\n", currentTimestep, currentSlice);
    updateDataGL();
    updateGL();
    break;

  case Qt::Key_Left:
    recedeSlice();
    fprintf(stderr, "%d.%d\n", currentTimestep, currentSlice);
    updateDataGL();
    updateGL();
    break;

  case Qt::Key_Right:
    advanceSlice();
    fprintf(stderr, "%d.%d\n", currentTimestep, currentSlice);
    updateDataGL();
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

void CGLWidget::updateDataGL()
{
  // trackSuperLevelset();

  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, W, H, 0, GL_LUMINANCE, GL_FLOAT, 
      &data[W*H*nphi*currentTimestep + W*H*currentSlice]);
  // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, W, H, 0, GL_RGBA, GL_FLOAT, &d[0]);
  
  glBindTexture(GL_TEXTURE_2D, texLabels);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, W, H, 0, GL_RGB, GL_UNSIGNED_BYTE, 
      &labelsRGB[3*(W*H*nphi*currentTimestep + W*H*currentSlice)]);
}

void CGLWidget::initializeGL()
{
  glewInit();
  
  glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

  glGenTextures(1, &texLabels);
  glBindTexture(GL_TEXTURE_2D, texLabels);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

  updateDataGL();

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

  glScalef(0.5004, 1, 1);
  glRotatef(-90, 0, 0, 1);
  glTranslatef(-0.5, -0.5, 0);
  glColor4f(1, 1, 1, 1);
  // glBindTexture(GL_TEXTURE_2D, tex);
  glBindTexture(GL_TEXTURE_2D, texLabels);
  glEnable(GL_TEXTURE_2D);
  glBegin(GL_QUADS);
  glVertex2f(0, 0); glTexCoord2f(0, 0);
  glVertex2f(1, 0); glTexCoord2f(1, 0);
  glVertex2f(1, 1); glTexCoord2f(1, 1);
  glVertex2f(0, 1); glTexCoord2f(0, 1);
  glEnd();
  glDisable(GL_TEXTURE_2D);

  CHECK_GLERROR();
}

