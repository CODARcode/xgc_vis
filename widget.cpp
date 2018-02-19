#include <QMouseEvent>
#include <QFileDialog>
#include <QInputDialog>
#include <QDebug>
#include <fstream>
#include <iostream>
#include <queue>
#include "widget.h"
#include "volren.cuh"

#ifdef __APPLE__
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#define CHECK_GLERROR()\
{\
  GLenum err = glGetError();\
  if (err != GL_NO_ERROR) {\
    const GLubyte *errString = gluErrorString(err);\
    qDebug("[%s line %d] GL Error: %s\n",\
            __FILE__, __LINE__, errString);\
  }\
}

CGLWidget::CGLWidget(const QGLFormat& fmt, QWidget *parent, QGLWidget *sharedWidget)
  : QGLWidget(fmt, parent, sharedWidget), 
    _fovy(30.f), _znear(0.1f), _zfar(10.f), 
    _eye(0, 0, 2.5), _center(0, 0, 0), _up(0, 1, 0), 
    rc(NULL),
    toggle_mesh(true), toggle_wireframe(false), toggle_extrema(false), toggle_labels(false), 
    current_slice(0)
{
}

CGLWidget::~CGLWidget()
{
  rc_destroy_ctx(&rc);
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

  case Qt::Key_E:
    toggle_extrema = !toggle_extrema;
    updateGL();
    break;

  case Qt::Key_L:
    toggle_labels = !toggle_labels;
    updateGL();
    break;

  case Qt::Key_Right:
    current_slice = (current_slice + 1) % nPhi;
    updateGL();
    break;

  case Qt::Key_Left:
    current_slice = (current_slice - 1 + nPhi) % nPhi;
    updateGL();
    break;

  case Qt::Key_P: {
    QString filename = QFileDialog::getSaveFileName(this, "save to png", "./", "*.png");
    if (filename.isEmpty()) return;
    QImage img = grabFrameBuffer(false);
    img.save(filename);
  }

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
  glewInit();
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

  { // volren
    rc_create_ctx(&rc);
    rc_set_stepsize(rc, 0.5);
  }

  CHECK_GLERROR(); 
}

void CGLWidget::resizeGL(int w, int h)
{
  _trackball.reshape(w, h); 
  glViewport(0, 0, w, h);

  rc_set_viewport(rc, 0, 0, w, h);

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

  // renderSinglePlane();
  renderMultiplePlanes();

  CHECK_GLERROR(); 
}

void CGLWidget::renderExtremum() 
{
  glPointSize(4.f);
  
  glColor3f(1, 0, 0);
  glBegin(GL_POINTS);
  for (int i=0; i<maximum[current_slice].size(); i++) {
    int k = maximum[current_slice][i];
    glVertex2f(coords[k*2], coords[k*2+1]);
  }
  glEnd();
  
  glColor3f(0, 1, 0);
  glBegin(GL_POINTS);
  for (int i=0; i<minimum[current_slice].size(); i++) {
    int k = minimum[current_slice][i];
    glVertex2f(coords[k*2], coords[k*2+1]);
  }
  glEnd();
}

void CGLWidget::renderLabels()
{
  glPointSize(4.f);

  const std::vector<size_t> &labels = all_labels[current_slice];

  glBegin(GL_POINTS);
  for (int i=0; i<nNodes; i++) {
    int label = labels[i];
    QColor c = label_colors[label];
    glColor3ub(c.red(), c.green(), c.blue());
    glVertex2f(coords[i*2], coords[i*2+1]);
  }
  glEnd();
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

  glScalef(0.4, 0.4, 0.4);
  for (int i=0; i<nPhi; i++) {
    glPushMatrix();
    glRotatef(360.f/nPhi*i, 0, 1, 0);

#if 1
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    glVertexPointer(2, GL_FLOAT, 0, f_vertices.data());
    glColorPointer(3, GL_FLOAT, 0, f_colors.data());
    glDrawArrays(GL_TRIANGLES, 0, f_vertices.size()/2);

    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
#endif

#if 0
    glPointSize(4.f);
    
    glColor4f(1, 0, 0, 0.8);
    glBegin(GL_POINTS);
    for (int j=0; j<maximum[i].size(); j++) {
      int k = maximum[i][j];
      glVertex2f(coords[k*2], coords[k*2+1]);
    }
    glEnd();
    
    glColor4f(0, 1, 0, 0.8);
    glBegin(GL_POINTS);
    for (int j=0; j<minimum[i].size(); j++) {
      int k = minimum[i][j];
      glVertex2f(coords[k*2], coords[k*2+1]);
    }
    glEnd();
#endif
    glPopMatrix();
  }
}

void CGLWidget::renderSinglePlane()
{
  if (f_vertices.size() == 0) return;

  glColor3f(0, 0, 0);

  glScalef(0.4, 0.4, 0.4);
  glTranslatef(-1.7, 0, 0);

  if (toggle_mesh) {
    if (toggle_wireframe) {
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // render wireframe
    }
    else {
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    glVertexPointer(2, GL_FLOAT, 0, f_vertices.data());
    glColorPointer(3, GL_FLOAT, 0, &f_colors[current_slice*nTriangles*3*3] );
    glDrawArrays(GL_TRIANGLES, 0, f_vertices.size()/2);

    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
  }

  if (toggle_labels)
    renderLabels();

  if (toggle_extrema)
    renderExtremum();

  CHECK_GLERROR();
}
  
void CGLWidget::setTriangularMesh(int nNodes_, int nTriangles_, int nPhi_, double *coords_, int *conn_)
{
  nNodes = nNodes_; 
  nTriangles = nTriangles_;
  nPhi = nPhi_;
  coords = coords_;
  conn = conn_;

  f_vertices.clear();

  for (int i=0; i<nTriangles; i++) {
    int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    f_vertices.push_back(coords[i0*2]);
    f_vertices.push_back(coords[i0*2+1]);
    f_vertices.push_back(coords[i1*2]);
    f_vertices.push_back(coords[i1*2+1]);
    f_vertices.push_back(coords[i2*2]);
    f_vertices.push_back(coords[i2*2+1]);
  }

  // init nodeGraph
  nodeGraph.clear();
  nodeGraph.resize(nNodes);
  for (int i=0; i<nTriangles; i++) {
    int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    if (i0 == i1 || i1 == i2 || i2 == i0) 
      fprintf(stderr, "WARN : triangle %d, conn={%d, %d, %d}\n", i, i0, i1, i2);
    if (i0 == i1 && i1 == i2)
      fprintf(stderr, "FATAL: triangle %d, conn={%d, %d, %d}\n", i, i0, i1, i2);
    nodeGraph[i0].insert(i1);
    nodeGraph[i0].insert(i2);
    nodeGraph[i1].insert(i2);
    nodeGraph[i1].insert(i0);
    nodeGraph[i2].insert(i0);
    nodeGraph[i2].insert(i1);
  }
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

void CGLWidget::constructDiscreteGradient(double *dpot) // based on MS theory
{
  // 0-cell: vertices
  // 1-cell: edges
  // 2-cell: cells
  
  typedef struct {
    unsigned char index;
    double value;
    std::vector<int> hCells;
  } msCell;

}

struct XGCData { // either single slice or all slices
  std::vector<std::set<int> > *nodeGraph;
  double *dpot;
  int nNodes, nPhi;
};

static double value(size_t v, void *d)
{
  XGCData *data = (XGCData*)d;
  // fprintf(stderr, "%d, %f\n", v, data->dpot[v]);
  return data->dpot[v];
}

static size_t neighbors(size_t v, size_t *nbrs, void *d)
{
  XGCData *data = (XGCData*)d;
  std::set<int>& nodes = (*data->nodeGraph)[v];
  
  int i = 0;
  for (std::set<int>::iterator it = nodes.begin(); it != nodes.end(); it ++) {
    nbrs[i] = *it;
    i ++;
  }

  return i;
}

static size_t neighbors3D(size_t v, size_t *nbrs, void *d)
{
  XGCData *data = (XGCData*)d;
 
  const size_t nNodes = data->nNodes;
  const size_t nPhi = data->nPhi;

  size_t plane = v / nNodes;
  size_t node = v % nNodes;

  // fprintf(stderr, "nnodes=%d, nphi=%d, plane=%d, node=%d\n", nNodes, nPhi, plane, node);

  std::set<int>& nodes2D = (*data->nodeGraph)[node];
  
  int i = 0;
  for (std::set<int>::iterator it = nodes2D.begin(); it != nodes2D.end(); it ++) {
    size_t nbr2D = *it;
    nbrs[i++] = ((plane - 1 + nPhi) % nPhi) * nNodes + nbr2D;
    nbrs[i++] = plane * nNodes + nbr2D;
    nbrs[i++] = ((plane + 1 + nPhi) % nPhi) * nNodes + nbr2D;
  }

  return i;
}

static double volumePriority(ctNode *node, void *d)
{
  XGCData *data = static_cast<XGCData*>(d);
  ctArc *arc = ctNode_leafArc(node);

  ctNode *lo = arc->lo,
         *hi = arc->hi;

  double val_lo = value(lo->i, d),
         val_hi = value(hi->i, d);

  int count = 0;

  std::queue<size_t> Q;
  std::set<size_t> visited;

  Q.push(lo->i);
  visited.insert(lo->i); 

  while (!Q.empty()) {
    size_t p = Q.front(); 
    Q.pop();
    double val_p = value(p, d);

    std::set<int> &neighbors = (*data->nodeGraph)[p];
    for (std::set<int>::iterator it = neighbors.begin(); it != neighbors.end(); it ++) {
      const int neighbor = *it; 

      if (visited.find(neighbor) == visited.end()) { // not found
        double val_q = value(neighbor, d);
        if (val_q >= val_lo && val_q < val_hi) {
          Q.push(neighbor);
          visited.insert(neighbor);
          count ++;
        }
      }
    }
  }

  fprintf(stderr, "volume=%d\n", count);
  return count;
}

static void printContourTree(ctBranch* b)
{
  fprintf(stderr, "%d, %d, %p\n", b->extremum, b->saddle, b->children.head);
  
  for (ctBranch* c = b->children.head; c != NULL; c = c->nextChild) 
    printContourTree(c);
}

int CGLWidget::flood2D(size_t seed, size_t id, std::vector<size_t> &labels, double min, double max, void *d)
{
  fprintf(stderr, "seed=%lu, id=%lu, min=%.10e, max=%.10e\n", seed, id, min, max);

  XGCData *data = (XGCData*)d;

  int count = 1;
  std::queue<size_t> Q;
  Q.push(seed);
  labels[seed] = id;

  while (!Q.empty()) {
    size_t p = Q.front(); 
    Q.pop();

    size_t nbrs[32]; //TODO
    size_t nnbrs = neighbors(p, nbrs, d);

    for (int i=0; i<nnbrs; i++) {
      const size_t neighbor = nbrs[i];

      if (labels[neighbor] != id) {
        double val = value(neighbor, d);
        if (val >= min && val <= max) {
          // fprintf(stderr, "true for %d, val=%.10e, min=%.10e, max=%.10e\n", neighbor, val, min, max);
          Q.push(neighbor);
          labels[neighbor] = id;
          count ++;
        }
      }
    }
  }

  return count;
}
  
void CGLWidget::extractStreamers(int plane, ctBranch *root, std::map<ctBranch*, size_t>& branchSet, int nStreamers, double percentage, void *d)
{
  std::vector<std::pair<const ctBranch*, double> > branchPersistences;
  for (auto &kv : branchSet) {
    branchPersistences.push_back(std::make_pair<const ctBranch*, double>(kv.first, value(kv.first->extremum, d) - value(kv.first->saddle, d)));
  }
  branchPersistences.push_back(std::make_pair<const ctBranch*, double>(root, value(root->saddle, d) - value(root->extremum, d)));

  std::sort(branchPersistences.begin(), branchPersistences.end(), 
    [](std::pair<const ctBranch*, double> p0, std::pair<const ctBranch*, double> p1) {
      return p0.second < p1.second;
    });

  std::vector<size_t> labels(nNodes, 0);
 
  static int id = 1;
 
  // minimum streamers
  for (int i=0; i<std::min((size_t)nStreamers, branchPersistences.size()); i++) {
    size_t extremum;
    if (branchPersistences[i].first == root)
      extremum = branchPersistences[i].first->saddle; 
    else 
      extremum = branchPersistences[i].first->extremum;

    double val_extremum = value(extremum, d);
    const int myid = 1; // id++; // -i
    // label_colors[myid] = QColor(rand()%256, rand()%256, rand()%256);
    label_colors[myid] = QColor(0, 255, 0);
    int count = flood2D(extremum, myid, labels, val_extremum, percentage*val_extremum, d); // presumably the val is less than 0
    fprintf(stderr, "count=%d\n", count);
  }

  // maximum streamers
  for (int i=branchPersistences.size()-1; i>branchPersistences.size()-nStreamers-1; i--) {
    size_t extremum = branchPersistences[i].first->extremum; // maximum
    double val_extremum = value(extremum, d);
    const int myid = 2; // id++;
    // label_colors[myid] = QColor(rand()%256, rand()%256, rand()%256);
    label_colors[myid] = QColor(255, 0, 0);
    int count = flood2D(extremum, myid, labels, val_extremum*percentage, val_extremum, d);
    fprintf(stderr, "count=%d\n", count);
  }

  all_labels[plane] = labels;
}

void CGLWidget::addExtremumFromBranchDecomposition(int plane, ctBranch *root, ctBranch *b, void *d)
{
  double val_extremum = value(b->extremum, d), 
         val_saddle = value(b->saddle, d);

  if (b == root) {
    if (val_extremum >= val_saddle) {
      maximum[plane].push_back(b->extremum);
      minimum[plane].push_back(b->saddle);
    } else {
      maximum[plane].push_back(b->saddle);
      minimum[plane].push_back(b->extremum);
    }
  } else {
    if (val_extremum >= val_saddle) 
      maximum[plane].push_back(b->extremum);
    else 
      minimum[plane].push_back(b->extremum);
  }

  for (ctBranch *c = b->children.head; c != NULL; c = c->nextChild) {
    addExtremumFromBranchDecomposition(plane, root, c, d);
  }
}

void CGLWidget::simplifyBranchDecompositionByThreshold(ctBranch *b, double threshold, void *d)
{
  // XGCData *data = (XGCData*)d;

RESTART: 
  for (ctBranch *c=b->children.head; c!=NULL; c=c->nextChild)
    if (fabs(value(c->extremum, d) - value(c->saddle, d)) < threshold) {
      ctBranchList_remove(&b->children, c);
      goto RESTART;
    } else {
      simplifyBranchDecompositionByThreshold(c, threshold, d);
    }
}

std::set<ctBranch*> CGLWidget::simplifyBranchDecompositionByNumbers(ctBranch* rootBranch, std::map<ctBranch*, size_t> &branchSet, int nLimit, void *d) // no more than a number
{
  float minPriority; // default is persistence
  ctBranch *minPriorityBranch = NULL; 

  std::list<ctBranch*> branchList;
  for (auto &kv : branchSet) {
    branchList.push_back(kv.first);
  }

  std::set<size_t> removed_labels;
  std::set<ctBranch*> removed_branches;

  branchList.sort(
      [&d](ctBranch* b0, ctBranch* b1) {
        return fabs(value(b0->extremum, d) - value(b0->saddle, d)) < fabs(value(b1->extremum, d) - value(b1->saddle, d));
        // return value(b0->extremum, d) - value(b0->saddle, d) < value(b1->extremum, d) - value(b1->saddle, d);
      });

  while (minPriorityBranch != rootBranch && branchSet.size() > nLimit) {
    minPriority = 1e38f; 
    minPriorityBranch = rootBranch; 

    std::list<ctBranch*>::iterator it_erase = branchList.end();
    for (std::list<ctBranch*>::iterator it = branchList.begin(); it != branchList.end(); it ++) {
      ctBranch *b = *it;
      if (b->children.head == NULL) {
        minPriority = fabs(value(b->extremum, d) - value(b->saddle, d));
        minPriorityBranch = b;
        it_erase = it;
        break;
      }
    }

    if (it_erase != branchList.end())
      branchList.erase(it_erase);

    // qDebug() << "priority =" << minPriority << "branch =" << minPriorityBranch; 

    if (minPriorityBranch == rootBranch) break; 
    ctBranch *parentBranch = minPriorityBranch->parent;

    removed_branches.insert(minPriorityBranch);
    removed_labels.insert(branchSet[minPriorityBranch]);

    ctBranchList_remove(&parentBranch->children, minPriorityBranch); 
    branchSet.erase(minPriorityBranch); 
  }

  return removed_branches;
}

void CGLWidget::buildContourTree(int plane, double *dpot_)
{
  fprintf(stderr, "building contour tree for plane %d, nNodes=%d\n", plane, nNodes);

  std::vector<size_t> totalOrder; 
  for (int i=0; i<nNodes; i++) 
    totalOrder.push_back(i);

  XGCData data;
  data.nodeGraph = &nodeGraph;
  data.dpot = dpot_ + plane * nNodes;
  data.nNodes = nNodes;
  data.nPhi = nPhi;

  std::stable_sort(totalOrder.begin(), totalOrder.end(),
      [&data](size_t v0, size_t v1) {
        return data.dpot[v0] < data.dpot[v1];
        // if (fabs(data.dpot[v0] - data.dpot[v1]) < 1e-5) return v0 < v1; 
        // else return data.dpot[v0] < data.dpot[v1];
      });

  ctContext *ctx = ct_init(
      nNodes, 
      &totalOrder.front(),  
      &value,
      &neighbors, 
      &data);

  // ct_priorityFunc(ctx, volumePriority);

  ct_sweepAndMerge(ctx);
  ctBranch *root = ct_decompose(ctx);
  ctBranch **branchMap = ct_branchMap(ctx);

  // build branchSet
  std::map<ctBranch*, size_t> branchSet;
  size_t nBranches = 0;
  for (int i=0; i<nNodes; i++) {
    if (branchSet.find(branchMap[i]) == branchSet.end()) {
      size_t branchId = nBranches ++;
      branchSet[branchMap[i]] = branchId;
    }
  }

  for (auto &kv : branchSet) {
    label_colors[kv.second] = QColor(rand()%256, rand()%256, rand()%256);
  }

  // simplifyBranchDecompositionByThreshold(root, 3e-10, &data);
  std::set<ctBranch*> removed_branches = simplifyBranchDecompositionByNumbers(root, branchSet, nLimit, &data);
  // addExtremumFromBranchDecomposition(plane, root, root, &data);
  // extractStreamers(plane, root, branchSet, 80, 0.1, &data);

  double diff = 0, total = 0;
  size_t count_simplified = 0;

  std::vector<size_t> labels(nNodes); 
  for (size_t i=0; i<nNodes; i++) {
    ctBranch *b = branchMap[i];
    ctBranch *bp = NULL;
    while (removed_branches.find(b) != removed_branches.end()) {
      bp = b;
      b = b->parent;
    }
    labels[i] = branchSet[b];
      
    double &src = data.dpot[i];
    total += src*src;
    if (bp != NULL) {
      double dst = value(bp->saddle, &data);
      // fprintf(stderr, "src=%f, dst=%f\n", src, dst);
      diff += src*src - dst*dst;
      src = dst;
      count_simplified ++;
    }
  }
  all_labels[plane] = labels;

  fprintf(stderr, "total=%f, diff=%f, loss=%f, count_simplified=%zu\n", 
      total, diff, diff/total, count_simplified);

#if 0 // 2D topology segmentation
  std::vector<size_t> labels(nNodes, 0);
  buildSegmentation(root, labels, &data);
  // for (int i=0; i<nNodes; i++) 
  //   fprintf(stderr, "%d, %d\n", i, labels[i]);

  all_labels[plane] = labels;

  // printContourTree(root);
#else // threshold-based streamer extraction
  // extractStreamersFromExtremum(plane, dpot_, 0.1);
#endif

  ct_cleanup(ctx);
}

void CGLWidget::buildContourTree3D(double *dpot_)
{
  fprintf(stderr, "building contour tree over 3D...\n");

  std::vector<size_t> totalOrder; 
  for (int i=0; i<nNodes*nPhi; i++) 
    totalOrder.push_back(i);

  XGCData data;
  data.nodeGraph = &nodeGraph;
  data.dpot = dpot_;
  data.nNodes = nNodes;
  data.nPhi = nPhi;
  
  std::stable_sort(totalOrder.begin(), totalOrder.end(),
      [&data](size_t v0, size_t v1) {
        return data.dpot[v0] < data.dpot[v1];
        // if (fabs(data.dpot[v0] - data.dpot[v1]) < 1e-5) return v0 < v1; 
        // else return data.dpot[v0] < data.dpot[v1];
      });
  
  ctContext *ctx = ct_init(
      nNodes * nPhi, 
      &totalOrder.front(),  
      &value,
      &neighbors3D, 
      &data);
  
  ct_sweepAndMerge(ctx);
  ctBranch *root = ct_decompose(ctx);
  ctBranch **branchMap = ct_branchMap(ctx);
  ctArc **arcMap = ct_arcMap(ctx);

#if 0
  std::map<ctArc*, size_t> arcSet;
  size_t nArcs = 0;
  for (int i=0; i<nNodes*nPhi; i++) {
    if (arcSet.find(arcMap[i]) == arcSet.end()) {
      size_t arcId = nArcs ++;
      // arcSet.insert(arcMap[i], arcId);
      arcSet[arcMap[i]] = arcId;
      fprintf(stderr, "arcId=%d, %p, hi=%p, lo=%p\n", arcId, arcMap[i], arcMap[i]->hi, arcMap[i]->lo);
      label_colors[arcId] = QColor(rand()%256, rand()%256, rand()%256);
    }
  }
 
  for (int i=0; i<nPhi; i++) {
    std::vector<size_t> labels(nNodes, 0);
    for (int j=0; j<nNodes; j++) {
      size_t arcId = arcSet[arcMap[i*nNodes + j]];
      labels[j] = arcId;
    }
    all_labels[i] = labels;
  }
#else
  std::map<ctBranch*, size_t> branchSet;
  size_t nBranches = 0;
  for (int i=0; i<nNodes*nPhi; i++) {
    if (branchSet.find(branchMap[i]) == branchSet.end()) {
      size_t branchId = nBranches ++;
      branchSet[branchMap[i]] = branchId;
    }
  }
  for (auto &kv : branchSet) 
    label_colors[kv.second] = QColor(rand()%256, rand()%256, rand()%256);
 
  std::set<ctBranch*> removed_branches = simplifyBranchDecompositionByNumbers(root, branchSet, nLimit, &data);
 
  for (int i=0; i<nPhi; i++) {
    std::vector<size_t> labels(nNodes, 0);
    for (int j=0; j<nNodes; j++) {
      ctBranch *b = branchMap[i*nNodes + j];
      while (removed_branches.find(b) != removed_branches.end())
        b = b->parent;
      size_t branchId = branchSet[b];
      labels[j] = branchId;
    }
    all_labels[i] = labels;
  }
#endif

#if 0
  // simplifyBranchDecompositionByThreshold(root, nLimit, &data);
  
  std::vector<size_t> labels(nNodes*nPhi, 0);
  buildSegmentation3D(root, labels, &data);
  // TODO

  for (int i=0; i<nPhi; i++) {
    std::vector<size_t> l(nNodes);
    std::copy(labels.begin() + i*nNodes, labels.begin() + (i+1)*nNodes - 1, l.begin());
    all_labels[i] = l;
  }
  
  ct_cleanup(ctx); 
  // TODO
#endif

  fprintf(stderr, "built contour tree over 3D.\n");
}

void CGLWidget::buildSegmentation3D(ctBranch *b, std::vector<size_t> &labels, void *d)
{
  static int id = 1; // FIXME
  const int currentId = ++ id;

  label_colors[currentId] = QColor(rand()%256, rand()%256, rand()%256);

  XGCData *data = (XGCData*)d;

  const double val_extremum = value(b->extremum, d), 
               val_saddle = value(b->saddle, d);
  const double val_hi = std::max(val_extremum, val_saddle),
               val_lo = std::min(val_extremum, val_saddle);

  int count = 1;

  std::queue<size_t> Q;
  Q.push(b->extremum);
  labels[b->extremum] = currentId; 

  while (!Q.empty()) {
    size_t p = Q.front();
    Q.pop();

    size_t nbrs[128]; // TODO
    size_t nnbrs = neighbors3D(p, nbrs, d);

    for (int i=0; i<nnbrs; i++) {
      const size_t neighbor = nbrs[i];

      if (labels[neighbor] != currentId) {
        double val = value(neighbor, d);
        bool qualify = true;
        if (val_extremum > val_saddle) {
          if (val > val_extremum || val <= val_saddle) qualify = false;
        } else {
          if (val < val_extremum || val >= val_saddle) qualify = false;
        }
        
        if (qualify) {
          Q.push(neighbor);
          labels[neighbor] = currentId;
          count ++;
        }
      }
    }
  }

  // fprintf(stderr, "#count=%d\n", count);

  for (ctBranch *c = b->children.head; c != NULL; c = c->nextChild) 
    buildSegmentation3D(c, labels, d);
}

void CGLWidget::buildSegmentation(ctBranch *b, std::vector<size_t> &labels, void *d)
{
  static int id = 1; // FIXME
  const int currentId = ++ id;

  label_colors[currentId] = QColor(rand()%256, rand()%256, rand()%256);

  XGCData *data = (XGCData*)d;

  const double val_extremum = value(b->extremum, d), 
               val_saddle = value(b->saddle, d);
  const double val_hi = std::max(val_extremum, val_saddle),
               val_lo = std::min(val_extremum, val_saddle);

  int count = 1;

  std::queue<size_t> Q;
  Q.push(b->extremum);
  labels[b->extremum] = currentId; 

  while (!Q.empty()) {
    size_t p = Q.front();
    Q.pop();

    std::set<int> &neighbors = nodeGraph[p];
    for (std::set<int>::iterator it = neighbors.begin(); it != neighbors.end(); it ++) {
      const int neighbor = *it;

      if (labels[neighbor] != currentId) {
        double val = value(neighbor, d);
        bool qualify = true;
        if (val_extremum > val_saddle) {
          if (val > val_extremum || val <= val_saddle) qualify = false;
        } else {
          if (val < val_extremum || val >= val_saddle) qualify = false;
        }
        
        if (qualify) {
          Q.push(neighbor);
          labels[neighbor] = currentId;
          count ++;
        }
      }
    }
  }

  // fprintf(stderr, "#count=%d\n", count);

  for (ctBranch *c = b->children.head; c != NULL; c = c->nextChild) 
    buildSegmentation(c, labels, d);
}

void CGLWidget::extractExtremum(int plane, double *dpot_)
{
  double *dpot = dpot_ + plane * nNodes;

  for (int i=0; i<nNodes; i++) {
    std::set<int> &neighbors = nodeGraph[i];
    bool local_max = true, local_min = true;

    for (std::set<int>::iterator it = neighbors.begin(); it != neighbors.end(); it ++) {
      const int j = *it;
      if (dpot[i] >= dpot[j]) local_min = false;
      if (dpot[i] <= dpot[j]) local_max = false;
    }

    if (local_max)
      maximum[plane].push_back(i); 
    else if (local_min)
      minimum[plane].push_back(i);
  }
}

void CGLWidget::extractStreamersFromExtremum(int plane, double *dpot, double percentage)
{
  XGCData data;
  data.nodeGraph = &nodeGraph;
  data.dpot = dpot + plane * nNodes;
  data.nNodes = nNodes;
  data.nPhi = nPhi;

  struct Streamer {
    size_t id;
    size_t seed;
    double min, max; 
  };

  int count = 0;
  std::vector<Streamer> streamers;
  for (int i=0; i<maximum[plane].size(); i++) {
    Streamer s; 
    s.id = 1 + count ++;
    s.seed = maximum[plane][i];
    s.max = value(s.seed, &data);
    s.min = 0.1 * s.max;
    streamers.push_back(s);
  }

  std::vector<size_t> labels(nNodes, 0);
  for (int i=0; i<streamers.size(); i++) {
    int count = flood2D(streamers[i].seed, streamers[i].id, labels, streamers[i].min, streamers[i].max, &data);
    label_colors[streamers[i].id] = QColor(rand()%256, rand()%256, rand()%256);
    fprintf(stderr, "count=%d\n", count);
  }

  all_labels[plane] = labels;
}

void CGLWidget::updateDpot(double *dpot)
{
  // const float min = -100, max = 100;
  const float min = -50, max = 50;
  // const float min = -1e-9, max = 1e-9;

  f_colors.clear();

  for (int plane = 0; plane < nPhi; plane ++) {
    for (int i=0; i<nTriangles; i++) {
      int v[3] = {conn[i*3], conn[i*3+1], conn[i*3+2]};

      for (int j=0; j<3; j++) {
        float val = clamp_normalize(min, max, (float)dpot[plane*nNodes + v[j]]);
        f_colors.push_back(val);
        f_colors.push_back(1-val);
        f_colors.push_back(0);
      }
    }
  }
}

void CGLWidget::setData(double *dpot)
{
  // updateDpot(dpot);
  // buildContourTree3D(dpot);

#if 0
  for (int i=0; i<nPhi; i++) {
  // for (int i=0; i<nPhi; i++) {
    buildContourTree(i, dpot); 
    extractExtremum(i, dpot); // dpot + i*nNodes);
  }
#endif

  updateDpot(dpot);
}
