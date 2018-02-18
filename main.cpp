#include <mpi.h>
#include <adios.h>
#include <adios_read.h>
#include <cassert>
#include <cfloat>
#include <chrono>
#include <queue>
#include <stack>
#include <QApplication>
#include "widget.h"
#include "volren.cuh"

extern "C"
{
extern void adios_read_bp_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran);
}

bool readValueInt(ADIOS_FILE *fp, const char *nm, int *val)
{
  ADIOS_VARINFO *avi = adios_inq_var(fp, nm);
  if (avi == NULL || avi->type != adios_integer) return false;
  *val = *((int*)avi->value);
  return true;
}


template <typename T> int readScalars(ADIOS_FILE *fp, const char *var, T **arr)
{
  ADIOS_VARINFO *avi = adios_inq_var(fp, var);
  if (avi == NULL) return false;

  adios_inq_var_stat(fp, avi, 0, 0);
  adios_inq_var_blockinfo(fp, avi);
  adios_inq_var_meshinfo(fp, avi);
    
  int nt = 1;
  uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {0, 0, 0, 0};
  
  for (int i = 0; i < avi->ndim; i++) {
    st[i] = 0;
    sz[i] = avi->dims[i];
    nt = nt * sz[i];
  }
  // fprintf(stderr, "%d, %d, %d, %d\n", sz[0], sz[1], sz[2], sz[3]);
  
  ADIOS_SELECTION *sel = adios_selection_boundingbox(avi->ndim, st, sz);
  assert(sel->type == ADIOS_SELECTION_BOUNDINGBOX);

  *arr = (T*)malloc(sizeof(double)*nt);

  adios_schedule_read_byid(fp, sel, avi->varid, 0, 1, *arr);
  int retval = adios_perform_reads(fp, 1);
  
  adios_selection_delete(sel);
  return avi->ndim;
}

bool readTriangularMesh(ADIOS_FILE *fp, int &nNodes, int &nTriangles, double **coords, int **conn, int **nextNode) 
{
  // read number of nodes and triangles
  readValueInt(fp, "n_n", &nNodes);
  readValueInt(fp, "n_t", &nTriangles);
  
  // read coordinates
  readScalars<double>(fp, "/coordinates/values", coords);

  // read triangles
  readScalars<int>(fp, "/cell_set[0]/node_connect_list", conn);

  // read nextNode
  readScalars<int>(fp, "nextnode", nextNode);
  
  fprintf(stderr, "mesh loaded: nNodes=%d, nTriangles=%d\n", 
      nNodes, nTriangles);

  return true;
}

void writeVTK(const int nNodes, const int nTriangles, double *coords, int *conn) 
{

}

template <typename T>
inline T min3(T x, T y, T z) {
  return std::min(std::min(x, y), z);
}

template <typename T>
inline T max3(T x, T y, T z) {
  return std::max(std::max(x, y), z);
}

struct AABB {
  int id = 0;
  double A[2] = {FLT_MAX, FLT_MAX}, B[2] = {FLT_MIN, FLT_MIN};
  double C[2] = {0, 0}; // centroid

  bool contains(const double *X) {
    return X[0] >= A[0] && X[0] < B[0] && X[1] >= A[1] && X[1] < B[1];
  }

  void updateCentroid() {
    C[0] = (A[0] + B[0]) / 2;
    C[1] = (A[1] + B[1]) / 2;
  }

  void print() const {
    fprintf(stderr, "A={%f, %f}, B={%f, %f}, centroid={%f, %f}\n", 
        A[0], A[1], B[0], B[1], C[0], C[1]);
  }
};

struct QuadNode {
  QuadNode *parent = NULL;
  QuadNode *children[4] = {NULL};
  AABB aabb;
  std::vector<AABB*> elements; 

  bool isLeaf() const {
    return elements.size() > 0;
    // return children[0] == NULL && children[1] == NULL 
    //   && children[2] == NULL && children[3] == NULL;
  }

  void updateBounds() {
    for (int i=0; i<elements.size(); i++) {
      AABB *aabb1 = elements[i];
      aabb.A[0] = std::min(aabb.A[0], aabb1->A[0]);
      aabb.A[1] = std::min(aabb.A[1], aabb1->A[1]);
      aabb.B[0] = std::max(aabb.B[0], aabb1->B[0]);
      aabb.B[1] = std::max(aabb.B[1], aabb1->B[1]);
    }
    aabb.updateCentroid();
  }

  void print() const {
    fprintf(stderr, "parent=%p, A={%f, %f}, B={%f, %f}, centroid={%f, %f}, children={%p, %p, %p, %p}, element=%d\n", 
        parent, aabb.A[0], aabb.A[1], aabb.B[0], aabb.B[1], aabb.C[0], aabb.C[1], 
        children[0], children[1], children[2], children[3], 
        elements.empty() ? -1 : elements[0]->id);
  }
};

void traverseQuadNode(QuadNode *q)
{
  q->print();
  for (int j=0; j<4; j++) 
    if (q->children[j] != NULL) 
      traverseQuadNode(q->children[j]);
}

void subdivideQuadNode(QuadNode *q) 
{
  if (q->elements.size() <= 1) {
    q->updateBounds();
    return;
  }

  // fprintf(stderr, "subdividing %p, parent=%p, #elements=%zu\n", q, q->parent, q->elements.size());
  for (int j=0; j<4; j++) {
    q->children[j] = new QuadNode;
    q->children[j]->parent = q;
  }

  // left-bottom
  q->children[0]->aabb.A[0] = q->aabb.A[0];
  q->children[0]->aabb.A[1] = q->aabb.A[1];
  q->children[0]->aabb.B[0] = q->aabb.C[0];
  q->children[0]->aabb.B[1] = q->aabb.C[1];
  q->children[0]->aabb.updateCentroid();
 
  // right-bottom
  q->children[1]->aabb.A[0] = q->aabb.C[0];
  q->children[1]->aabb.A[1] = q->aabb.A[1];
  q->children[1]->aabb.B[0] = q->aabb.B[0];
  q->children[1]->aabb.B[1] = q->aabb.C[1];
  q->children[1]->aabb.updateCentroid();
  
  // right-top
  q->children[2]->aabb.A[0] = q->aabb.C[0];
  q->children[2]->aabb.A[1] = q->aabb.C[1];
  q->children[2]->aabb.B[0] = q->aabb.B[0];
  q->children[2]->aabb.B[1] = q->aabb.B[1];
  q->children[2]->aabb.updateCentroid();
  
  // left-top
  q->children[3]->aabb.A[0] = q->aabb.A[0];
  q->children[3]->aabb.A[1] = q->aabb.C[1];
  q->children[3]->aabb.B[0] = q->aabb.C[0];
  q->children[3]->aabb.B[1] = q->aabb.B[1];
  q->children[3]->aabb.updateCentroid();

  for (int i=0; i<q->elements.size(); i++) {
    for (int j=0; j<4; j++) {
      if (q->children[j]->aabb.contains(q->elements[i]->C)) {
        q->children[j]->elements.push_back(q->elements[i]);
        break;
      }
    }
  }

  if (q->parent != NULL) 
    q->updateBounds();
  q->elements.clear();

  for (int j=0; j<4; j++) {
    if (q->children[j]->elements.empty()) {
      delete q->children[j];
      q->children[j] = NULL;
    } else {
      subdivideQuadNode(q->children[j]);
    }
  }
}

bool insideTriangle(const double *p, const double *p1, const double *p2, const double *p3)
{
  double alpha = ((p2[1] - p3[1])*(p[0] - p3[0]) + (p3[0] - p2[0])*(p[1] - p3[1])) /
          ((p2[1] - p3[1])*(p1[0] - p3[0]) + (p3[0] - p2[0])*(p1[1] - p3[1]));
  double beta = ((p3[1] - p1[1])*(p[0] - p3[0]) + (p1[0] - p3[0])*(p[1] - p3[1])) /
         ((p2[1] - p3[1])*(p1[0] - p3[0]) + (p3[0] - p2[0])*(p1[1] - p3[1]));
  double gamma = 1.0 - alpha - beta;
  // fprintf(stderr, "barycentric: %f, %f, %f\n", alpha, beta, gamma);
  return alpha >= 0 && beta >= 0 && gamma >= 0;
}

int locatePointBruteForce(const double *X, QuadNode *q, int nNodes, int nTriangles, const double *coords, const int *conn)
{
  int rtn = -1;
  for (int id=0; id<nTriangles; id++) {
    const int i0 = conn[id*3], i1 = conn[id*3+1], i2 = conn[id*3+2];
    bool succ = insideTriangle(X, &coords[i0*2], &coords[i1*2], &coords[i2*2]);
    if (succ) {
      rtn = id;
      break;
    }
  }
  return rtn;
}

int locatePointNonRecursive(const double *X, QuadNode *q, int nNodes, int nTriangles, const double *coords, const int *conn)
{
  std::stack<QuadNode*> S;
  S.push(q);

  while (!S.empty()) {
    QuadNode *q = S.top();
    S.pop();

    // fprintf(stderr, "checking %p\n", q);
    
    if (q->isLeaf()) {
      const int id = q->elements[0]->id;
      const int i0 = conn[id*3], i1 = conn[id*3+1], i2 = conn[id*3+2];
      int succ = insideTriangle(X, &coords[i0*2], &coords[i1*2], &coords[i2*2]);
      if (succ) return id;
    } else if (q->aabb.contains(X)) {
      for (int j=0; j<4; j++)
        if (q->children[j] != NULL)
          S.push(q->children[j]);
    }
  }

  return -1;
}

int locatePointRecursive(const double *X, QuadNode *q, int nNodes, int nTriangles, const double *coords, const int *conn)
{
  if (q->aabb.contains(X)) {
    if (q->isLeaf()) {
      const int id = q->elements[0]->id;
      const int i0 = conn[id*3], i1 = conn[id*3+1], i2 = conn[id*3+2];
      int succ = insideTriangle(X, &coords[i0*2], &coords[i1*2], &coords[i2*2]);
      if (succ) {
        // fprintf(stderr, "leaf node %d contains! triangle check=%d\n", id, result);
        return id;
      }
    } else {
      for (int j=0; j<4; j++) {
        if (q->children[j] != NULL) {
          int result = locatePointRecursive(X, q->children[j], nNodes, nTriangles, coords, conn);
          if (result >= 0) return result;
        }
      }
    }
  }
  return -1;
}

int convertBVH(QuadNode* r, QuadNodeD** rd_, const int *conn, const double *coords) {
  std::map<QuadNode*, int> nodeMap;
  std::map<int, QuadNode*> nodeReverseMap;
  std::stack<QuadNode*> S;
  S.push(r);

  int quadNodeCount = 0;
  int maxStackSize = 0;
  while (!S.empty()) {
    maxStackSize = std::max(maxStackSize, static_cast<int>(S.size()));

    QuadNode *q = S.top();
    S.pop();

    int nodeId = quadNodeCount ++;
    nodeMap[q] = nodeId;
    nodeReverseMap[nodeId] = q;

    for (int j=0; j<4; j++) 
      if (q->children[j] != NULL)
        S.push(q->children[j]);
  }

  *rd_ = (QuadNodeD*)malloc(sizeof(QuadNodeD)*quadNodeCount);
  QuadNodeD *rd = *rd_;

  for (int i=0; i<quadNodeCount; i++) {
    QuadNode q = *nodeReverseMap[i];
    QuadNodeD &d = rd[i];

    // parent
    if (q.parent == NULL) d.parentId = -1; // root
    else d.parentId = nodeMap[q.parent];

    // children
    for (int j=0; j<4; j++)
      if (q.children[j] == NULL) d.childrenIds[j] = -1;
      else d.childrenIds[j] = nodeMap[q.children[j]];

    // bounds
    d.Ax = q.aabb.A[0];
    d.Ay = q.aabb.A[1];
    d.Bx = q.aabb.B[0];
    d.By = q.aabb.B[1];
    // fprintf(stderr, "%f, %f, %f, %f\n", d.Ax, d.Ay, d.Bx, d.By);

    // triangle
    if (q.isLeaf()) {
      const int id = q.elements[0]->id;
      d.triangleId = id;
      d.i0 = conn[id*3];
      d.i1 = conn[id*3+1];
      d.i2 = conn[id*3+2];
      d.x0 = coords[d.i0*2];
      d.y0 = coords[d.i0*2+1];
      d.x1 = coords[d.i1*2];
      d.y1 = coords[d.i1*2+1];
      d.x2 = coords[d.i2*2];
      d.y2 = coords[d.i2*2+1];
    } else 
      d.triangleId = -1;
  }

  fprintf(stderr, "quadNodeCount=%d, maxStackSize=%d\n", quadNodeCount, maxStackSize);
  return quadNodeCount;
}

void createBVH(int nNodes, int nTriangles, const double *coords, const int *conn)
{
  QuadNode *root = new QuadNode;

  // global bounds
  AABB &aabb = root->aabb;
  for (int i=0; i<nNodes; i++) {
    aabb.A[0] = std::min(aabb.A[0], coords[2*i]);
    aabb.A[1] = std::min(aabb.A[1], coords[2*i+1]);
    aabb.B[0] = std::max(aabb.B[0], coords[2*i]);
    aabb.B[1] = std::max(aabb.B[1], coords[2*i+1]);
  }
  aabb.updateCentroid();
  // aabb.print();

  std::vector<AABB> triangles(nTriangles);
  for (int i=0; i<nTriangles; i++) {
    const int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    double x0 = coords[i0*2], x1 = coords[i1*2], x2 = coords[i2*2], 
           y0 = coords[i0*2+1], y1 = coords[i1*2+1], y2 = coords[i2*2+1];

    triangles[i].C[0] = (x0 + x1 + x2) / 3;
    triangles[i].C[1] = (y0 + y1 + y2) / 3;
    triangles[i].A[0] = min3(x0, x1, x2); 
    triangles[i].A[1] = min3(y0, y1, y2); 
    triangles[i].B[0] = max3(x0, x1, x2); 
    triangles[i].B[1] = max3(y0, y1, y2);
    triangles[i].id = i;

    root->elements.push_back(&triangles[i]);
  }

  subdivideQuadNode(root);
  // traverseQuadNode(root);

  QuadNodeD *rd;
  convertBVH(root, &rd, conn, coords);

#if 1
  fprintf(stderr, "BVH built.\n");
  typedef std::chrono::high_resolution_clock clock;

  const double X[2] = {2.0, 1.03};
  auto t0 = clock::now();
  int r0 = locatePointBruteForce(X, root, nNodes, nTriangles, coords, conn);
  auto t1 = clock::now();
  int r1 = locatePointRecursive(X, root, nNodes, nTriangles, coords, conn);
  auto t2 = clock::now();
  int r2 = locatePointNonRecursive(X, root, nNodes, nTriangles, coords, conn);
  auto t3 = clock::now();
  float alpha, beta, gamma;
  int r3 = QuadNodeD_locatePoint(rd, X[0], X[1], alpha, beta, gamma);
  auto t4 = clock::now();
  int r4 = QuadNodeD_locatePoint_recursive(rd, rd, X[0], X[1], alpha, beta, gamma);
  auto t5 = clock::now();
  fprintf(stderr, "r0=%d, r1=%d, r2=%d, r3=%d, r4=%d\n", r0, r1, r2, r3, r4);

  float tt0 = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
  float tt1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
  float tt2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count();
  float tt3 = std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count();
  float tt4 = std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count();
  fprintf(stderr, "tt0=%f, tt1=%f, tt2=%f, tt3=%f, tt4=%f\n", tt0, tt1, tt2, tt3, tt4);
#endif
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  ADIOS_FILE *meshFP = adios_read_open_file("xgc.mesh.bp", ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  ADIOS_FILE *varFP = adios_read_open_file("xgc.3d.bp", ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
    
  adios_read_bp_reset_dimension_order(meshFP, 0);
  adios_read_bp_reset_dimension_order(varFP, 0);

  int nNodes, nTriangles, nPhi = 1;
  readValueInt(varFP, "nphi", &nPhi);

  fprintf(stderr, "reading mesh...\n");
  double *coords; 
  int *conn, *nextNode;
  readTriangularMesh(meshFP, nNodes, nTriangles, &coords, &conn, &nextNode);
  // fprintf(stderr, "nNodes=%d, nTriangles=%d, nPhi=%d\n", 
  //     nNodes, nTriangles, nPhi);

  createBVH(nNodes, nTriangles, coords, conn);

  // twisting mesh using nextNode
#if 0
  for (int i=0; i<nNodes*3; i++) {
    // conn[i] = nextNode[nextNode[nextNode[nextNode[nextNode[nextNode[nextNode[nextNode[nextNode[conn[i]]]]]]]]]];
    conn[i] = nextNode[conn[i]];
  }
#endif

#if 0
  const double p[2] = {10, 10}, p0[2] = {0, 0}, p1[2] = {10, 0}, p2[2] = {0, 10};
  fprintf(stderr, "awegawegawegaweg, %d\n", 
      insideTriangle(p, p0, p1, p2));
#endif

  fprintf(stderr, "reading data...\n");
  double *dpot;
  readScalars<double>(varFP, "dpot", &dpot);

  fprintf(stderr, "starting GUI...\n");
  QApplication app(argc, argv);
  QGLFormat fmt = QGLFormat::defaultFormat();
  fmt.setSampleBuffers(true);
  fmt.setSamples(16); 
  QGLFormat::setDefaultFormat(fmt); 
 
  const size_t nLimit = atoi(argv[1]);
  fprintf(stderr, "nLimit=%zu\n", nLimit);

  CGLWidget *widget = new CGLWidget;
  widget->setTriangularMesh(nNodes, nTriangles, nPhi, coords, conn);
  widget->setNLimit(nLimit);
  widget->setData(dpot);
  widget->show();
  app.exec();

  free(coords);
  free(conn);

  MPI_Finalize();
}
