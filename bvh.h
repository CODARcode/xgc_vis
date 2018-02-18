#ifndef _BVH_H
#define _BVH_H

#include <algorithm>
#include <vector>
#include <cfloat>
#include <cstdlib>

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

struct QuadNodeD {
  // tree
  int parentId;
  int childrenIds[4];

  // bounds
  float Ax, Ay, Bx, By;

  // triangle
  int triangleId; // -1 if the node if not leaf
  int i0, i1, i2;
  float x0, y0, x1, y1, x2, y2;
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

void traverseQuadNode(QuadNode *q);

void subdivideQuadNode(QuadNode *q);

bool insideTriangle(const double *p, const double *p1, const double *p2, const double *p3);

int locatePointBruteForce(const double *X, QuadNode *q, int nNodes, int nTriangles, const double *coords, const int *conn);

int locatePointNonRecursive(const double *X, QuadNode *q, int nNodes, int nTriangles, const double *coords, const int *conn);

int locatePointRecursive(const double *X, QuadNode *q, int nNodes, int nTriangles, const double *coords, const int *conn);

int convertBVH(QuadNode* r, QuadNodeD** rd_, const int *conn, const double *coords);

void deleteBVH(QuadNode *q);

void buildBVH(int nNodes, int nTriangles, const double *coords, const int *conn);

#endif
