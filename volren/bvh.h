#ifndef _BVH_H
#define _BVH_H

#include <algorithm>
#include <vector>
#include <cstdio>
#include <cfloat>
#include <cstdlib>
#include "aabb.h"
#include "bvh.cuh"

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

std::vector<BVHNodeD> convertBVH(QuadNode* r, const std::vector<int> &conn, const std::vector<double> &coords);

void deleteBVH(QuadNode *q);

std::vector<BVHNodeD> buildBVHGPU(int nNodes, int nTriangles, const std::vector<double> &coords, const std::vector<int> &conn);

#endif
