#include "bvh.h"
#include <cassert>
#include <iostream>
#include <chrono>
#include <queue>
#include <stack>
#include <map>

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

std::vector<BVHNodeD> convertBVH(QuadNode* r, const int *conn, const double *coords) {
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

  std::vector<BVHNodeD> rd(quadNodeCount);

  for (int i=0; i<quadNodeCount; i++) {
    QuadNode q = *nodeReverseMap[i];
    BVHNodeD &d = rd[i];

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
  return rd; 
  // return quadNodeCount;
}

void deleteBVH(QuadNode *q) {
  for (int j=0; j<4; j++) {
    if (q->children[j] != NULL) 
      delete q->children[j];
  }
  delete q;
}

std::vector<BVHNodeD> buildBVHGPU(int nNodes, int nTriangles, const double *coords, const int *conn)
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
  aabb.print();

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

  std::vector<BVHNodeD> rd = convertBVH(root, conn, coords);

#if 1
  fprintf(stderr, "BVH built.\n");
  typedef std::chrono::high_resolution_clock clock;

  const double X[2] = {2.3, -0.4};
  auto t0 = clock::now();
  int r0 = locatePointBruteForce(X, root, nNodes, nTriangles, coords, conn);
  auto t1 = clock::now();
  int r1 = locatePointRecursive(X, root, nNodes, nTriangles, coords, conn);
  auto t2 = clock::now();
  int r2 = locatePointNonRecursive(X, root, nNodes, nTriangles, coords, conn);
  auto t3 = clock::now();
#if 0 
  float alpha, beta, gamma;
  int r3 = BVHNodeD_locatePoint(rd, X[0], X[1], alpha, beta, gamma);
  auto t4 = clock::now();
  int r4 = BVHNodeD_locatePoint_recursive(rd, rd, X[0], X[1], alpha, beta, gamma);
  auto t5 = clock::now();
  fprintf(stderr, "r0=%d, r1=%d, r2=%d, r3=%d, r4=%d\n", r0, r1, r2, r3, r4);
#endif
  fprintf(stderr, "x={%f, %f}, r0=%d, r1=%d, r2=%d\n", X[0], X[1], r0, r1, r2);

  float tt0 = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
  float tt1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
  float tt2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count();
  // float tt3 = std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count();
  // float tt4 = std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count();
  // fprintf(stderr, "tt0=%f, tt1=%f, tt2=%f, tt3=%f, tt4=%f\n", tt0, tt1, tt2, tt3, tt4);
  fprintf(stderr, "tt0=%f, tt1=%f, tt2=%f\n", tt0, tt1, tt2);
#endif
  deleteBVH(root);

  return rd;
}

