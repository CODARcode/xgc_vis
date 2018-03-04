#ifndef _KDBVH_H
#define _KDBVH_H

#include "io/xgcMesh.h"
#include "volren/aabb.h"

struct KDBVHNode {
  KDBVHNode *parent = NULL, *lchild = NULL, *hchild = NULL;
  int split_dir = 0; // 0: x, 1: y
  float split_plane = 0.f;
  AABB aabb;
};

void buildKDBVHRecursively(const XGCMesh &m, KDBVHNode *n, std::vector<int> &tids, int dir=0);

void buildKDBVH(const XGCMesh &m);

#endif
