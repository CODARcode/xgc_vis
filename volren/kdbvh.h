#ifndef _KDBVH_H
#define _KDBVH_H

#include "io/xgcMesh.h"
#include "volren/aabb.h"
#include "volren/bvh.cuh"

struct KDBVHNode {
  KDBVHNode *parent = NULL, *lchild = NULL, *hchild = NULL;
  int split_dir = 0; // 0: x, 1: y
  int tid = -1; // only valid for leaf node
  AABB aabb;
};

void buildKDBVHRecursively(const XGCMesh &m, KDBVHNode *n, std::vector<int> &tids, int dir=0);
void deleteKDBVHNodeRecursively(KDBVHNode *n);

std::vector<BVHNodeD> buildKDBVHGPU(const XGCMesh &m);

#endif
