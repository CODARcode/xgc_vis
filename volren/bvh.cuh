#ifndef _BVH_CUH
#define _BVH_CUH

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


#endif
