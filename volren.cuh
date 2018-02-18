struct QuadNodeD {
  // tree
  int parentId;
  int childrenIds[4];
  bool isLeaf;

  // bounds
  float Ax, Ay, Bx, By;

  // triangle
  int i0, i1, i2;
  float x0, y0, x1, y1, x2, y2;
};
