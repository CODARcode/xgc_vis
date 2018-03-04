#include "kdbvh.h"

void buildKDBVHRecursively(const XGCMesh &m, KDBVHNode *node, std::vector<int> &tids, int dir)
{
  size_t size = tids.size();
  size_t half_size = size/2; 
 
  if (size > 1) {
    // sort triangles by centroids
    std::stable_sort(tids.begin(), tids.end(), [&m, dir](int t0, int t1) { return m.centroidsf[t0*2+dir] < m.centroidsf[t1*2+dir]; });
    std::vector<int> lo(tids.begin(), tids.begin() + half_size);
    std::vector<int> hi(tids.begin() + half_size, tids.end());

    // build children
    const int split_dir = (dir+1)%2;

    KDBVHNode *l = new KDBVHNode, *h = new KDBVHNode;
    buildKDBVHRecursively(m, l, lo, split_dir);
    buildKDBVHRecursively(m, h, hi, split_dir);

    // update AABB
    node->aabb.A[0] = std::min(l->aabb.A[0], h->aabb.A[0]);
    node->aabb.A[1] = std::min(l->aabb.A[1], h->aabb.A[1]);
    node->aabb.B[0] = std::max(l->aabb.B[0], h->aabb.B[0]);
    node->aabb.B[1] = std::max(l->aabb.B[1], h->aabb.B[1]);
 
    node->lchild = l; 
    node->hchild = h;
  } else { // size=1
    fprintf(stderr, "leaf node, tid=%d\n", tids[0]);
   
    const int i = tids[0];
    const int i0 = m.conn[i*3], i1 = m.conn[i*3+1], i2 = m.conn[i*3+2];
    float x0 = m.coords[i0*2], x1 = m.coords[i1*2], x2 = m.coords[i2*2],
          y0 = m.coords[i0*2+1], y1 = m.coords[i1*2+1], y2 = m.coords[i2*2+1];

    node->aabb.A[0] = min3(x0, x1, x2);
    node->aabb.B[0] = max3(x0, x1, x2);
    node->aabb.A[1] = min3(y0, y1, y2);
    node->aabb.B[1] = max3(y0, y1, y2);
  }
}

void buildKDBVH(const XGCMesh &m)
{
  std::vector<int> tids(m.nNodes);
  for (int i=0; i<m.nNodes; i++) tids[i] = i;

  KDBVHNode *root = new KDBVHNode;
  buildKDBVHRecursively(m, root, tids);
}
