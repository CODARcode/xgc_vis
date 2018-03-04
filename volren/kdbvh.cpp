#include "kdbvh.h"
#include <stack>

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
    // fprintf(stderr, "leaf node, tid=%d\n", tids[0]);
    const int i = tids[0];
    const int i0 = m.conn[i*3], i1 = m.conn[i*3+1], i2 = m.conn[i*3+2];
    float x0 = m.coords[i0*2], x1 = m.coords[i1*2], x2 = m.coords[i2*2],
          y0 = m.coords[i0*2+1], y1 = m.coords[i1*2+1], y2 = m.coords[i2*2+1];

    node->aabb.A[0] = min3(x0, x1, x2);
    node->aabb.B[0] = max3(x0, x1, x2);
    node->aabb.A[1] = min3(y0, y1, y2);
    node->aabb.B[1] = max3(y0, y1, y2);

    node->tid = i;
  }
}

void deleteKDBVHNodeRecursively(KDBVHNode *n)
{
  delete n->lchild; 
  delete n->hchild;
  delete n;
}

std::vector<BVHNodeD> buildKDBVHGPU(const XGCMesh &m)
{
  std::vector<int> tids(m.nNodes);
  for (int i=0; i<m.nNodes; i++) tids[i] = i;

  KDBVHNode *r = new KDBVHNode;
  buildKDBVHRecursively(m, r, tids);

  std::map<KDBVHNode*, int> nodeMap; 
  std::map<int, KDBVHNode*> nodeReverseMap;
  std::stack<KDBVHNode*> S;
  S.push(r);

  int kdNodeCount = 0;
  int maxStackSize = 0;

  while (!S.empty()) {
    maxStackSize = std::max(maxStackSize, static_cast<int>(S.size()));

    KDBVHNode *n = S.top();
    S.pop(); 

    int nodeId = kdNodeCount ++;
    nodeMap[n] = nodeId;
    nodeReverseMap[nodeId] = n;

    if (n->lchild) S.push(n->lchild);
    if (n->hchild) S.push(n->hchild);
  }

  std::vector<BVHNodeD> rd(kdNodeCount);

  for (int i=0; i<kdNodeCount; i++) {
    KDBVHNode n = *nodeReverseMap[i];
    BVHNodeD &d = rd[i];

    // parent
    if (n.parent == NULL) d.parentId = -1;
    else d.parentId = nodeMap[n.parent];

    // children
    d.childrenIds[0] = nodeMap[n.lchild];
    d.childrenIds[1] = nodeMap[n.hchild];
    
    // bounds
    d.Ax = n.aabb.A[0];
    d.Ay = n.aabb.A[1];
    d.Bx = n.aabb.B[0];
    d.By = n.aabb.B[1];
    // fprintf(stderr, "%f, %f, %f, %f\n", d.Ax, d.Ay, d.Bx, d.By);

    if (n.lchild == NULL || n.hchild == NULL) {
      const int id = n.tid; 
      d.triangleId = id;
      d.i0 = m.conn[id*3];
      d.i1 = m.conn[id*3+1];
      d.i2 = m.conn[id*3+2];
      d.x0 = m.coords[d.i0*2];
      d.y0 = m.coords[d.i0*2+1];
      d.x1 = m.coords[d.i1*2];
      d.y1 = m.coords[d.i1*2+1];
      d.x2 = m.coords[d.i2*2];
      d.y2 = m.coords[d.i2*2+1];
    } else 
      d.triangleId = -1;
  }

  fprintf(stderr, "nNodes=%d, kdNodeCount=%d, maxStackSize=%d\n", m.nNodes, kdNodeCount, maxStackSize);
  deleteKDBVHNodeRecursively(r);

  return rd;
}
