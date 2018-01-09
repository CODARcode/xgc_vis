#include <fstream>
#include <iostream>
#include <queue>
#include "xgcBlobExtractor.h"
#include "json.hpp"

using json = nlohmann::json;

XGCBlobExtractor::XGCBlobExtractor(int nNodes_, int nTriangles_, int nPhi_, double *coords_, int *conn_) :
  nNodes(nNodes_), nTriangles(nTriangles_), nPhi(nPhi_), 
  coords(coords_, coords_ + nNodes_*2), 
  conn(conn_, conn_ + 3*nTriangles_)
{
  // init nodeGraph
  nodeGraph.clear();
  nodeGraph.resize(nNodes);
  for (int i=0; i<nTriangles; i++) {
    int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    nodeGraph[i0].insert(i1);
    nodeGraph[i1].insert(i2);
    nodeGraph[i2].insert(i0);
  }
}

void XGCBlobExtractor::constructDiscreteGradient(double *dpot) // based on MS theory
{
  // 0-cell: vertices
  // 1-cell: edges
  // 2-cell: cells
  
  typedef struct {
    unsigned char index;
    double value;
    std::vector<int> hCells;
  } msCell;

}

struct XGCData { // either single slice or all slices
  std::vector<std::set<int> > *nodeGraph;
  double *dpot;
  int nNodes, nPhi;
};

static double value(size_t v, void *d)
{
  XGCData *data = (XGCData*)d;
  // fprintf(stderr, "%d, %f\n", v, data->dpot[v]);
  return data->dpot[v];
}

static size_t neighbors(size_t v, size_t *nbrs, void *d)
{
  XGCData *data = (XGCData*)d;
  std::set<int>& nodes = (*data->nodeGraph)[v];
  
  int i = 0;
  for (std::set<int>::iterator it = nodes.begin(); it != nodes.end(); it ++) {
    nbrs[i] = *it;
    i ++;
  }

  return i;
}

static size_t neighbors3D(size_t v, size_t *nbrs, void *d)
{
  XGCData *data = (XGCData*)d;
 
  const size_t nNodes = data->nNodes;
  const size_t nPhi = data->nPhi;

  size_t plane = v / nNodes;
  size_t node = v % nNodes;

  // fprintf(stderr, "nnodes=%d, nphi=%d, plane=%d, node=%d\n", nNodes, nPhi, plane, node);

  std::set<int>& nodes2D = (*data->nodeGraph)[node];
  
  int i = 0;
  for (std::set<int>::iterator it = nodes2D.begin(); it != nodes2D.end(); it ++) {
    size_t nbr2D = *it;
    nbrs[i++] = ((plane - 1 + nPhi) % nPhi) * nNodes + nbr2D;
    nbrs[i++] = plane * nNodes + nbr2D;
    nbrs[i++] = ((plane + 1 + nPhi) % nPhi) * nNodes + nbr2D;
  }

  return i;
}

static double volumePriority(ctNode *node, void *d)
{
  XGCData *data = static_cast<XGCData*>(d);
  ctArc *arc = ctNode_leafArc(node);

  ctNode *lo = arc->lo,
         *hi = arc->hi;

  double val_lo = value(lo->i, d),
         val_hi = value(hi->i, d);

  int count = 0;

  std::queue<size_t> Q;
  std::set<size_t> visited;

  Q.push(lo->i);
  visited.insert(lo->i); 

  while (!Q.empty()) {
    size_t p = Q.front(); 
    Q.pop();
    double val_p = value(p, d);

    std::set<int> &neighbors = (*data->nodeGraph)[p];
    for (std::set<int>::iterator it = neighbors.begin(); it != neighbors.end(); it ++) {
      const int neighbor = *it; 

      if (visited.find(neighbor) == visited.end()) { // not found
        double val_q = value(neighbor, d);
        if (val_q >= val_lo && val_q < val_hi) {
          Q.push(neighbor);
          visited.insert(neighbor);
          count ++;
        }
      }
    }
  }

  fprintf(stderr, "volume=%d\n", count);
  return count;
}

static void printContourTree(ctBranch* b)
{
  fprintf(stderr, "%d, %d, %p\n", b->extremum, b->saddle, b->children.head);
  
  for (ctBranch* c = b->children.head; c != NULL; c = c->nextChild) 
    printContourTree(c);
}

json branch2json(ctBranch* b, std::map<ctBranch*, size_t> &branchSet, void *d)
{
  json j;

  j["id"] = branchSet[b];
  j["extremum"] = b->extremum;
  j["extremum_val"] = value(b->extremum, d);
  j["saddle"] = b->saddle;
  j["saddle_val"] = value(b->saddle, d);
  if (branchSet.find(b->parent) != branchSet.end()) 
    j["parent"] = branchSet[b->parent];

  std::vector<size_t> children;
  for (ctBranch *c = b->children.head; c != NULL; c = c->nextChild)
    children.push_back(branchSet[c]);
  if (children.size() > 0)
    j["children"] = children;

  if (branchSet.find(b->nextChild) != branchSet.end()) 
    j["nextChild"] = branchSet[b->nextChild];
  
  if (branchSet.find(b->prevChild) != branchSet.end()) 
    j["prevChild"] = branchSet[b->prevChild];

  return j;
}

void unserializeBranch(const std::string& str)
{

}

void XGCBlobExtractor::simplifyBranchDecompositionByThreshold(ctBranch *b, double threshold, void *d)
{
  // XGCData *data = (XGCData*)d;

RESTART: 
  for (ctBranch *c=b->children.head; c!=NULL; c=c->nextChild)
    if (fabs(value(c->extremum, d) - value(c->saddle, d)) < threshold) {
      ctBranchList_remove(&b->children, c);
      goto RESTART;
    } else {
      simplifyBranchDecompositionByThreshold(c, threshold, d);
    }
}

#if 0
void XGCBlobExtractor::simplifyBranchDecompositionByNumbers(ctBranch *b, int nLimit, void *d) // no more than a number
{
  float minPriority; // default is persistence
  ctBranch *minPriorityBranch = NULL; 

  qDebug() << "calculating priorities..."; 

  QList< QPair<float, ctBranch*> > m_branchPriorityList; 
  foreach (ctBranch *c, m_reverseBranchMap.keys()) {
    // m_branchPriorityList << qMakePair(static_cast<float>(volumePriority(c, blk)), c); 
    m_branchPriorityList << qMakePair(static_cast<float>(fabs(value(c->extremum, blk) - value(c->saddle, blk))), c); 
  }
  qStableSort(m_branchPriorityList); 

  qDebug() << "simplifying branch decomposition..."; 
  while (minPriorityBranch != m_rootBranch && m_reverseBranchMap.size() > nLimit) {
    minPriority = 1e38f; 
    minPriorityBranch = m_rootBranch; 

    for (QList< QPair<float, ctBranch*> >::iterator itor=m_branchPriorityList.begin(); itor!=m_branchPriorityList.end(); itor++) {
      if (itor->second->children.head == NULL) {
        minPriority = itor->first; 
        minPriorityBranch = itor->second; 
        m_branchPriorityList.erase(itor); 
        break; 
      }
    }

    qDebug() << "priority =" << minPriority << "branch =" << minPriorityBranch; 

    if (minPriorityBranch == m_rootBranch) break; 
    ctBranch *parentBranch = minPriorityBranch->parent; 

    ctBranchList_remove(&parentBranch->children, minPriorityBranch); 
    m_reverseBranchMap.remove(minPriorityBranch); 
  }
}
#endif

void XGCBlobExtractor::buildContourTree2D(int plane)
{
  fprintf(stderr, "building contour tree for plane %d, nNodes=%d\n", plane, nNodes);

  std::vector<size_t> totalOrder; 
  for (int i=0; i<nNodes; i++) 
    totalOrder.push_back(i);

  XGCData data;
  data.nodeGraph = &nodeGraph;
  data.dpot = dpot + plane * nNodes;
  data.nNodes = nNodes;
  data.nPhi = nPhi;

  std::sort(totalOrder.begin(), totalOrder.end(),
      [&data](size_t v0, size_t v1) {
        return data.dpot[v0] < data.dpot[v1];
        // if (fabs(data.dpot[v0] - data.dpot[v1]) < 1e-5) return v0 < v1; 
        // else return data.dpot[v0] < data.dpot[v1];
      });

  ctContext *ctx = ct_init(
      nNodes, 
      &totalOrder.front(),  
      &value,
      &neighbors, 
      &data);

  // ct_priorityFunc(ctx, volumePriority);

  ct_sweepAndMerge(ctx);
  ctBranch *root = ct_decompose(ctx);
  ctBranch **map = ct_branchMap(ctx);

  simplifyBranchDecompositionByThreshold(root, 20, &data);

  std::vector<size_t> labels(nNodes, 0);
  buildSegmentation(root, labels, &data);
  // for (int i=0; i<nNodes; i++) 
  //   fprintf(stderr, "%d, %d\n", i, labels[i]);

  all_labels[plane] = labels;

  // printContourTree(root);

  ct_cleanup(ctx);
}

void XGCBlobExtractor::buildContourTree3D()
{
  fprintf(stderr, "building contour tree over 3D...\n");

  std::vector<size_t> totalOrder; 
  for (int i=0; i<nNodes*nPhi; i++) 
    totalOrder.push_back(i);

  XGCData data;
  data.nodeGraph = &nodeGraph;
  data.dpot = dpot;
  data.nNodes = nNodes;
  data.nPhi = nPhi;
  
  std::sort(totalOrder.begin(), totalOrder.end(),
      [&data](size_t v0, size_t v1) {
        return data.dpot[v0] < data.dpot[v1];
        // if (fabs(data.dpot[v0] - data.dpot[v1]) < 1e-5) return v0 < v1; 
        // else return data.dpot[v0] < data.dpot[v1];
      });
  
  ctContext *ctx = ct_init(
      nNodes * nPhi, 
      &totalOrder.front(),  
      &value,
      &neighbors3D, 
      &data);
  
  ct_sweepAndMerge(ctx);
  ctBranch *root = ct_decompose(ctx);
  ctBranch **branchMap = ct_branchMap(ctx);
  ctArc **arcMap = ct_arcMap(ctx);

#if 0
  std::map<ctArc*, size_t> arcSet;
  size_t nArcs = 0;
  for (int i=0; i<nNodes*nPhi; i++) {
    if (arcSet.find(arcMap[i]) == arcSet.end()) {
      size_t arcId = nArcs ++;
      // arcSet.insert(arcMap[i], arcId);
      arcSet[arcMap[i]] = arcId;
      fprintf(stderr, "arcId=%d, %p, hi=%p, lo=%p\n", arcId, arcMap[i], arcMap[i]->hi, arcMap[i]->lo);
      label_colors[arcId] = QColor(rand()%256, rand()%256, rand()%256);
    }
  }
 
  for (int i=0; i<nPhi; i++) {
    std::vector<size_t> labels(nNodes, 0);
    for (int j=0; j<nNodes; j++) {
      size_t arcId = arcSet[arcMap[i*nNodes + j]];
      labels[j] = arcId;
    }
    all_labels[i] = labels;
  }
#else
  size_t nBranches = 0;
  for (int i=0; i<nNodes*nPhi; i++) {
    if (branchSet.find(branchMap[i]) == branchSet.end()) {
      size_t branchId = nBranches ++;
      branchSet[branchMap[i]] = branchId;
      // fprintf(stderr, "branchId=%d, %p\n", branchId, branchMap[i]);
      // label_colors[branchId] = QColor(rand()%256, rand()%256, rand()%256);
    }
  }
 
  for (int i=0; i<nPhi; i++) {
    std::vector<size_t> labels(nNodes, 0);
    for (int j=0; j<nNodes; j++) {
      size_t branchId = branchSet[branchMap[i*nNodes + j]];
      labels[j] = branchId;
    }
    all_labels[i] = labels;
  }
#endif

#if 1 // simplifications
  simplifyBranchDecompositionByThreshold(root, 20, &data);
  
  std::vector<size_t> labels(nNodes*nPhi, 0);
  buildSegmentation3D(root, labels, &data);
  // TODO

  for (int i=0; i<nPhi; i++) {
    std::vector<size_t> l(nNodes);
    std::copy(labels.begin() + i*nNodes, labels.begin() + (i+1)*nNodes - 1, l.begin());
    all_labels[i] = l;
  }
  
  ct_cleanup(ctx); 
  // TODO
#endif

  fprintf(stderr, "built contour tree over 3D.\n");
}

void XGCBlobExtractor::buildSegmentation3D(ctBranch *b, std::vector<size_t> &labels, void *d)
{
  static int id = 1; // FIXME
  const int currentId = ++ id;

  // label_colors[currentId] = QColor(rand()%256, rand()%256, rand()%256);

  XGCData *data = (XGCData*)d;

  const double val_extremum = value(b->extremum, d), 
               val_saddle = value(b->saddle, d);
  const double val_hi = std::max(val_extremum, val_saddle),
               val_lo = std::min(val_extremum, val_saddle);

  int count = 1;

  std::queue<size_t> Q;
  Q.push(b->extremum);
  labels[b->extremum] = currentId; 

  while (!Q.empty()) {
    size_t p = Q.front();
    Q.pop();

    size_t nbrs[128]; // TODO
    size_t nnbrs = neighbors3D(p, nbrs, d);

    for (int i=0; i<nnbrs; i++) {
      const size_t neighbor = nbrs[i];

      if (labels[neighbor] != currentId) {
        double val = value(neighbor, d);
        bool qualify = true;
        if (val_extremum > val_saddle) {
          if (val > val_extremum || val <= val_saddle) qualify = false;
        } else {
          if (val < val_extremum || val >= val_saddle) qualify = false;
        }
        
        if (qualify) {
          Q.push(neighbor);
          labels[neighbor] = currentId;
          count ++;
        }
      }
    }
  }

  // fprintf(stderr, "#count=%d\n", count);

  for (ctBranch *c = b->children.head; c != NULL; c = c->nextChild) 
    buildSegmentation3D(c, labels, d);
}

void XGCBlobExtractor::buildSegmentation(ctBranch *b, std::vector<size_t> &labels, void *d)
{
  static int id = 1; // FIXME
  const int currentId = ++ id;

  // label_colors[currentId] = QColor(rand()%256, rand()%256, rand()%256);

  XGCData *data = (XGCData*)d;

  const double val_extremum = value(b->extremum, d), 
               val_saddle = value(b->saddle, d);
  const double val_hi = std::max(val_extremum, val_saddle),
               val_lo = std::min(val_extremum, val_saddle);

  int count = 1;

  std::queue<size_t> Q;
  Q.push(b->extremum);
  labels[b->extremum] = currentId; 

  while (!Q.empty()) {
    size_t p = Q.front();
    Q.pop();

    std::set<int> &neighbors = nodeGraph[p];
    for (std::set<int>::iterator it = neighbors.begin(); it != neighbors.end(); it ++) {
      const int neighbor = *it;

      if (labels[neighbor] != currentId) {
        double val = value(neighbor, d);
        bool qualify = true;
        if (val_extremum > val_saddle) {
          if (val > val_extremum || val <= val_saddle) qualify = false;
        } else {
          if (val < val_extremum || val >= val_saddle) qualify = false;
        }
        
        if (qualify) {
          Q.push(neighbor);
          labels[neighbor] = currentId;
          count ++;
        }
      }
    }
  }

  // fprintf(stderr, "#count=%d\n", count);

  for (ctBranch *c = b->children.head; c != NULL; c = c->nextChild) 
    buildSegmentation(c, labels, d);
}

void XGCBlobExtractor::extractExtremum(int plane, double *dpot_)
{
  double *dpot = dpot_ + plane * nNodes;

  for (int i=0; i<nNodes; i++) {
    std::set<int> &neighbors = nodeGraph[i];
    bool local_max = true, local_min = true;

    for (std::set<int>::iterator it = neighbors.begin(); it != neighbors.end(); it ++) {
      const int j = *it;
      if (dpot[i] >= dpot[j]) local_min = false;
      if (dpot[i] <= dpot[j]) local_max = false;
    }

    if (local_max)
      maximum[plane].push_back(i); 
    else if (local_min)
      minimum[plane].push_back(i);
  }
}

void XGCBlobExtractor::setData(double *dpot_)
{
  dpot = dpot_;
#if 0
  const float min = -100, max = 100;

  f_colors.clear();

  for (int plane = 0; plane < nPhi; plane ++) {
    for (int i=0; i<nTriangles; i++) {
      int v[3] = {conn[i*3], conn[i*3+1], conn[i*3+2]};

      for (int j=0; j<3; j++) {
        float val = clamp_normalize(min, max, (float)dpot[plane*nNodes + v[j]]);
        f_colors.push_back(val);
        f_colors.push_back(1-val);
        f_colors.push_back(0);
      }
    }
  }

  buildContourTree3D(dpot);

  for (int i=0; i<nPhi; i++) {
    // buildContourTree(i, dpot); 
    extractExtremum(i, dpot); // dpot + i*nNodes);
  }
#endif
}

void XGCBlobExtractor::dumpLabels(const std::string& filename) 
{
  FILE *fp = fopen(filename.c_str(), "wb");
  for (int i=0; i<nPhi; i++) 
    fwrite(all_labels[i].data(), sizeof(size_t), nNodes, fp);
  fclose(fp);
}


void XGCBlobExtractor::dumpBranchDecompositions(const std::string& filename) {
  XGCData data;
  data.nodeGraph = &nodeGraph;
  data.dpot = dpot;
  data.nNodes = nNodes;
  data.nPhi = nPhi;

  std::ofstream ofs(filename, std::ofstream::out);
  json jresults;
  for (std::map<ctBranch*, size_t>::iterator it = branchSet.begin();  it != branchSet.end();  it ++) {
    jresults[it->second] = branch2json(it->first, branchSet, &data);
  }
  ofs << jresults.dump();
  ofs.close();
}

void XGCBlobExtractor::dumpMesh(const std::string& filename) {
  std::ofstream ofs(filename, std::ofstream::out);
  json j;

  j["nPhi"] = nPhi;
  j["nNodes"] = nNodes;
  j["nTriangles"] = nTriangles;

  j["coords"] = coords;
  j["conn"] = conn;

  ofs << j.dump();
  ofs.close();
}
