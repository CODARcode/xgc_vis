#ifndef _XGC_BLOB_EXTRACTOR_H
#define _XGC_BLOB_EXTRACTOR_H

#include <tourtre.h>
#include <cmath>
#include <set>
#include <map>

class XGCBlobExtractor
{
public:
  XGCBlobExtractor() {}
  ~XGCBlobExtractor() {}
 
  void setMesh(int nNodes, int nTriangles, int nPhi, double *coords, int *conn); 
  void setData(double *dpot);

public: 
  void buildContourTree2D(int plane);
  void simplifyBranchDecompositionByThreshold(ctBranch *b, double threshold, void *);
  void simplifyBranchDecompositionByNumbers(ctBranch *b, int nLimit, void *);
  void buildSegmentation(ctBranch *b, std::vector<size_t> &labels, void*); 

  void buildContourTree3D(); 
  void buildSegmentation3D(ctBranch *b, std::vector<size_t> &labels, void*); 

  void extractExtremum(int plane, double *dpot);
  void constructDiscreteGradient(double *dpot);

private: // mesh
  double *coords; 
  int *conn;
  int nNodes, nTriangles, nPhi;

private: // data
  double *dpot;

private: // analysis
  std::vector<std::set<int> > nodeGraph; // node->{neighbor nodes}
  std::map<int, std::vector<int> > maximum, minimum;
  std::map<int, std::vector<size_t> > all_labels;
}; 

#endif
