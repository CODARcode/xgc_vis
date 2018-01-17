#ifndef _XGC_BLOB_EXTRACTOR_H
#define _XGC_BLOB_EXTRACTOR_H

#include <tourtre.h>
#include <cmath>
#include <vector>
#include <set>
#include <map>

class XGCBlobExtractor
{
public:
  XGCBlobExtractor(int nNodes, int nTriangles, int nPhi, double *coords, int *conn); 
  ~XGCBlobExtractor() {}
 
  void setData(double *dpot);
  void setPersistenceThreshold(double threshold) {persistenceThreshold = threshold;}

  std::vector<int>& getLabels(int plane) {return all_labels[plane];}

  void dumpMesh(const std::string& filename);
  void dumpLabels(const std::string& filename);
  void dumpBranchDecompositions(const std::string& filename);

public: 
  void buildContourTree2D(int plane);
  void simplifyBranchDecompositionByThreshold(ctBranch *b, double threshold, void *);
  void simplifyBranchDecompositionByNumbers(ctBranch* rootBranch, std::map<ctBranch*, size_t> &branchSet, int nLimit, void *d);
  void buildSegmentation(ctBranch *b, std::vector<int> &labels, void*); 

  void buildContourTree3D(); 
  void buildSegmentation3D(ctBranch *b, std::vector<int> &labels, void*); 

  void extractExtremum(int plane, double *dpot);
  void constructDiscreteGradient(double *dpot);

  int flood2D(size_t seed, int id, std::vector<int> &labels, double min, double max, void *d);
  void extractStreamers(int plane, ctBranch *root, std::map<ctBranch*, size_t>& branchSet, int nStreamers, double percentage, void *d);

private: // mesh
  std::vector<double> coords;
  std::vector<int> conn;
  // double *coords; 
  // int *conn;
  int nNodes, nTriangles, nPhi;

private: // data
  double *dpot;

private: // parameters
  double persistenceThreshold;

private: // analysis
  std::vector<std::set<int> > nodeGraph; // node->{neighbor nodes}
  std::map<int, std::vector<int> > maximum, minimum;
  std::map<int, std::vector<int> > all_labels;
  std::map<ctBranch*, size_t> branchSet;
}; 

#endif
