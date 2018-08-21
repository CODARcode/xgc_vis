#ifndef _XGC_BLOB_EXTRACTOR_H
#define _XGC_BLOB_EXTRACTOR_H

#include <tourtre.h>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <json.hpp>
// #include <ftk/transition/transition.h>
#include <ftk/graph/graph.hh>

class XGCBlobExtractor
{
  using json = nlohmann::json;

public:
  XGCBlobExtractor(int nNodes, int nTriangles, double *coords, int *conn); 
  ~XGCBlobExtractor() {}

  int getNNodes() const {return nNodes;}
  int getNPhi() const {return nPhi;}
  size_t getTimestep() const {return timestep;}

  void setData(size_t timestep, int nPhi, double *dpot);
  const double *getData() const {return dpot;}
  void setPersistenceThreshold(double threshold) {persistenceThreshold = threshold;}

  std::vector<int>& getLabels(int plane) {return all_labels[plane];}
  std::vector<int> getFlattenedLabels(int plane);

  json jsonfyMesh() const;
  json jsonfyBranches(std::map<ctBranch*, size_t> &branchSet, size_t top=0); 

  void dumpLabels(const std::string& filename);
  void dumpBranches(const std::string& filename, std::map<ctBranch*, size_t>&, size_t top=0);

public:
  void buildContourTree2DAll();

  std::map<ctBranch*, size_t> buildContourTree2D(int plane);
  void simplifyBranchDecompositionByThreshold(ctBranch *b, double threshold, void *);
  void simplifyBranchDecompositionByNumbers(ctBranch* rootBranch, std::map<ctBranch*, size_t> &branchSet, int nLimit, void *d);
  void buildSegmentation(ctBranch *b, std::vector<int> &labels, void*); 

  void buildContourTree3D(); 
  void buildSegmentation3D(ctBranch *b, std::vector<int> &labels, void*); 

  void extractExtremum(int plane, double *dpot);
  void constructDiscreteGradient(double *dpot);

  int flood2D(size_t seed, int id, std::vector<int> &labels, double min, double max, void *d);
  void extractStreamers(int plane, ctBranch *root, std::map<ctBranch*, size_t>& branchSet, int nStreamers, double percentage, void *d);

public:
  void findConnectedComponents(
      const std::vector<size_t> &totalOrder, 
      size_t i, 
      std::set<size_t> &cc,
      size_t &lowest);
  void buildJoinTree(
      const std::vector<size_t> &totalOrder, 
      const std::vector<size_t> &ranks);

private: // mesh
  std::vector<double> coords;
  std::vector<int> conn;
  // double *coords; 
  // int *conn;
  int nNodes, nTriangles, nPhi;

private: // data
  size_t timestep;
  double *dpot;

private: // parameters
  double persistenceThreshold;

private: // analysis
  std::vector<std::set<size_t> > nodeGraph; // node->{neighbor nodes}
  std::map<int, std::vector<int> > maximum, minimum;
  std::map<int, std::vector<int> > all_labels, all_signs;
  // std::map<ctBranch*, size_t> branchSet;

public:
  // FeatureTransitionMatrix relateFeatures(const std::vector<int> &labels0, const std::vector<int> &signs0, 
  //     const std::vector<int> &labels1, const std::vector<int> &sings1);
}; 

#endif
