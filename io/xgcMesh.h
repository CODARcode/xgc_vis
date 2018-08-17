#ifndef _XGCMESH_H
#define _XGCMESH_H

#include "def.h"
#include <mpi.h>
#include <adios.h>
#include <adios_read.h>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <cfloat>
#include <json.hpp>

struct XGCMesh {
  int nNodes, nTriangles, nPhi, iPhi;
  int *conn = NULL, *nextNode = NULL;
  int *neighbors = NULL;
  double *coords = NULL;
  double *psi = NULL;
  float *centroidsf = NULL;
  float *psif = NULL;
  float *dispf = NULL; // displacement derived from nextNode
  float *invdetf = NULL; // inversed determinant of triangles
  float psi_min = FLT_MAX, psi_max = -FLT_MAX;
  int psi_min_node; 
  float psi_min_x, psi_min_y;
  double coords_min_x, coords_min_y, coords_max_x, coords_max_y, coords_centroid_x, coords_centroid_y;

#if WITH_ADIOS
  void readMeshFromADIOS(const std::string& filename, ADIOS_READ_METHOD readMethod, MPI_Comm comm);
#endif

#if WITH_H5
  void readMeshFromH5(const std::string& filename);
#endif

  int findPsiSaddle(); 

  using json = nlohmann::json;
  json jsonfyMesh() const;
  json jsonfyMeshInfo() const;

  std::set<size_t> getNodeNeighbors2D(size_t i) const {return nodeNeighbors[i];}
  std::set<size_t> getNodeNeighbors3D(size_t v) const {
    size_t plane = v / nNodes; 
    size_t node = v % nNodes;
    const std::set<size_t> &neighbors2D = nodeNeighbors[node];

    std::set<size_t> neighbors;
    for (auto &neighbor2D : neighbors2D) {
      neighbors.insert( ((plane - 1 + nPhi) % nPhi) * nNodes + neighbor2D);
      neighbors.insert( plane * nNodes + neighbor2D );
      neighbors.insert( ((plane + 1 + nPhi) % nPhi) * nNodes + neighbor2D);
    }
    return neighbors;
  }

  std::list<std::list<double> > marchingTriangles(double*, double isoval);
  std::vector<double> sampleScalarsAlongPsiContour(double *scalar, int nSamples, double isoval);
  std::vector<double> testMarchingTriangles(double *scalar, double isoval);

  ~XGCMesh();
  
private:
  void preprocessMesh();

  void buildNeighbors();
  void buildNodeGraph();

  std::vector<std::set<size_t> > nodeNeighbors; // nodeGraph; // node->{neighbor nodes}
};

#endif
