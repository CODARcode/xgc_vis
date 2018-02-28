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
  int nNodes, nTriangles, nPhi;
  int *conn = NULL, *nextNode = NULL;
  int *neighbors = NULL;
  double *coords = NULL;
  double *psi = NULL;
  float *psif = NULL;
  float *dispf = NULL; // displacement derived from nextNode
  float *invdetf = NULL; // inversed determinant of triangles
  float psi_min = FLT_MAX, psi_max = -FLT_MAX;

  std::vector<std::set<size_t> > nodeGraph; // node->{neighbor nodes}

  void readMeshFromADIOS(const std::string& filename, ADIOS_READ_METHOD readMethod, MPI_Comm comm);

  using json = nlohmann::json;
  json jsonfyMesh() const;
  json jsonfyMeshInfo() const;

  void deriveSinglePrecisionPsi();
  void deriveInversedDeterminants();
  void deriveDisplacements();

  void buildNeighbors();
  void buildNodeGraph();

  std::list<std::list<double> > marchingTriangles(double*, double isoval);

  std::vector<double> sampleScalarsAlongPsiContour(double *scalar, int nSamples, double isoval);

  ~XGCMesh();
};

#endif
