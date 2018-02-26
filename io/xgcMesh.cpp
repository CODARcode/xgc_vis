#include "def.h"
#include <unistd.h>
#include <mpi.h>
#include <adios_error.h>
#include <cassert>
#include <cfloat>
#include <iostream>
#include <map>
#include "io/xgcMesh.h"
#include "io/bp_utils.hpp"

void XGCMesh::readMeshFromADIOS(const std::string& filename, ADIOS_READ_METHOD readMethod, MPI_Comm comm)
{
  adios_read_init_method(readMethod, comm, "");
  ADIOS_FILE *fp = NULL;
  while (1) {
    fp = adios_read_open_file(filename.c_str(), ADIOS_READ_METHOD_BP, comm); // always use ADIOS_READ_METHOD_BP for mesh
    if (fp == NULL) {
      fprintf(stderr, "failed to open mesh: %s, will retry in 5 seconds.\n", filename.c_str()); 
      sleep(5);
    } else break;
  }
  adios_read_bp_reset_dimension_order(fp, 0);
  
  readValueInt(fp, "n_n", &nNodes);
  readValueInt(fp, "n_t", &nTriangles);
  readScalars<double>(fp, "/coordinates/values", &coords);
  readScalars<int>(fp, "/cell_set[0]/node_connect_list", &conn);
  readScalars<int>(fp, "nextnode", &nextNode);
  readScalars<double>(fp, "psi", &psi);
 
  adios_read_finalize_method (ADIOS_READ_METHOD_BP);
  adios_read_close(fp);
  
  deriveSinglePrecisionPsi();
  deriveDisplacements();
  deriveInversedDeterminants();
}

void XGCMesh::deriveSinglePrecisionPsi() {
  psif = (float*)realloc(psif, sizeof(float)*nNodes);
  for (int i=0; i<nNodes; i++) {
    psif[i] = psi[i];
    psi_min = std::min(psi_min, psif[i]);
    psi_max = std::max(psi_max, psif[i]);
  }
  // fprintf(stderr, "psi_min=%f, psi_max=%f\n", psi_min, psi_max);
}

void XGCMesh::deriveInversedDeterminants() {
  invdetf = (float*)realloc(invdetf, sizeof(float)*nTriangles);
  for (int i=0; i<nTriangles; i++) {
    const int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    float x0 = coords[i0*2], x1 = coords[i1*2], x2 = coords[i2*2],
          y0 = coords[i0*2+1], y1 = coords[i1*2+1], y2 = coords[i2*2+1];
    float det = (y1-y2)*(x0-x2) + (x2-x1)*(y0-y2);
    invdetf[i] = 1.f / det;
  }
}

void XGCMesh::deriveDisplacements() {
  dispf = (float*)realloc(dispf, sizeof(float)*nNodes*2);
  for (int i=0; i<nNodes; i++) {
    int j = nextNode[i];
    dispf[i*2] = coords[j*2] - coords[i*2];
    dispf[i*2+1] = coords[j*2+1] - coords[i*2+1];
    // fprintf(stderr, "%d, %d, %f, %f\n", i, j, dispf[i*2], dispf[i*2+1]);
  }
}

void XGCMesh::buildNeighbors()
{
  typedef std::tuple<int, int> edgeType;
  std::map<edgeType, int> edgeTriangleMap;

  auto makeEdge = [](int i0, int i1) {return std::make_tuple(i0, i1);};

  auto addEdge = [&makeEdge, &edgeTriangleMap](int triangleId, int i0, int i1) {
    // fprintf(stderr, "triangleId=%d, adding edge {%d, %d}\n", triangleId, i0, i1);
    edgeTriangleMap[ makeEdge(i0, i1) ] = triangleId;
  };

  for (int i=0; i<nTriangles; i++) {
    const int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    addEdge(i, i0, i1);
    addEdge(i, i1, i2);
    addEdge(i, i2, i0);
  }

  neighbors = (int*)realloc(neighbors, nTriangles*3*sizeof(int));

  for (int i=0; i<nTriangles; i++) {
    const int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    const edgeType edges[3] = {makeEdge(i0, i2), makeEdge(i2, i1), makeEdge(i1, i0)};

    for (int j=0; j<3; j++) {
      const edgeType& e = edges[j];
      if (edgeTriangleMap.find(e) != edgeTriangleMap.end()) neighbors[i*3+j] = edgeTriangleMap[e];
      else neighbors[i*3+j] = -1;
    }

    // fprintf(stderr, "triangleId=%d, neighbors={%d, %d, %d}\n", i, 
    //     neighbors[i*3], neighbors[i*3+1], neighbors[i*3+2]);
  }
}

void XGCMesh::buildNodeGraph()
{
  nodeGraph.clear();
  nodeGraph.resize(nNodes);
  for (int i=0; i<nTriangles; i++) {
    int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];

    nodeGraph[i0].insert(i1);
    nodeGraph[i0].insert(i2);
    nodeGraph[i1].insert(i0);
    nodeGraph[i1].insert(i2);
    nodeGraph[i2].insert(i0);
    nodeGraph[i2].insert(i1);
  }
}

void XGCMesh::marchingTriangles(double *scalar, double isoval)
{
  auto checkZero = [scalar, isoval](int i0, int i1, double &alpha) {
    double f0 = scalar[i0], f1 = scalar[i1];
    alpha = (isoval - f0) / (f1 - f0);
    return alpha >= 0 && alpha < 1;
  };

  struct IntersectedTriangle {
    double alpha0, alpha1, alpha2;
  };

  std::set<IntersectedTriangle> interstectedTriangles;

  double alpha; 
  for (int i=0; i<nTriangles; i++) {
    int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    checkZero(i0, i1, alpha); 
  }
}

XGCMesh::~XGCMesh() {
  free(conn);
  free(psi);
  free(nextNode);
  free(coords);
  free(dispf);
  free(invdetf);
}
