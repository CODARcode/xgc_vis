#include "def.h"
#include <unistd.h>
#include <mpi.h>
#include <adios_error.h>
#include <cassert>
#include <cfloat>
#include <iostream>
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

XGCMesh::~XGCMesh() {
  free(conn);
  free(psi);
  free(nextNode);
  free(coords);
  free(dispf);
  free(invdetf);
}
