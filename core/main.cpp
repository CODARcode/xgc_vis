#include <mpi.h>
#include <adios.h>
#include <adios_read.h>
#include <cassert>
#include <iostream>
#include "core/bp_utils.hpp"
#include "core/xgcBlobExtractor.h"

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  ADIOS_FILE *meshFP = adios_read_open_file("xgc.mesh.bp", ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  ADIOS_FILE *varFP = adios_read_open_file("xgc.3d.bp", ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
    
  adios_read_bp_reset_dimension_order(meshFP, 0);
  adios_read_bp_reset_dimension_order(varFP, 0);

  int nNodes, nTriangles, nPhi;
  readValueInt(varFP, "nphi", &nPhi);

  fprintf(stderr, "reading mesh...\n");
  double *coords; 
  int *conn;
  readTriangularMesh(meshFP, nNodes, nTriangles, &coords, &conn);
  // fprintf(stderr, "nNodes=%d, nTriangles=%d, nPhi=%d\n", 
  //     nNodes, nTriangles, nPhi);

  fprintf(stderr, "reading data...\n");
  double *dpot;
  readScalars<double>(varFP, "dpot", &dpot);

  fprintf(stderr, "starting analysis..\n");
  XGCBlobExtractor *extractor = new XGCBlobExtractor;
  extractor->setMesh(nNodes, nTriangles, nPhi, coords, conn);
  extractor->setData(dpot);
  extractor->buildContourTree3D();

  fprintf(stderr, "dumping results..\n");
  // extractor->dumpInfo("xgc.info");
  extractor->dumpLabels("xgc.labels");
  extractor->dumpBranchDecompositions("xgc.branches");

  delete extractor;
  free(dpot);
  free(coords);
  free(conn);

  MPI_Finalize();
  return 0;
}
