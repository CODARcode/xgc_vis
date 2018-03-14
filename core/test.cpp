#include <mpi.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include "io/xgcMesh.h"
#include "io/xgcData.h"
#include "core/xgcLevelSetAnalysis.h"

int main(int argc, char **argv)
{
  XGCMesh m;
  m.readMeshFromADIOS(argv[1], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);

  ADIOS_FILE *varFP = adios_read_open_file(argv[2], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  XGCData d;
  d.readDpotFromADIOS(m, varFP);

  XGCLevelSetAnalysis::thresholdingByPercentageOfTotalEnergy(m, d, 0.98);

  return 0;
}
