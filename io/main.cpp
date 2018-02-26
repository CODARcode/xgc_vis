#include "io/xgcMesh.h"
#include "io/xgcData.h"

#include <vtkDataSet.h>
#include <vtkDataSetWriter.h>

int main(int argc, char **argv)
{
  XGCMesh m;
  m.readMeshFromADIOS(argv[1], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);

  ADIOS_FILE *varFP = adios_read_open_file(argv[2], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  XGCData d;
  d.readDpotFromADIOS(m, varFP);

  vtkDataSet* grid = d.convert2DSliceToVTK(m);
  
  vtkDataSetWriter *wrt = vtkDataSetWriter::New();
  wrt->SetFileTypeToBinary();
  wrt->SetFileName("myxgc.vtk");
  wrt->SetInputData(grid);
  wrt->Write();
  wrt->Delete();

  return 0;
}
