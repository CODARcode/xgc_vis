#include "io/xgcMesh.h"
#include "io/xgcData.h"

#include <vtkDataSet.h>
#include <vtkDataSetWriter.h>
#include <vtkContourFilter.h>
#include <vtkContourGrid.h>

int main(int argc, char **argv)
{
  XGCMesh m;
  m.readMeshFromADIOS(argv[1], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);

  ADIOS_FILE *varFP = adios_read_open_file(argv[2], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  XGCData d;
  d.readDpotFromADIOS(m, varFP);

  m.buildNeighbors();
  m.buildNodeGraph();
  m.marchingTriangles(m.psi, 0.2);

#if 1
  vtkDataSet *grid = d.convert2DSliceToVTK(m);

  // vtkContourGrid *contourGrid = vtkContourGrid::New();
  // contourGrid->SetInputData(grid);
  // contourGrid->GenerateValues(1, 0.2, 0.2);
  // contourGrid->Update();
  
  vtkContourFilter *contourFilter = vtkContourFilter::New();
  contourFilter->SetInputData(grid);
  contourFilter->GenerateValues(4, 0.2, 0.4); 

  vtkDataSetWriter *wrt = vtkDataSetWriter::New();
  wrt->SetFileTypeToBinary();
  wrt->SetFileName("myxgc.vtk");
  // wrt->SetInputData(grid);
  wrt->SetInputData(contourFilter->GetOutput());
  // wrt->SetInputData(contourGrid->GetOutput());
  wrt->Write();
  wrt->Delete();
#endif

  return 0;
}
