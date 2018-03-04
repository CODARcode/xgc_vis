#include "io/xgcMesh.h"
#include "io/xgcData.h"
#include "io/xgcEq.h"

#if WITH_VTK
#include <vtkDataSet.h>
#include <vtkDataSetWriter.h>
#include <vtkContourFilter.h>
#include <vtkContourGrid.h>
#endif

int main(int argc, char **argv)
{
  XGCEq eq;
  eq.parseFromFile(argv[1]);

#if 0
  XGCMesh m;
  m.readMeshFromADIOS(argv[1], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);

  ADIOS_FILE *varFP = adios_read_open_file(argv[2], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  XGCData d;
  d.readDpotFromADIOS(m, varFP);
#endif
  // fprintf(stderr, "psi_min_node=%d, coords={%f, %f}\n", m.psi_min_node, m.psi_min_x, m.psi_min_y);


#if 0
  m.buildNeighbors();
  m.buildNodeGraph();
  // m.marchingTriangles(m.psi, 0.2);
  // m.marchingTriangles(d.dpot, 50);
  m.sampleScalarsAlongPsiContour(d.dpot, 360, 0.2);
#endif

#if WITH_VTK
  vtkDataSet *grid = d.convert2DSliceToVTK(m);

  // vtkContourGrid *contourGrid = vtkContourGrid::New();
  // contourGrid->SetInputData(grid);
  // contourGrid->GenerateValues(1, 0.2, 0.2);
  // contourGrid->Update();
  
  vtkContourFilter *contourFilter = vtkContourFilter::New();
  contourFilter->SetInputData(grid);
  contourFilter->GenerateValues(1, 0.2, 0.2); 
  contourFilter->Update();

  auto *contours = contourFilter->GetOutput(); // vtkPolyData

  for (int i=0; i<contours->GetNumberOfCells(); i++) {
    auto *line = contours->GetCell(i);
    int id0 = line->GetPointId(0), 
        id1 = line->GetPointId(1);

    std::cout << *line;
  }
#endif

#if 0
  vtkDataSetWriter *wrt = vtkDataSetWriter::New();
  // wrt->SetFileTypeToBinary();
  wrt->SetFileName("myxgc.vtk");
  // wrt->SetInputData(grid);
  wrt->SetInputData(contourFilter->GetOutput());
  // wrt->SetInputData(contourGrid->GetOutput());
  wrt->Write();
  wrt->Delete();
#endif

  return 0;
}
