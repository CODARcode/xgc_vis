#include "xgcData.h"
#include "io/bp_utils.hpp"
#include <iostream>
#include <string>

#if WITH_VTK
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkLongArray.h>
#include <vtkGenericCell.h>
#include <vtkDataSetWriter.h>
#include <vtkPointData.h>
#include <vtkPoints2D.h>
#endif


XGCData::~XGCData() {
  free(dpot);
  free(dpotf);
  free(graddpotf);
}

void XGCData::deriveSinglePrecisionDpot(const XGCMesh& m) {
  dpotf = (float*)realloc(dpotf, sizeof(float)*m.nNodes*m.nPhi);
  for (int i=0; i<m.nNodes*m.nPhi; i++) 
    dpotf[i] = dpot[i];
}

void XGCData::deriveGradient(const XGCMesh& m) {
  graddpotf = (float*)realloc(graddpotf, sizeof(float)*m.nTriangles*m.nPhi*2);
  for (int j=0; j<m.nPhi; j++) 
    for (int i=0; i<m.nTriangles; i++) {
      const int i0 = m.conn[i*3], i1 = m.conn[i*3+1], i2 = m.conn[i*3+2];
      float x0 = m.coords[i0*2], x1 = m.coords[i1*2], x2 = m.coords[i2*2],
            y0 = m.coords[i0*2+1], y1 = m.coords[i1*2+1], y2 = m.coords[i2*2+1];
      double f0 = dpot[j*m.nNodes+i0], f1 = dpot[j*m.nNodes+i1], f2 = dpot[j*m.nNodes+i2];
      double invdet = m.invdetf[i]; 

      graddpotf[j*m.nTriangles*2]   = ((y1-y2)*f0 + (y2-y0)*f1 + (y1-y2)*f2) * invdet;
      graddpotf[j*m.nTriangles*2+1] = ((x2-x1)*f0 + (x0-x2)*f1 + (x2-x1)*f2) * invdet;
    }
}

void XGCData::readDpotFromADIOS(XGCMesh &m, ADIOS_FILE *fp)
{
  adios_read_bp_reset_dimension_order(fp, 0);

  readValueInt(fp, "nphi", &m.nPhi);
    
  ADIOS_VARINFO *avi = adios_inq_var(fp, "dpot");
  assert(avi != NULL);

  adios_inq_var_stat(fp, avi, 0, 0);
  adios_inq_var_blockinfo(fp, avi);
  adios_inq_var_meshinfo(fp, avi);

  uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {static_cast<uint64_t>(m.nPhi), static_cast<uint64_t>(m.nNodes), 1, 1};
  ADIOS_SELECTION *sel = adios_selection_boundingbox(avi->ndim, st, sz);

  assert(sel->type == ADIOS_SELECTION_BOUNDINGBOX);

  if (dpot == NULL) 
    dpot = (double*)malloc(sizeof(double)*m.nPhi*m.nNodes);

  adios_schedule_read_byid(fp, sel, avi->varid, 0, 1, dpot);
  adios_perform_reads(fp, 1);
  adios_selection_delete(sel);
}

#if WITH_VTK
vtkDataSet* XGCData::convert2DSliceToVTK(XGCMesh &m) 
{
  vtkUnstructuredGrid *grid = vtkUnstructuredGrid::New();

  vtkPoints *pts = vtkPoints::New();
  // vtkPoints2D *pts = vtkPoints2D::New();
  pts->SetNumberOfPoints(m.nNodes);

  for (int i=0; i<m.nNodes; i++)
    pts->SetPoint(i, m.coords[i*2], m.coords[i*2+1], 0); 

  for (int i=0; i<m.nTriangles; i++) {
    vtkIdType ids[3] = {m.conn[i*3], m.conn[i*3+1], m.conn[i*3+2]};
    grid->InsertNextCell(VTK_TRIANGLE, 3, ids);
  }

  grid->SetPoints(pts);
  // pts->Delete();

#if 0
  vtkDataArray *dpotArray = vtkDoubleArray::New();
  dpotArray->SetName("dpot");
  dpotArray->SetNumberOfComponents(1);
  dpotArray->SetNumberOfTuples(m.nNodes);
  memcpy(dpotArray->GetVoidPointer(0), dpot, sizeof(double)*m.nNodes);
#endif
  vtkDataArray *psiArray = vtkDoubleArray::New();
  psiArray->SetName("psi");
  psiArray->SetNumberOfComponents(1);
  psiArray->SetNumberOfTuples(m.nNodes);
  memcpy(psiArray->GetVoidPointer(0), m.psi, sizeof(double)*m.nNodes);

  // grid->GetPointData()->AddArray(dpotArray);
  grid->GetPointData()->AddArray(psiArray);
  grid->GetPointData()->SetActiveScalars("psi");

  return grid;
}
#endif
