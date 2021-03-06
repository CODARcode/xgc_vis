#include "io/xgcData.h"
#include "io/bp_utils.hpp"
#include "common/base64.h"
#include <iostream>
#include <string>

#if WITH_H5
#include <hdf5.h>
#endif

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

using json = nlohmann::json;

XGCData::~XGCData() {
}

void XGCData::deriveSinglePrecisionDpot(const XGCMesh& m) {
  dpotf_min = FLT_MAX; dpotf_max = -FLT_MAX;
  dpotf.resize(m.nNodes*m.nPhi);
  for (int i=0; i<m.nNodes*m.nPhi; i++) {
    dpotf[i] = dpot[i];
    dpotf_min = std::min(dpotf_min, dpotf[i]);
    dpotf_max = std::max(dpotf_max, dpotf[i]);
  }
}

void XGCData::deriveGradient(const XGCMesh& m) {
  graddpotf.resize(m.nTriangles*m.nPhi*2);
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

#if WITH_ADIOS
void XGCData::readDpotFromADIOS(XGCMesh &m, ADIOS_FILE *fp)
{
  adios_read_bp_reset_dimension_order(fp, 0);

  readValueInt(fp, "nphi", &m.nPhi);
  readValueInt(fp, "iphi", &m.iPhi);
    
  ADIOS_VARINFO *avi = adios_inq_var(fp, "dpot");
  assert(avi != NULL);

  adios_inq_var_stat(fp, avi, 0, 0);
  adios_inq_var_blockinfo(fp, avi);
  adios_inq_var_meshinfo(fp, avi);

  uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {static_cast<uint64_t>(m.nPhi), static_cast<uint64_t>(m.nNodes), 1, 1};
  ADIOS_SELECTION *sel = adios_selection_boundingbox(avi->ndim, st, sz);

  assert(sel->type == ADIOS_SELECTION_BOUNDINGBOX);

  dpot.resize(m.nPhi*m.nNodes);

  adios_schedule_read_byid(fp, sel, avi->varid, 0, 1, &dpot[0]);
  adios_perform_reads(fp, 1);
  adios_selection_delete(sel);
}
#endif

#if WITH_H5
void XGCData::readDpotFromH5(XGCMesh &m, const std::string& filename)
{
  hid_t h5fid = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  
  hid_t h5id_nphi = H5Dopen2(h5fid, "/nphi", H5P_DEFAULT);
  H5Dread(h5id_nphi, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m.nPhi);
  H5Dclose(h5id_nphi);
  
  hid_t h5id_iphi = H5Dopen2(h5fid, "/iphi", H5P_DEFAULT);
  H5Dread(h5id_iphi, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m.iPhi);
  H5Dclose(h5id_iphi);

  // dpot
  std::vector<double> dpot1(m.nPhi*m.nNodes);
  hid_t h5id_dpot = H5Dopen2(h5fid, "/dpot", H5P_DEFAULT);
  H5Dread(h5id_dpot, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dpot1[0]);
  H5Dclose(h5id_dpot);
 
  dpot.resize(m.nPhi*m.nNodes);
  for (int i=0; i<m.nPhi; i++) 
    for (int j=0; j<m.nNodes; j++)
      dpot[i*m.nNodes + j] = dpot1[j*m.nPhi + i];
 
  // narmalized dpot
  std::vector<double> dneOverne0_(m.nPhi*m.nNodes);
  hid_t h5id_dneOverne0 = H5Dopen2(h5fid, "/dneOverne0", H5P_DEFAULT);
  H5Dread(h5id_dneOverne0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dneOverne0_[0]);
  H5Dclose(h5id_dneOverne0);

#if 1
  dneOverne0.resize(m.nPhi*m.nNodes);
  for (int i=0; i<m.nPhi; i++) 
    for (int j=0; j<m.nNodes; j++)
      dneOverne0[i*m.nNodes + j] = dneOverne0_[j*m.nPhi + i];
#endif

  H5Fclose(h5fid);
}
#endif

json XGCData::jsonfyDataInfo(const XGCMesh&) const 
{
  json j;
  j["dpot_min"] = dpotf_min; 
  j["dpot_max"] = dpotf_max;

  return j;
}

json XGCData::jsonfySingleSliceData(const XGCMesh& m) const 
{
  json j = jsonfyDataInfo(m);
  j["dpot"] = base64_encode((unsigned char*)dpotf.data(), sizeof(float)*m.nNodes);

  return j;
}

json XGCData::jsonfyData(const XGCMesh& m) const 
{
  json j = jsonfyDataInfo(m);
  j["dpot"] = base64_encode((unsigned char*)dpotf.data(), sizeof(float)*m.nNodes*m.nPhi);

  return j;
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

std::vector<double> XGCData::sampleAlongPsiContourPolar(const XGCMesh &m, double isoval)
{
  std::vector<double> contour = sampleAlongPsiContour(m, isoval);
  std::vector<double> output;
  for (int i=0; i<contour.size()/3; i++) {
    double X = contour[i*3] - m.coords_centroid_x, Y = contour[i*3+1] - m.coords_centroid_y, val = contour[i*3+2];
    output.push_back(atan2(Y, X));
    output.push_back(val);
  }
  return output;
}

std::vector<double> XGCData::sampleAlongPsiContour(const XGCMesh &m, double isoval)
{
  struct float3 {
    double x, y, z;
    double& operator[](int i) {if (i==0) return x; else if (i==1) return y; else return z;}
  };

  std::map<int, float3> candidateTriangles;

  auto findZero = [&m, isoval](int i0, int i1, double &alpha) {
    double f0 = m.psi[i0], f1 = m.psi[i1];
    alpha = (isoval - f0) / (f1 - f0);
    bool b = alpha >= 0 && alpha < 1;
    if (!b) alpha = std::nan("");
    return b;
  };

  auto findIntersection = [&findZero, &candidateTriangles, &m, isoval](int triangleId, const int i[3]) {
    float3 a;
    bool b[3] = {findZero(i[0], i[1], a[0]), findZero(i[1], i[2], a[1]), findZero(i[2], i[0], a[2])};
    if (b[0] || b[1] || b[2])
      candidateTriangles[triangleId] = a;
  };

  for (int i=0; i<m.nTriangles; i++)
    findIntersection(i, &m.conn[i*3]);

  // traversal
  auto traceContourTrianglesOnSingleDirection = [this, &m, &candidateTriangles](std::list<double>& contour, int seed, bool forward) {
    int current = seed;
    while (candidateTriangles.find(current) != candidateTriangles.end()) {
      float3 intersections = candidateTriangles[current];
      candidateTriangles.erase(current);
  
      int edge;
      double alpha;
      for (int j=0; j<3; j++) {
        if ((!std::isnan(intersections[j])) && candidateTriangles.find(m.neighbors[current*3+j]) != candidateTriangles.end()) {
          edge = j; 
          alpha = intersections[j];
          break;
        }
      }
      
      int n0 = m.conn[current*3+edge], n1 = m.conn[current*3+(edge+1)%3];
      double X = (1-alpha) * m.coords[n0*2] + alpha * m.coords[n1*2], 
             Y = (1-alpha) * m.coords[n0*2+1] + alpha * m.coords[n1*2+1], 
             val = (1-alpha) * dpot[n0] + alpha * dpot[n1];
      if (forward) {contour.push_back(X); contour.push_back(Y); contour.push_back(val);}
      else {contour.push_back(val); contour.push_front(Y); contour.push_front(X);}

      current = m.neighbors[current*3+edge];
    }
  };
    
  int seed = candidateTriangles.begin()->first;
  std::list<double> contour;
  traceContourTrianglesOnSingleDirection(contour, seed, true);
  return std::vector<double>(std::begin(contour), std::end(contour));
}
