#ifndef _XGCDATA_H
#define _XGCDATA_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <json.hpp>
#include "io/xgcMesh.h"

struct XGCData {
  double *dpot = NULL;
  float *dpotf = NULL, *graddpotf = NULL;
  float dpotf_min = FLT_MAX, dpotf_max = -FLT_MAX;
  
  void deriveSinglePrecisionDpot(const XGCMesh& m);
  void deriveGradient(const XGCMesh& m);

  ~XGCData();
  void readDpotFromADIOS(XGCMesh &m, ADIOS_FILE *fp);

  using json = nlohmann::json;
  json jsonfyData(const XGCMesh&) const; 
  json jsonfyDataInfo(const XGCMesh&) const;

#if WITH_VTK
  struct vtkDataSet* convert2DSliceToVTK(XGCMesh& m);
#endif
  
  std::vector<double> sampleAlongPsiContour(const XGCMesh &m, double isoval);
};

#endif
