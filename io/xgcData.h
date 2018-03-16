#ifndef _XGCDATA_H
#define _XGCDATA_H

#include "def.h"
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
  void deriveTotalEnergy();

  ~XGCData();

#if WITH_ADIOS
  void readDpotFromADIOS(XGCMesh &m, ADIOS_FILE *fp);
#endif

#if WITH_H5
  void readDpotFromH5(XGCMesh& m, const std::string& filename);
#endif

  using json = nlohmann::json;
  json jsonfyData(const XGCMesh&) const; 
  json jsonfySingleSliceData(const XGCMesh&) const; 
  json jsonfyDataInfo(const XGCMesh&) const;

#if WITH_VTK
  struct vtkDataSet* convert2DSliceToVTK(XGCMesh& m);
#endif

  std::vector<double> sampleAlongPsiContour(const XGCMesh &m, double isoval);
  std::vector<double> sampleAlongPsiContourPolar(const XGCMesh &m, double isoval);
};

#endif
