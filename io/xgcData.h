#ifndef _XGCDATA_H
#define _XGCDATA_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include "io/xgcMesh.h"

struct XGCData {
  double *dpot = NULL;
  float *dpotf = NULL, *graddpotf = NULL;
  
  void deriveSinglePrecisionDpot(const XGCMesh& m);
  void deriveGradient(const XGCMesh& m);

  ~XGCData();
  void readDpotFromADIOS(XGCMesh &m, ADIOS_FILE *fp);

#if WITH_VTK
  struct vtkDataSet* convert2DSliceToVTK(XGCMesh& m);
#endif
};

#endif
