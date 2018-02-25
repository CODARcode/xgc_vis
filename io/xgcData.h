#ifndef _XGCDATA_H
#define _XGCDATA_H

#include <string>
#include "io/xgcMesh.h"

struct XGCData {
  double *dpot = NULL;
  float *dpotf = NULL, *graddpotf = NULL;

  void deriveSinglePrecisionDpot(const XGCMesh& m);
  void deriveGradient(const XGCMesh& m);

  ~XGCData();
};

#endif
