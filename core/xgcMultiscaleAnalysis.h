#ifndef _XGCMULTISCALEANALYSIS_H
#define _XGCMULTISCALEANALYSIS_H

#include "io/xgcMesh.h"
#include "io/xgcData.h"

struct XGCMultiscaleAnalysis {
  static void extractSurface(const XGCMesh &m, const XGCData& d, double isovalue);
};

#endif
