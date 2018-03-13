#ifndef _XGC_LEVEL_SET_ANALYSIS
#define _XGC_LEVEL_SET_ANALYSIS

#include "io/xgcMesh.h"
#include "io/xgcData.h"

struct XGCLevelSetAnalysis {
  static void thresholdByPercentageOfTotalEnergy(const XGCMesh &m, const XGCData &d, double percent);
};

#endif
