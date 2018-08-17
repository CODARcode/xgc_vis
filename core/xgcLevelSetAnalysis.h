#ifndef _XGC_LEVEL_SET_ANALYSIS
#define _XGC_LEVEL_SET_ANALYSIS

#include "io/xgcMesh.h"
#include "io/xgcData.h"

struct XGCLevelSetAnalysis {
  static std::vector<std::set<size_t> > thresholdingByPercentageOfTotalEnergy(const XGCMesh &m, const XGCData &d, double percent);

  static std::vector<std::set<size_t> > extractSuperLevelSet2D(const XGCMesh &m, const XGCData&d, double isoval);
};

#endif
