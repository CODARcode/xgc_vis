#include "core/xgcLevelSetAnalysis.h"
#include <ftk/algorithms/cca.h>
#include <set>
#include <queue>

std::vector<std::set<size_t> > XGCLevelSetAnalysis::extractSuperLevelSetOfEnergy2D
  (const XGCMesh&m, const XGCData &d, double percent)
{
  // compute total energy 
  double total_energy = 0;
  for (int i=0; i<m.nNodes/**m.nPhi*/; i++) 
    total_energy += d.dpot[i] * d.dpot[i];
  fprintf(stderr, "total_energy=%f\n", total_energy);

  // compute threshold
  double threshold_energy = total_energy * percent;
  fprintf(stderr, "threshold_energy=%f\n", threshold_energy);

  // sorting 
  std::vector<double> total_order(m.nNodes/**m.nPhi*/, 0);
  for (int i=0; i<total_order.size(); i++) total_order[i] = i;
  std::stable_sort(total_order.begin(), total_order.end(), 
      [&d](size_t v0, size_t v1) { return d.dpot[v0] < d.dpot[v1]; });

  // summation
  std::set<size_t> qualified;
  double energy = 0;
  for (int i=0; i<total_order.size(); i++) {
    size_t j = total_order[i];
    energy += d.dpot[j] * d.dpot[j];
    if (energy >= threshold_energy) 
      break;
    else qualified.insert(j);
  }

  // free total_order
  total_order.clear();
  
  auto neighbors = std::bind(&XGCMesh::getNodeNeighbors2D, &m, std::placeholders::_1);
  return ftk::extractConnectedComponents<size_t>(neighbors, qualified);
}


std::vector<std::set<size_t> > XGCLevelSetAnalysis::extractSuperLevelSet2D(const XGCMesh &m, const XGCData &d, double isoval)
{
  auto neighbors = std::bind(&XGCMesh::getNodeNeighbors2D, &m, std::placeholders::_1);
  auto criterion = [&d, isoval](size_t i) {return d.dpot[i] >= isoval;};

  return ftk::extractConnectedComponents<size_t>(m.nNodes, neighbors, criterion);
}

std::vector<std::set<size_t> > XGCLevelSetAnalysis::extractSuperLevelSet3D(const XGCMesh &m, const XGCData &d, double isoval)
{
  auto neighbors = std::bind(&XGCMesh::getNodeNeighbors3D, &m, std::placeholders::_1);
  auto criterion = [&d, isoval](size_t i) {return d.dpot[i] >= isoval;};

  return ftk::extractConnectedComponents<size_t>(m.nNodes*m.nPhi, neighbors, criterion);
}
