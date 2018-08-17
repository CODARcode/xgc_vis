#include "core/xgcLevelSetAnalysis.h"
#include <set>
#include <queue>

template <class IdType>
std::vector<std::set<IdType> > connectedComponentAnalysis(
    const std::function<std::set<IdType>(size_t) >& neighbors,
    std::set<IdType> seeds)
{
  // extract connected components
  std::vector<std::set<IdType> > components;

  std::set<IdType> Q;

  while (!seeds.empty()) {
    Q.insert(*seeds.begin());

    std::set<IdType> visited;
    while (!Q.empty()) {
      IdType current = *Q.begin();
      Q.erase(current);
      visited.insert(current);

      for (auto neighbor : neighbors(current)) {
        if (seeds.find(neighbor) != seeds.end()
            && visited.find(neighbor) == visited.end()
            && Q.find(neighbor) == Q.end()) {
          Q.insert(neighbor);
        }
      }
    }

    for (auto v : visited)
      seeds.erase(v);

    components.push_back(visited);
  }

  return components;
}

template <class IdType>
std::vector<std::set<IdType> > connectedComponentAnalysis(
    IdType nNodes,
    const std::function<std::set<IdType>(size_t) >& neighbors,
    const std::function<bool(IdType)>& criterion)
{
  // find seeds
  std::set<IdType> seeds;
  for (IdType i=0; i<nNodes; i++) 
    if (criterion(i)) 
      seeds.insert(i);

  return connectedComponentAnalysis(neighbors, seeds);
}


#if 0
std::vector<std::set<size_t> > XGCLevelSetAnalysis::thresholdingByPercentageOfTotalEnergy
  (const XGCMesh &m, const XGCData &d, double percent)
{
  // compute total energy 
  double total_energy = 0;
  for (int i=0; i<m.nNodes*m.nPhi; i++) 
    total_energy += d.dpot[i] * d.dpot[i];
  fprintf(stderr, "total_energy=%f\n", total_energy);

  // compute threshold
  double threshold_energy = total_energy * percent;
  fprintf(stderr, "threshold_energy=%f\n", threshold_energy);

  // sorting 
  std::vector<double> total_order(m.nNodes*m.nPhi, 0);
  for (int i=0; i<total_order.size(); i++) total_order[i] = i;
  std::stable_sort(total_order.begin(), total_order.end(), 
      [&d](size_t v0, size_t v1) { return d.dpot[v0] < d.dpot[v1]; });

  // summation
  std::set<size_t> seeds;
  double energy = 0;
  for (int i=0; i<total_order.size(); i++) {
    size_t j = total_order[i];
    energy += d.dpot[j] * d.dpot[j];
    if (energy >= threshold_energy) { 
      fprintf(stderr, "energy=%f, i=%d, n=%zu, nNodes=%d, nPhi=%d\n", energy, i, total_order.size(), m.nNodes, m.nPhi); 
      break;
    }
    else seeds.insert(j);
  }

  // free total_order
  total_order.clear();

  auto find_neighbors = [&m](size_t v) {
    size_t plane = v / m.nNodes; 
    size_t node = v % m.nNodes;
    const std::set<size_t> &local_neighbors = m.nodeGraph[node];

    std::set<size_t> neighbors;
    for (auto &local_neighbor : local_neighbors) {
      neighbors.insert( ((plane - 1 + m.nPhi) % m.nPhi) * m.nNodes + local_neighbor);
      neighbors.insert( plane * m.nNodes + local_neighbor );
      neighbors.insert( ((plane + 1 + m.nPhi) % m.nPhi) * m.nNodes + local_neighbor);
    }

    return neighbors;
  };

  // connected components analysis
  int current_component_id = 0;
  std::map<int, std::vector<size_t> > components;

  std::queue<size_t> Q;
  
  while (!candidates.empty()) {
    Q.push(*candidates.begin());
    while (!Q.empty()) {
      size_t current = Q.front();
      Q.pop();
      candidates.erase(current);

      // fprintf(stderr, "current_component_id=%d, current=%zu, size_candidate=%zu\n", current_component_id, current, candidates.size());
      components[current_component_id].push_back(current);

      for (auto neighbor : find_neighbors(current)) { // FIXME: neighbor ids
        if (candidates.find(neighbor) != candidates.end()) { // neighbor is above threshold and connected to the current
          if (neighbor == 281) fprintf(stderr, "f, current=%zu\n", current);
          Q.push(neighbor);
        }
      }
    }

    current_component_id ++;
  }
  
  // std::vector<int> segmentation(m.nNodes*m.nPhi, 0);
}
#endif


std::vector<std::set<size_t> > XGCLevelSetAnalysis::extractSuperLevelSet2D(const XGCMesh &m, const XGCData &d, double isoval)
{
  auto neighbors = [&m](size_t i) {return m.nodeGraph[i];};
  auto criterion = [&d, isoval](size_t i) {return d.dpot[i] >= isoval;};

  fprintf(stderr, "start cca..\n");
  std::vector<std::set<size_t> > components = connectedComponentAnalysis<size_t>(m.nNodes, neighbors, criterion);
  fprintf(stderr, "#components=%lu\n", components.size());

  return components;
}
