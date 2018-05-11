#include "core/xgcLevelSetAnalysis.h"
#include <set>
#include <queue>

void XGCLevelSetAnalysis::thresholdingByPercentageOfTotalEnergy(const XGCMesh &m, const XGCData &d, double percent)
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
  std::set<size_t> candidates;
  double energy = 0;
  for (int i=0; i<total_order.size(); i++) {
    size_t j = total_order[i];
    energy += d.dpot[j] * d.dpot[j];
    if (energy >= threshold_energy) { 
      fprintf(stderr, "energy=%f, i=%d, n=%zu, nNodes=%d, nPhi=%d\n", energy, i, total_order.size(), m.nNodes, m.nPhi); 
      break;
    }
    else candidates.insert(j);
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

std::map<int, std::set<size_t> > XGCLevelSetAnalysis::extractSuperLevelSet2D(const XGCMesh &m, const XGCData &d, double isoval)
{
  // find candidates
  std::set<size_t> candidates;
  for (int i=0; i<m.nNodes; i++) 
    if (d.dpot[i] >= isoval) candidates.insert(i);
  // fprintf(stderr, "#candidates=%lu\n", candidates.size());

  // connected components analysis
  int current_component_id = 0;
  std::map<int, std::set<size_t> > components;

  std::set<size_t> Q;
  // std::queue<size_t> Q;
  
  while (!candidates.empty()) {
    Q.insert(*candidates.begin());

    std::set<size_t> visited;
    while (!Q.empty()) {
      size_t current = *Q.begin();
      Q.erase(current);
      visited.insert(current);

      // fprintf(stderr, "current_component_id=%d, current=%zu, size_visited=%zu\n", current_component_id, current, visited.size());

      for (auto neighbor : m.nodeGraph[current]) { // neighbor ids
        if (candidates.find(neighbor) != candidates.end() && visited.find(neighbor) == visited.end() && Q.find(neighbor) == Q.end()) { 
          // fprintf(stderr, "pusing %lu\n", neighbor);
          Q.insert(neighbor);
        }
      }
    }
    
    for (auto v : visited) 
      candidates.erase(v);

    components[current_component_id++] = visited;
    // fprintf(stderr, "component_id=%d, size=%lu\n", current_component_id, visited.size());
  }

  return components;
  // fprintf(stderr, "#components=%lu\n", components.size());
}
