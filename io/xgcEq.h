#ifndef _EQDATA_ 
#define _EQDATA_ // equilibrium data

#include "def.h"
#include <string>
#include <vector>

struct XGCEq {
  void parseFromFile(const std::string& filename);
  void saveToNCFile(const std::string& filename);
  void printInfo() const;

  size_t mr, mz, mpsi;
  double min_r, max_r, min_z, max_z; 
  double axis_r, axis_z, axis_b;
  double x_psi, x_r, x_z;

  std::vector<double> psigrid, rgrid, zgrid;
  std::vector<double> I;
  std::vector<double> psi_rz; // 2D array
};

#endif
