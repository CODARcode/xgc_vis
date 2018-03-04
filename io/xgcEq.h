#ifndef _EQDATA_ 
#define _EQDATA_ // equilibrium data

#include <string>
#include <vector>

struct XGCEqData {
  void parseFromFile(const std::string& filename);

  int mr, mz, mpsi;
  double min_r, max_r, min_z, max_z; 
  double axis_r, axis_z, axis_b;
  double x_psi, x_r, x_z;

  std::vector<double> psi_grid;
  std::vector<double> I;
  std::vector<double> psi_rz;
};

#endif
