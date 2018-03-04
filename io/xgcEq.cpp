#include "xgcEq.h"
#include <fstream>
#include <cassert>

void XGCEq::parseFromFile(const std::string& filename)
{
  std::ifstream ifs(filename, std::ifstream::in);
  
  // ignore the first line
  ifs.ignore(1000, '\n');

  ifs >> mr >> mz >> mpsi;
  ifs >> min_r >> max_r >> min_z >> max_z;
  ifs >> axis_r >> axis_z >> axis_b;
  ifs >> x_psi >> x_r >> x_z;
  
  psi_grid.resize(mpsi);
  for (int i=0; i<mpsi; i++) ifs >> psi_grid[i];

  I.resize(mpsi);
  for (int i=0; i<mpsi; i++) ifs >> I[i];

  psi_rz.resize(mr*mz);
  for (int i=0; i<mr*mz; i++) ifs >> psi_rz[i];

  int end_flag;
  ifs >> end_flag;
  assert(end_flag == -1);

  ifs.close();

  printInfo();
}

void XGCEq::printInfo() const 
{
  fprintf(stderr, "mr=%d, mz=%d, mpsi=%d\n", mr, mz, mpsi);
  fprintf(stderr, "min_r=%f, max_r=%f, min_z=%f, max_z=%f\n", 
      min_r, max_r, min_z, max_z);
  fprintf(stderr, "axis_r=%f, axis_z=%f, axis_b=%f\n", 
      axis_r, axis_z, axis_b);
  fprintf(stderr, "x_psi=%f, x_r=%f, x_z=%f\n", 
      x_psi, x_r, x_z);
}
