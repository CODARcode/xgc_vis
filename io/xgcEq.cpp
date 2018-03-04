#include "xgcEq.h"
#include <fstream>

void XGCEqData::parseFromFile(const std::string& filename)
{
  std::ifstream ifs(filename, std::ifstream::in);
  
  // ignore the first line
  ifs.ignore(1000, '\n');

  ifs >> mr >> mz >> mpsi;

  fprintf(stderr, "mr=%d, mz=%d, mpsi=%d\n", mr, mz, mpsi);

  ifs.close();
}
