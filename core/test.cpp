#include <mpi.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include "core/xgcBlobExtractor.h"

int main(int argc, char **argv)
{
  const int nNodes = 10, nTriangles = 10;
  double *coords = (double*)malloc(sizeof(double)*nNodes); // doesn't matter
  int conn[] = {
    3, 4, 9, 
    3, 1, 4, 
    1, 2, 4, 
    2, 8, 4, 
    2, 0, 8, 
    2, 7, 0, 
    2, 5, 7, 
    2, 1, 5, 
    1, 3, 5, 
    3, 6, 5};
  double psi[nNodes];
  for (int i=0; i<nNodes; i++) 
    psi[i] = i;

  XGCBlobExtractor *ex = new XGCBlobExtractor(nNodes, nTriangles, coords, conn);
  ex->setData(0, 1, psi);

  ex->buildContourTree2D(0);
  return 0;
}
