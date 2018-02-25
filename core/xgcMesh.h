#ifndef _XGCMESH_H
#define _XGCMESH_H

struct XGCMesh {
  int nNodes, nTriangles, nPhi;
  int *conn = NULL, *nextNode = NULL;
  double *coords = NULL;
  double *psi = NULL;
  float *psif = NULL;
  float *dispf = NULL; // displacement derived from nextNode
  float *invdetf = NULL; // inversed determinant of triangles
  float psi_min = FLT_MAX, psi_max = -FLT_MAX;

  void deriveSinglePrecisionPsi();
  void deriveInversedDeterminants();
  void deriveDisplacements();

  ~XGCMesh();
};

#endif
