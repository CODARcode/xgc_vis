#include "xgcData.h"
#include <iostream>
#include <string>

XGCData::~XGCData() {
  free(dpot);
  free(dpotf);
  free(graddpotf);
}

void XGCData::deriveSinglePrecisionDpot(const XGCMesh& m) {
  dpotf = (float*)realloc(dpotf, sizeof(float)*m.nNodes*m.nPhi);
  for (int i=0; i<m.nNodes*m.nPhi; i++) 
    dpotf[i] = dpot[i];
}

void XGCData::deriveGradient(const XGCMesh& m) {
  graddpotf = (float*)realloc(graddpotf, sizeof(float)*m.nTriangles*m.nPhi*2);
  for (int j=0; j<m.nPhi; j++) 
    for (int i=0; i<m.nTriangles; i++) {
      const int i0 = m.conn[i*3], i1 = m.conn[i*3+1], i2 = m.conn[i*3+2];
      float x0 = m.coords[i0*2], x1 = m.coords[i1*2], x2 = m.coords[i2*2],
            y0 = m.coords[i0*2+1], y1 = m.coords[i1*2+1], y2 = m.coords[i2*2+1];
      double f0 = dpot[j*m.nNodes+i0], f1 = dpot[j*m.nNodes+i1], f2 = dpot[j*m.nNodes+i2];
      double invdet = m.invdetf[i]; 

      graddpotf[j*m.nTriangles*2]   = ((y1-y2)*f0 + (y2-y0)*f1 + (y1-y2)*f2) * invdet;
      graddpotf[j*m.nTriangles*2+1] = ((x2-x1)*f0 + (x0-x2)*f1 + (x2-x1)*f2) * invdet;
    }
}
