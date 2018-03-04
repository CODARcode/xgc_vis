#include <iostream>
#include <fstream>
#include "io/xgcMesh.h"
#include "io/xgcData.h"
#include "volren/bvh.h"
#include "volren/volren.cuh"
#include "volren/kdbvh.h"

int main(int argc, char **argv)
{
  XGCMesh m;
  m.readMeshFromADIOS(argv[1], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);

  ADIOS_FILE *varFP = adios_read_open_file(argv[2], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  XGCData d;
  d.readDpotFromADIOS(m, varFP);

  buildKDBVH(m);

#if 0
  float *framebuf = NULL; // FIXME
  ctx_rc *rc;
  rc_create_ctx(&rc);
  rc_set_stepsize(rc, 0.5);
  rc_set_viewport(rc, 0, 0, 1024, 768);

  rc_clear_output(rc);
  // rc_set_invmvpd(rc, imvmvpd);
  rc_render(rc);
  // rc_dump_output(rc, framebuf);
#endif
  return 0;
}
