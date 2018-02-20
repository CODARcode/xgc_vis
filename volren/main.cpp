#include <iostream>
#include <fstream>
#include "volren/bvh.h"
#include "volren/volren.cuh"

int main(int argc, char **argv)
{
  float *framebuf = NULL; // FIXME
  ctx_rc *rc;
  rc_create_ctx(&rc);
  rc_set_stepsize(rc, 0.5);
  rc_set_viewport(rc, 0, 0, 1024, 768);

  rc_clear_output(rc);
  // rc_set_invmvpd(rc, imvmvpd);
  rc_render(rc);
  // rc_dump_output(rc, framebuf);

  return 0;
}
