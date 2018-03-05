#ifndef _VOLREN_CUH
#define _VOLREN_CUH

#include "def.h"

#ifdef __cplusplus
extern "C" {
#endif

struct BVHNodeD;
struct cudaArray; 

static const int size_tf = 256;

struct ctx_rc {
  int viewport[4];
  float invmvp[16]; // inverse(proj*modelview)

  int *d_neighbors, *h_neighbors; // neighbor indices are for quadnodes.
  BVHNodeD *d_bvh, *h_bvh;
  int *d_viewport;
  float *d_invmvp;

  float *d_invdet, *h_invdet;
  float *d_disp, *h_disp;
  float *d_grad, *h_grad;
  float *d_psi, *h_psi;
  float psi_min, psi_max;
  
  bool toggle_psi_range;
  float psi_range_min, psi_range_max;

  bool toggle_angle_range;
  float angle_range_min, angle_range_max;

  bool toggle_shading;
  float Ka, Kd, Ks;
  float light_direction[3];

  bool toggle_slice_highlight;
  float slice_highlight_ratio;

  float *d_data, *h_data;
  int nNodes, nPhi, iPhi, nTriangles, nBVHNodes;

  // float *d_output;
  unsigned char *d_output_rgba8;
  void *h_output;

  float *h_tf, *d_tf; // fixed length 1024
  float *h_ptf, *d_ptf; // preintegrated tf
  
  float stepsize;
  float trans[2];

  int transform; 
  int shading;
};

void rc_create_ctx(ctx_rc **ctx); 
void rc_destroy_ctx(ctx_rc **ctx); 

void rc_bind_invdet(ctx_rc *ctx, int nTriangles, float *invdet); // determinant of triangles
void rc_bind_psi(ctx_rc *ctx, int nNodes, float *psi, float psi_min, float psi_max);
void rc_bind_disp(ctx_rc *ctx, int nNodes, float *disp); // displacements of nodes
void rc_bind_bvh(ctx_rc *ctx, int nBVHNodes, BVHNodeD *bvh);
void rc_bind_neighbors(ctx_rc *ctx, int nTriangles, int *neighbors);
void rc_bind_data(ctx_rc *ctx, int nNodes, int nTriangles, int nPhi, int iPhi, float *data, float *grad);

void rc_set_default_tf(ctx_rc *ctx);
void rc_set_tf(ctx_rc *ctx, float *tf); 
void rc_set_stepsize(ctx_rc *ctx, float stepsize); 
void rc_set_viewport(ctx_rc *ctx, int x, int y, int w, int h);
void rc_set_invmvpf(ctx_rc *ctx, float *invmvp); 
void rc_set_invmvpd(ctx_rc *ctx, double *invmvp); 
void rc_set_range(ctx_rc *ctx, float dpot_min, float dpot_max);
void rc_set_psi_range(ctx_rc *ctx, bool on, float psi_min, float psi_max);
void rc_set_angle_range(ctx_rc *ctx, bool on, float angle_min, float angle_max);
void rc_set_shading(ctx_rc *ctx, bool on, float Ka, float Kd, float Ks, float lx, float ly, float lz);
void rc_set_slice_highlight_ratio(ctx_rc *ctx, bool on, float ratio);

void rc_render(ctx_rc *ctx);
void rc_render_cpu(ctx_rc *ctx);

void rc_test_point_locator(ctx_rc *ctx, float x, float y);

void rc_clear_output(ctx_rc *ctx); 
// void rc_dump_output(ctx_rc *ctx, float *output); 
void rc_copy_output_to_host(ctx_rc *ctx); 

#ifdef __cplusplus
}
#endif

#endif
