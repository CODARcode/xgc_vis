#ifndef _VOLREN_CUH
#define _VOLREN_CUH

#ifdef __cplusplus
extern "C" {
#endif

struct QuadNodeD;
struct cudaArray; 

struct ctx_rc {
  int viewport[4];
  float projmatrix[16], mvmatrix[16], invmvp[16]; 

  QuadNodeD *d_bvh;
  int *d_viewport;
  float *d_invmvp;

  float *d_data;
  int nNodes, nPhi;

  float *d_output;
  void *h_output;
  cudaArray *d_tf; 
  int size_tf;
  float stepsize;
  float trans[2];

  int transform; 
  int shading;
};  

void rc_create_ctx(ctx_rc **ctx); 
void rc_destroy_ctx(ctx_rc **ctx); 

void rc_bind_transfer_function_array(cudaArray* array); 
void rc_bind_bvh(ctx_rc *ctx, int nQuadNodes, QuadNodeD *nodes);
void rc_bind_data(ctx_rc *ctx, int nNodes, int nPhi, const float *data);

void rc_set_stepsize(ctx_rc *ctx, float stepsize); 
void rc_set_viewport(ctx_rc *ctx, int x, int y, int w, int h);
void rc_set_invmvpf(ctx_rc *ctx, float *invmvp); 
void rc_set_invmvpd(ctx_rc *ctx, double *invmvp); 
void rc_set_range(ctx_rc *ctx, float a, float b); 

void rc_render(ctx_rc *ctx);

void rc_clear_output(ctx_rc *ctx); 
// void rc_dump_output(ctx_rc *ctx, float *output); 
void rc_copy_output_to_host(ctx_rc *ctx); 
void rc_copy_output_to_host_rgb8(ctx_rc *ctx); 

#ifdef __cplusplus
}
#endif

#endif
