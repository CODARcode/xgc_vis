#include "volren.cuh"
#include "bvh.cuh"
#include "common.cuh"

__device__ __host__
bool QuadNodeD_insideQuad(const QuadNodeD &q, float x, float y)
{
  return x >= q.Ax && x < q.Bx && y >= q.Ay && y < q.By;
}

__device__ __host__
bool QuadNodeD_insideTriangle(const QuadNodeD &q, float x, float y, float &alpha, float &beta, float &gamma) 
{
  alpha = ((q.y1 - q.y2)*(x - q.x2) + (q.x2 - q.x1)*(y - q.y2)) /
          ((q.y1 - q.y2)*(q.x0 - q.x2) + (q.x2 - q.x1)*(q.y0 - q.y2));
  beta = ((q.y2 - q.y0)*(x - q.x2) + (q.x0 - q.x2)*(y - q.y2)) /
         ((q.y1 - q.y2)*(q.x0 - q.x2) + (q.x2 - q.x1)*(q.y0 - q.y2));
  gamma = 1.0 - alpha - beta;
  // fprintf(stderr, "barycentric=%f, %f, %f\n", alpha, beta, gamma);
  return alpha >= 0 && beta >= 0 && gamma >= 0;
}

__device__ __host__
int QuadNodeD_locatePoint_recursive(const QuadNodeD *q, const QuadNodeD *nodes, float x, float y, float &alpha, float &beta, float &gamma)
{
  if (q->triangleId >= 0) { //leaf node
    bool succ = QuadNodeD_insideTriangle(*q, x, y, alpha, beta, gamma);
    if (succ) return q->triangleId;
  } else if (QuadNodeD_insideQuad(*q, x, y)) {
    for (int j=0; j<4; j++) {
      if (q->childrenIds[j] > 0) {
        int result = QuadNodeD_locatePoint_recursive(&nodes[q->childrenIds[j]], nodes, x, y, alpha, beta, gamma);
        if (result >= 0) return result;
      }
    }
  }
  return -1;
}

__device__ __host__
int QuadNodeD_locatePoint(QuadNodeD *nodes, float x, float y, float &alpha, float &beta, float &gamma)
{
  // float alpha, beta, gamma;
  static const int maxStackSize = 64;
  int stack[maxStackSize];
  int stackPos = 0;
  stack[stackPos++] = 0; // push root

  while (stackPos > 0) {
    const int i = stack[--stackPos]; // pop
    const QuadNodeD &q = nodes[i];

    // fprintf(stderr, "D_checking node %d, %f, %f, %f, %f\n", i, q.Ax, q.Ay, q.Bx, q.By);
    // fprintf(stderr, "D_checking node %d\n", i);

    if (q.triangleId >= 0) { // leaf node
      bool succ = QuadNodeD_insideTriangle(q, x, y, alpha, beta, gamma);
      if (succ) return q.triangleId;
    } else if (QuadNodeD_insideQuad(q, x, y)) { // non-leaf node
      for (int j=0; j<4; j++) {
        if (q.childrenIds[j] > 0)
          stack[stackPos++] = q.childrenIds[j];
      }
    }
  }
  return -1;
}

__device__ __host__
float QuadNodeD_sample(QuadNodeD *nodes, float x, float y, float *scalar) {
  float alpha, beta, gamma;
  int i = QuadNodeD_locatePoint(nodes, x, y, alpha, beta, gamma);
  const QuadNodeD &q = nodes[i];

  return alpha * scalar[q.i0] 
    + beta * scalar[q.i1]
    + gamma * scalar[q.i2];
}

texture<float4, 1, cudaReadModeElementType> texTransferFunc;

__constant__ int c_viewport[4];
__constant__ float c_invmvp[16]; 


__device__ __host__
float interpolateXGC(QuadNodeD *bvh, float3 pos, float *data)
{
  // compute cylindrical coordiates
  float r = sqrt(pos.x*pos.x + pos.y + pos.y);
  float theta = atan2(pos.y, pos.x);
  float z = pos.z;
  
  return 0;
}

template <int SHADING>
__device__ static void raycasting(
        float4 &dst,              // destination color
        int nPhi,                 // number of planes
        int nNodes,               // number of nodes 
        float *data,              // volume data in unstructured mesh
        QuadNodeD *bvh,
        float2 trans,             // range transformation 
        float3 rayO,              // ray origin 
        float3 rayD,              // ray direction
        float stepsize)           // stepsize
{
  float4 src = make_float4(0); 
  float3 pos; // actual position 
  float3 coords; // unnormalized tex coordinates 
  float  sample = 0.f; 
  float3 N, L = make_float3(1, 0, 0), V = rayD; 
  float3 Ka = make_float3(0.04), 
         Kd = make_float3(0.3), 
         Ks = make_float3(0.2); 
  // const float delta = 0.5f / dsz.x;   // for shading 

  const QuadNodeD &root = bvh[0];

  const float radius0 = root.Ax, radius1 = root.Bx, 
              z0 = root.Ay, z1 = root.By; 
  float tnear0, tfar0, tnear1, tfar1;

  bool b0 = intersectCylinder(rayO, rayD, tnear0, tfar0, radius0, z0, z1), 
       b1 = intersectCylinder(rayO, rayD, tnear1, tfar1, radius1, z0, z1);
  if (!(b0 || b1)) return;

  // tnear = ceilf((tnear - tnear0) / stepsize) * stepsize + tnear0; // align the ray with the global entry point
  float t = max(0.f, tnear0); // looking inside the volume

  while(1) { // raycasting
    pos = rayO + rayD*t;

    sample = interpolateXGC(bvh, pos, data); 

    // sample = QuadNodeD_sample(bvh, x, y, data);

    // sample = tex3Dtrans<DataType, readMode, TRANSFORM>(texVolume, trans, coords); 
    // src = tex1D(texTransferFunc, sample);
    // src = make_float4(sample, 1.0-sample, 0.0, 0.9);
    sample = pow(1.f - sample, 2.f); 
    src = make_float4(sample*2, 1.f-sample*2, 0.0, sample*0.4); 
   
#if 0
    if (SHADING) {
      float3 lit; 
      N = gradient(texVolume, coords, delta); 
      lit = cook(N, V, L, Ka, Kd, Ks); 
      src.x += lit.x; 
      src.y += lit.y; 
      src.z += lit.z; 
    }
#endif

    src.w = 1.f - pow(1.f - src.w, stepsize); // alpha correction  

    dst.x += (1.0 - dst.w) * src.x * src.w;
    dst.y += (1.0 - dst.w) * src.y * src.w;
    dst.z += (1.0 - dst.w) * src.z * src.w;
    dst.w += (1.0 - dst.w) * src.w;
    
    t += stepsize; 
    
    if (t>tfar1) break; // no early ray termination in compositing mode
    // if (t>tfar || dst.w>opacityThreshold) break;
  }
}

template <int SHADING>
__global__ static void raycasting_kernel(
        float *output, 
        float2 trans, 
        float stepsize, 
        float3 dsz, 
        float3 st, 
        float3 sz, 
        float3 gst, 
        float3 gsz)
{
  uint x = blockIdx.x*blockDim.x + threadIdx.x;
  uint y = blockIdx.y*blockDim.y + threadIdx.y;

  if (x >= c_viewport[2] || y>= c_viewport[3]) return;
  
  float coord[4], obj0[4], obj1[4]; 
  coord[0] = (x-c_viewport[0])*2.f / c_viewport[2] - 1.f; 
  coord[1] = (y-c_viewport[1])*2.f / c_viewport[3] - 1.f; 
  coord[2] = -1.0; 
  coord[3] = 1.0;

  mulmatvec(c_invmvp, coord, obj0); 
  coord[2] = 1.0; 
  mulmatvec(c_invmvp, coord, obj1); 
  if (obj0[3] == 0.f || obj1[3] == 0.f) return; 

  for (int i=0; i<3; i++)
      obj0[i] /= obj0[3], obj1[i] /= obj1[3]; 

  float3 rayO = make_float3(obj0[0], obj0[1], obj0[2]), 
         rayD = normalize(make_float3(obj1[0]-obj0[0], obj1[1]-obj0[1], obj1[2]-obj0[2]));

  float4 dst = make_float4(0.f); 

#if 0
      raycasting<unsigned char, cudaReadModeNormalizedFloat, BLOCKING, TRANSFORM, SHADING>(
                  dst, 
                  texVolumeUchar, trans, 
                  rayO, rayD, stepsize,  
                  dsz, st, sz, gst, gsz); 
#endif

  // GL_ONE_MINUS_DST_ALPHA, GL_ONE
  float w0 = 1-output[(y*c_viewport[2]+x)*4+3]; //, w1 = 1; make the compiler happy :)

  output[(y*c_viewport[2]+x)*4+0] += w0* dst.x;
  output[(y*c_viewport[2]+x)*4+1] += w0* dst.y;
  output[(y*c_viewport[2]+x)*4+2] += w0* dst.z;
  output[(y*c_viewport[2]+x)*4+3] += w0* dst.w;
}


/////////////////////////////
extern "C" {

void rc_render(ctx_rc *ctx)
{
  const dim3 blockSize(16, 16); 
  const dim3 gridSize = dim3(iDivUp(ctx->viewport[2], blockSize.x), iDivUp(ctx->viewport[3], blockSize.y));

  cudaMemcpyToSymbol(c_viewport, ctx->viewport, sizeof(int)*4);
  cudaMemcpyToSymbol(c_invmvp, ctx->invmvp, sizeof(float)*16);

#if 0
  raycasting_kernel<RCKERNEL_UCHAR, RCBLOCKING_DISABLED, RCTRANSFORM_DISABLED, RCSHADING_NONE><<<gridSize, blockSize>>>(
          ctx->d_output, 
          make_float2(ctx->trans[0], ctx->trans[1]), ctx->stepsize,  
          make_float3(ctx->dsz[0], ctx->dsz[1], ctx->dsz[2]), 
          make_float3(ctx->st[0], ctx->st[1], ctx->st[2]), 
          make_float3(ctx->sz[0], ctx->sz[1], ctx->sz[2]), 
          make_float3(ctx->gst[0], ctx->gst[1], ctx->gst[2]), 
          make_float3(ctx->gsz[0], ctx->gsz[1], ctx->gsz[2])); 
#endif

  checkLastCudaError("[rc_render]");
}

void rc_bind_transfer_function_array(cudaArray* array)
{
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>(); 

  texTransferFunc.normalized = true; 
  texTransferFunc.filterMode = cudaFilterModeLinear; 
  texTransferFunc.addressMode[0] = cudaAddressModeClamp; 
  cudaBindTextureToArray(texTransferFunc, array, channelDesc); 

  checkLastCudaError("[rc_bind_transfer_function_array]");
}

void rc_bind_data(ctx_rc *ctx, int nNode, int nPhi, const float *data)
{
  ctx->nNode = nNode;
  ctx->nPhi = nPhi;
  if (ctx->d_data == NULL)
  cudaMalloc((void**)&ctx->d_data, sizeof(float)*nNode*nPhi);
  cudaMemcpy(ctx->d_data, data, sizeof(float)*nNode*nPhi, cudaMemcpyHostToDevice);
}

void rc_create_ctx(ctx_rc **ctx)
{
  *ctx = (ctx_rc*)malloc(sizeof(ctx_rc));
  memset(*ctx, 0, sizeof(ctx_rc));

  cudaMalloc((void**)&((*ctx)->d_output), sizeof(float)*4096*4096); 
}

void rc_destroy_ctx(ctx_rc **ctx)
{
  // TODO: free any resources
  free(*ctx); 
  *ctx = NULL; 
}

void rc_set_viewport(ctx_rc *ctx, int x, int y, int w, int h)
{
  ctx->viewport[0] = x; 
  ctx->viewport[1] = y; 
  ctx->viewport[2] = w; 
  ctx->viewport[3] = h; 
}

void rc_set_range(ctx_rc *ctx, float a, float b)
{
  float c = 1.f/(b-a);
  ctx->trans[0] = c; 
  ctx->trans[1] = -a*c; 
}

void rc_set_stepsize(ctx_rc *ctx, float stepsize)
{
  ctx->stepsize = stepsize;
}

void rc_set_invmvpf(ctx_rc *ctx, float *invmvp)
{
  memcpy(ctx->invmvp, invmvp, sizeof(float)*16); 
}

void rc_set_invmvpd(ctx_rc *ctx, double *invmvp)
{
  for (int i=0; i<16; i++) {
    ctx->invmvp[i] = invmvp[i]; 
  }
}

void rc_clear_output(ctx_rc *ctx)
{
  cudaMemset(ctx->d_output, 0, 4*sizeof(float)*ctx->viewport[2]*ctx->viewport[3]);
}

void rc_dump_output(ctx_rc *ctx, float *output)
{
  cudaMemcpy(output, ctx->d_output, 4*sizeof(float)*ctx->viewport[2]*ctx->viewport[3], cudaMemcpyDeviceToHost); 
}

} // extern "C" 
/////////////////
