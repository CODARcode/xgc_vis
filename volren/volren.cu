#include "volren.cuh"
#include "bvh.cuh"
#include "common.cuh"

__device__ __host__
bool QuadNodeD_insideQuad(const QuadNodeD &q, float x, float y)
{
  return x >= q.Ax && x < q.Bx && y >= q.Ay && y < q.By;
}

__device__ __host__
bool QuadNodeD_insideTriangle(const QuadNodeD &q, float x, float y, float3 &lambda)
{
  lambda.x = ((q.y1 - q.y2)*(x - q.x2) + (q.x2 - q.x1)*(y - q.y2)) /
          ((q.y1 - q.y2)*(q.x0 - q.x2) + (q.x2 - q.x1)*(q.y0 - q.y2));
  lambda.y = ((q.y2 - q.y0)*(x - q.x2) + (q.x0 - q.x2)*(y - q.y2)) /
         ((q.y1 - q.y2)*(q.x0 - q.x2) + (q.x2 - q.x1)*(q.y0 - q.y2));
  lambda.z = 1.0 - lambda.x - lambda.y;
  // fprintf(stderr, "barycentric=%f, %f, %f\n", lambda.x, lambda.y, lambda.z);
  return lambda.x >= 0 && lambda.y >= 0 && lambda.z >= 0;
}

__device__ __host__
int QuadNodeD_locatePoint_recursive(const QuadNodeD *q, const QuadNodeD *nodes, float x, float y, float3 &lambda)
{
  if (q->triangleId >= 0) { //leaf node
    bool succ = QuadNodeD_insideTriangle(*q, x, y, lambda);
    if (succ) return q->triangleId;
  } else if (QuadNodeD_insideQuad(*q, x, y)) {
    for (int j=0; j<4; j++) {
      if (q->childrenIds[j] > 0) {
        int result = QuadNodeD_locatePoint_recursive(&nodes[q->childrenIds[j]], nodes, x, y, lambda);
        if (result >= 0) return result;
      }
    }
  }
  return -1;
}

__device__ __host__
int QuadNodeD_locatePoint(QuadNodeD *nodes, float x, float y, float3 &lambda)
{
  // float lambda.x, lambda.y, lambda.z;
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
      bool succ = QuadNodeD_insideTriangle(q, x, y, lambda);
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
float QuadNodeD_sample(QuadNodeD* bvh, int nid, float3 lambda, float *data) {
  const QuadNodeD &q = bvh[nid];

  return lambda.x * data[q.i0] 
    + lambda.y * data[q.i1]
    + lambda.z * data[q.i2];
}

texture<float4, 1, cudaReadModeElementType> texTransferFunc;

// __constant__ int c_viewport[4];
// __constant__ float c_invmvp[16]; 

__device__ __host__
float interpolateXGC(QuadNodeD *bvh, float3 pos, float *data)
{
  // compute cylindrical coordiates
  float r = sqrt(pos.x*pos.x + pos.y + pos.y);
  float theta = atan2(pos.y, pos.x);
  // float z = pos.z;
  
  return 0;
}

template <int SHADING>
__device__ static void rc(
        float4 &dst,              // destination color
        int nPhi,                 // number of planes
        int nNodes,               // number of nodes 
        float *data,              // volume data in unstructured mesh
        QuadNodeD *bvh,
        float2 trans,             // range transformation 
        float3 rayO,              // ray origin 
        float3 rayD,              // ray direction
        float stepsize, 
        float tnear, float tfar)
{
  const float pi = 3.141592654f;
  const float pi2 = 2*pi;
  float4 src;
  // float3 N, L = make_float3(1, 0, 0), V = rayD; 
  // float3 Ka = make_float3(0.04), 
  //        Kd = make_float3(0.3), 
  //        Ks = make_float3(0.2); 
  // const float delta = 0.5f / dsz.x;   // for shading 
  float3 pos;
  float t = tnear;

  int last_nid = -1, nid;

  while (t < tfar) {
    pos = rayO + rayD*t;

    // cylindar coordinates
    float r = pos.x*pos.x + pos.y*pos.y;
    float phi = atan2(pos.y, pos.x) + pi;
    float z = pos.z;

    float3 lambda;
#if 0 // coherent point locating
    if (last_nid > 0 && QuadNodeD_insideTriangle(bvh[last_nid], r, z, lambda)) // coherent
      nid = last_nid;
    else 
      nid = QuadNodeD_locatePoint(bvh, r, z, lambda);
#else  
    nid = QuadNodeD_locatePoint(bvh, r, z, lambda);
#endif

    if (nid != -1) {
      const float unitAngle = pi2/nPhi;

      int p0 = (int)(phi/unitAngle)%nPhi;
      int p1 = (p0+1)%nPhi;
    
      float alpha = (phi - unitAngle*p0) / unitAngle;
      float v0 = QuadNodeD_sample(bvh, nid, lambda, data + nNodes*p0); //  + nNodes*p0);
      float v1 = QuadNodeD_sample(bvh, nid, lambda, data + nNodes*p1); //  + nNodes*p1);

      float value = (1-alpha)*v0 + alpha*v1;
      
      float v = clamp(value*0.01, -0.5f, 0.5f);
      // sample = interpolateXGC(bvh, pos, data); 
      // sample = QuadNodeD_sample(bvh, x, y, data);
      // sample = tex3Dtrans<DataType, readMode, TRANSFORM>(texVolume, trans, coords); 
      // src = tex1D(texTransferFunc, sample);
      // src = make_float4(sample, 1.0-sample, 0.0, 0.9);
      // sample = pow(1.f - sample, 2.f); 
      // src = make_float4(sample*2, 1.f-sample*2, 0.0, sample*0.4); 
      // src = make_float4(lambda.x, lambda.y, lambda.z, 0.5);

      src = make_float4(v+0.5, 0.5-v, 0, min(1.f, v*v*10));
      // src = make_float4(phi/(pi*2), 1-phi/(pi*2), 0, 0.3);
      // src = make_float4(p1/8.0, 0.5, 0, 0.3);

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
    }
    
    t += stepsize; 
  }
  // dst.x = 1; dst.y = 0; dst.z = 0; dst.w = 1;
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
  const QuadNodeD &root = bvh[0];
  const float innerRadius = root.Ax, outerRadius = root.Bx;
  // const float innerRadius = 1.0, outerRadius = 1.2;
  const float z0 = root.Ay, z1 = root.By; 
  float tnear0, tfar0, tnear1, tfar1;
  bool b0 = intersectCylinder(rayO, rayD, tnear0, tfar0, outerRadius, z0, z1), 
       b1 = intersectCylinder(rayO, rayD, tnear1, tfar1, innerRadius, z0, z1);
  
#if 1
  if (b0 && (!b1))
    rc<SHADING>(dst, nPhi, nNodes, data, bvh, trans, rayO, rayD, stepsize, tnear0, tfar0);
  else if (b0 && b1) {
    rc<SHADING>(dst, nPhi, nNodes, data, bvh, trans, rayO, rayD, stepsize, tnear0, tnear1);
    rc<SHADING>(dst, nPhi, nNodes, data, bvh, trans, rayO, rayD, stepsize, tfar1, tfar0);
  }
#else
  if (b1)
    rc<SHADING>(dst, nPhi, nNodes, data, bvh, trans, rayO, rayD, stepsize, tnear1, tfar1);
#endif
}

template <int SHADING>
__global__ static void raycasting_kernel(
        float *output,
        int *viewport, 
        float *invmvp,
        int nPhi, 
        int nNode, 
        float *data, 
        QuadNodeD *bvh,
        float2 trans, 
        float stepsize)
{
  uint x = blockIdx.x*blockDim.x + threadIdx.x;
  uint y = blockIdx.y*blockDim.y + threadIdx.y;

  if (x >= viewport[2] || y>= viewport[3]) return;
  
  float coord[4], obj0[4], obj1[4]; 
  coord[0] = (x-viewport[0])*2.f / viewport[2] - 1.f; 
  coord[1] = (y-viewport[1])*2.f / viewport[3] - 1.f; 
  coord[2] = -1.0; 
  coord[3] = 1.0;

  mulmatvec(invmvp, coord, obj0); 
  coord[2] = 1.0; 
  mulmatvec(invmvp, coord, obj1); 
  if (obj0[3] == 0.f || obj1[3] == 0.f) return; 

  for (int i=0; i<3; i++)
      obj0[i] /= obj0[3], obj1[i] /= obj1[3]; 

  float3 rayO = make_float3(obj0[0], obj0[1], obj0[2]), 
         rayD = normalize(make_float3(obj1[0]-obj0[0], obj1[1]-obj0[1], obj1[2]-obj0[2]));
  float4 dst = make_float4(0.f); 

#if 1
  raycasting<SHADING>(dst, nPhi, nNode, data, bvh, trans, rayO, rayD, stepsize);
#else
  dst.x = 1; // rayD.x;
  dst.y = 0; // rayD.y;
  dst.z = 0; // rayD.z;
  dst.w = 1;
#endif

  // GL_ONE_MINUS_DST_ALPHA, GL_ONE
  float w0 = 1-output[(y*viewport[2]+x)*4+3]; //, w1 = 1; make the compiler happy :)

  output[(y*viewport[2]+x)*4+0] += w0* dst.x;
  output[(y*viewport[2]+x)*4+1] += w0* dst.y;
  output[(y*viewport[2]+x)*4+2] += w0* dst.z;
  output[(y*viewport[2]+x)*4+3] += w0* dst.w;
}


/////////////////////////////
extern "C" {

void rc_render(ctx_rc *ctx)
{
  const dim3 blockSize(16, 16); 
  const dim3 gridSize = dim3(iDivUp(ctx->viewport[2], blockSize.x), iDivUp(ctx->viewport[3], blockSize.y));

  cudaMemcpy(ctx->d_viewport, ctx->viewport, sizeof(int)*4, cudaMemcpyHostToDevice);
  // cudaMemcpyToSymbol(c_viewport, ctx->viewport, sizeof(int)*4);
  cudaMemcpy(ctx->d_invmvp, ctx->invmvp, sizeof(float)*16, cudaMemcpyHostToDevice);
  // cudaMemcpyToSymbol(c_invmvp, ctx->invmvp, sizeof(float)*16);
 
  raycasting_kernel<0><<<gridSize, blockSize>>>(
          ctx->d_output,
          ctx->d_viewport, 
          ctx->d_invmvp,
          ctx->nPhi, 
          ctx->nNodes,
          ctx->d_data, 
          ctx->d_bvh,
          make_float2(ctx->trans[0], ctx->trans[1]), 
          ctx->stepsize);
  checkLastCudaError("[rc_render]");
}

void rc_bind_bvh(ctx_rc *ctx, int nQuadNodes, QuadNodeD *bvh)
{
  if (ctx->d_bvh != NULL)
    cudaFree(ctx->d_bvh);

  cudaMalloc((void**)&ctx->d_bvh, sizeof(QuadNodeD)*nQuadNodes);
  cudaMemcpy(ctx->d_bvh, bvh, sizeof(QuadNodeD)*nQuadNodes, cudaMemcpyHostToDevice);
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

void rc_bind_data(ctx_rc *ctx, int nNodes, int nPhi, const float *data)
{
  ctx->nNodes = nNodes;
  ctx->nPhi = nPhi;
  if (ctx->d_data == NULL)
    cudaMalloc((void**)&ctx->d_data, sizeof(float)*nNodes*nPhi);
  cudaMemcpy(ctx->d_data, data, sizeof(float)*nNodes*nPhi, cudaMemcpyHostToDevice);
}

void rc_create_ctx(ctx_rc **ctx)
{
  cudaSetDevice(0);

  *ctx = (ctx_rc*)malloc(sizeof(ctx_rc));
  memset(*ctx, 0, sizeof(ctx_rc));

  const size_t max_npx = 4096*4096;

  cudaMalloc((void**)&((*ctx)->d_output), sizeof(float)*max_npx); 
  (*ctx)->h_output = malloc(sizeof(float)*max_npx);

  cudaMalloc((void**)&((*ctx)->d_viewport), sizeof(int)*4);
  cudaMalloc((void**)&((*ctx)->d_invmvp), sizeof(float)*16);
  
  checkLastCudaError("[rc_init]");
}

void rc_destroy_ctx(ctx_rc **ctx)
{
  // TODO: free any resources
  // cudaFree(ctx->d_output);
  free((*ctx)->h_output);
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

void rc_copy_output_to_host(ctx_rc *ctx)
{
  cudaMemcpy(ctx->h_output, ctx->d_output, 4*sizeof(float)*ctx->viewport[2]*ctx->viewport[3], cudaMemcpyDeviceToHost); 
}

void rc_copy_output_to_host_rgb8(ctx_rc *ctx)
{
  const size_t npx = ctx->viewport[2] * ctx->viewport[3];
  cudaMemcpy(ctx->h_output, ctx->d_output, 4*sizeof(float)*npx, cudaMemcpyDeviceToHost);
  float *ffb = (float*)ctx->h_output;
  unsigned char *ufb = (unsigned char*)ctx->h_output;

  for (int i=0; i<npx; i++) {
    ufb[i*3] = ffb[i*4]*255;
    ufb[i*3+1] = ffb[i*4+1]*255;
    ufb[i*3+2] = ffb[i*4+2]*255;
  }
}

} // extern "C" 
/////////////////
