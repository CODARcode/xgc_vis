#include "volren.cuh"
#include "bvh.cuh"
#include "common.cuh"
  
#if WITH_CUDA
texture<float4, 1, cudaReadModeElementType> texTransferFunc;
#endif


__device__ __host__
inline bool QuadNodeD_insideQuad(const QuadNodeD &q, float x, float y)
{
  return x >= q.Ax && x < q.Bx && y >= q.Ay && y < q.By;
}

__device__ __host__
inline bool QuadNodeD_insideTriangle(const QuadNodeD &q, float x, float y, float3 &lambda, float *invdet)
{
#if 0
  lambda.x = ((q.y1 - q.y2)*(x - q.x2) + (q.x2 - q.x1)*(y - q.y2)) /
          ((q.y1 - q.y2)*(q.x0 - q.x2) + (q.x2 - q.x1)*(q.y0 - q.y2));
  lambda.y = ((q.y2 - q.y0)*(x - q.x2) + (q.x0 - q.x2)*(y - q.y2)) /
         ((q.y1 - q.y2)*(q.x0 - q.x2) + (q.x2 - q.x1)*(q.y0 - q.y2));
#endif
  const float d = invdet[q.triangleId];
  lambda.x = ((q.y1 - q.y2)*(x - q.x2) + (q.x2 - q.x1)*(y - q.y2)) * d; 
  lambda.y = ((q.y2 - q.y0)*(x - q.x2) + (q.x0 - q.x2)*(y - q.y2)) * d;
  lambda.z = 1.0 - lambda.x - lambda.y;
  // fprintf(stderr, "barycentric=%f, %f, %f\n", lambda.x, lambda.y, lambda.z);
  return lambda.x >= 0 && lambda.y >= 0 && lambda.z >= 0;
}

__device__ __host__
inline int QuadNodeD_locatePoint_recursive(const QuadNodeD *q, const QuadNodeD *nodes, float x, float y, float3 &lambda, float *invdet)
{
  if (q->triangleId >= 0) { //leaf node
    bool succ = QuadNodeD_insideTriangle(*q, x, y, lambda, invdet);
    if (succ) return q->triangleId;
  } else if (QuadNodeD_insideQuad(*q, x, y)) {
    for (int j=0; j<4; j++) {
      if (q->childrenIds[j] > 0) {
        int result = QuadNodeD_locatePoint_recursive(&nodes[q->childrenIds[j]], nodes, x, y, lambda, invdet);
        if (result >= 0) return result;
      }
    }
  }
  return -1;
}

__device__ __host__
inline int QuadNodeD_locatePoint(QuadNodeD *nodes, float x, float y, float3 &lambda, float *invdet)
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
      bool succ = QuadNodeD_insideTriangle(q, x, y, lambda, invdet);
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
inline int QuadNodeD_locatePoint_coherent(QuadNodeD *bvh, int last_nid, float x, float y, float3 &lambda, float *invdet)
{
  // check if last_nid is valid
  if (last_nid<0) return QuadNodeD_locatePoint(bvh, x, y, lambda, invdet);

  // check if in the same triangle
  if (QuadNodeD_insideTriangle(bvh[last_nid], x, y, lambda, invdet)) return last_nid;

  // TODO: check if in the neighbors of last_nid

  return QuadNodeD_locatePoint(bvh, x, y, lambda, invdet);
}

__device__ __host__
inline float QuadNodeD_sample(int i0, int i1, int i2, float3 lambda, float *data) {
  return lambda.x * data[i0] + lambda.y * data[i1] + lambda.z * data[i2];
}

__device__ __host__
inline float2 QuadNodeD_sample2(int i0, int i1, int i2, float3 lambda, float *data) {
  return make_float2(lambda.x * data[i0*2] + lambda.y * data[i1*2] + lambda.z * data[i2*2],
      lambda.x * data[i0*2+1] + lambda.y * data[i1*2+1] + lambda.z * data[i2*2+1]);
}

__device__ __host__
inline float QuadNodeD_sample(QuadNodeD* bvh, int nid, float3 lambda, float *data) {
  const QuadNodeD &q = bvh[nid];
  return lambda.x * data[q.i0] + lambda.y * data[q.i1] + lambda.z * data[q.i2];
}

__device__ __host__
inline float2 QuadNodeD_sample2(QuadNodeD* bvh, int nid, float3 lambda, float *data) {
  const QuadNodeD &q = bvh[nid];
  return QuadNodeD_sample2(q.i0, q.i1, q.i2, lambda, data);
  // return make_float2(lambda.x * data[q.i0*2] + lambda.y * data[q.i1*2] + lambda.z * data[q.i2*2],
  //     lambda.x * data[q.i0*2+1] + lambda.y * data[q.i1*2+1] + lambda.z * data[q.i2*2+1]);
}

__device__ __host__
inline int interpolateXGC(float &value, QuadNodeD *bvh, float3 p, int nPhi, int nNodes, float *data, float *invdet)
{
  static const float pi = 3.141592654f;
  static const float pi2 = 2*pi;
  
  // cylindar coordinates
  float r = sqrt(p.x*p.x + p.y*p.y);
  float phi = atan2(p.y, p.x) + pi;
  float z = p.z;
  float3 lambda;
  
  int nid = QuadNodeD_locatePoint(bvh, r, z, lambda, invdet);
  if (nid == -1) return nid; 
  
  const QuadNodeD &q = bvh[nid];

  const float deltaAngle = pi2/nPhi;
  int p0 = (int)(phi/deltaAngle)%nPhi;
  int p1 = (p0+1)%nPhi;

  float alpha = (phi - deltaAngle*p0) / deltaAngle;
  // float v0 = QuadNodeD_sample(bvh, nid, lambda, data + nNodes*p0);
  // float v1 = QuadNodeD_sample(bvh, nid, lambda, data + nNodes*p1);
  float v0 = QuadNodeD_sample(q.i0, q.i1, q.i2, lambda, data + nNodes*p0); //  + nNodes*p0);
  float v1 = QuadNodeD_sample(q.i0, q.i1, q.i2, lambda, data + nNodes*p1); //  + nNodes*p1);
  // float dr0 = QuadNodeD

  value = (1-alpha)*v0 + alpha*v1;
  return nid;
}

__device__ __host__
inline int interpolateXGC2(float &value, QuadNodeD *bvh, float3 p, int nPhi, int nNodes, float *data, float *disp, float *invdet)
{
  static const float pi = 3.141592654f;
  static const float pi2 = 2*pi;
  
  // cylindar coordinates
  float r = sqrt(p.x*p.x + p.y*p.y);
  float phi = atan2(p.y, p.x) + pi;
  float z = p.z;
  float3 lambda;
  
  int nid = QuadNodeD_locatePoint(bvh, r, z, lambda, invdet);
  if (nid == -1) return nid; 
      
  const float deltaAngle = pi2/nPhi;

  int p0 = (int)(phi/deltaAngle)%nPhi;
  int p1 = (p0+1)%nPhi;
  float alpha = (phi - deltaAngle*p0) / deltaAngle;

  // interpolate disp
  const QuadNodeD &q = bvh[nid];
  float dx = lambda.x * disp[q.i0*2] + lambda.y * disp[q.i1*2] + lambda.z * disp[q.i2*2];
  float dy = lambda.x * disp[q.i0*2+1] + lambda.y * disp[q.i1*2+1] + lambda.z * disp[q.i2*2+1];
 
  float3 lambda0, lambda1;
  int nid0 = QuadNodeD_locatePoint_coherent(bvh, nid, r+dx*(1-alpha), z+dy*(1-alpha), lambda0, invdet);
  int nid1 = QuadNodeD_locatePoint_coherent(bvh, nid, r+dx*alpha, z+dy*alpha, lambda1, invdet);
  if (nid0 == -1 || nid1 == -1) {
    // fprintf(stderr, "nid=%d, nid0=%d, nid1=%d, dx=%f, dy=%f\n", nid, nid0, nid1, dx, dy);
    return -1;
  }

  float v0 = QuadNodeD_sample(bvh, nid0, lambda0, data + nNodes*p0); //  + nNodes*p0);
  float v1 = QuadNodeD_sample(bvh, nid1, lambda1, data + nNodes*p1); //  + nNodes*p1);

  // if (alpha<0 || alpha>=1) fprintf(stderr, "%f\n", alpha);
  
  value = (1-alpha)*v0 + alpha*v1;
  return nid;
}

__device__ __host__ 
static inline float4 value2color(float value, float *tf, float2 trans)
{
  const float x = clamp(value * trans.x + trans.y, 0.f, 1.f);
#if 0
  float v = x-0.5;
  return make_float4(x, 1-x, 0, fminf(0.999f, v*v*40));
#else 
  static const int n = 256;
  static const float delta = 1.f / (n-1);
  const int i = min((int)(x*(n-1)), n-2) , j = i + 1;
  const float beta = x - i*delta, alpha = 1 - beta;

  return make_float4(
      alpha * tf[i*4] + beta * tf[j*4], 
      alpha * tf[i*4+1] + beta * tf[j*4+1], 
      alpha * tf[i*4+2] + beta * tf[j*4+2], 
      alpha * tf[i*4+3] + beta * tf[j*4+3]);
#endif
}

  
template <int SHADING>
__device__ __host__ static inline void rc(
        float4 &dst,              // destination color
        int nPhi,                 // number of planes
        int nNodes,               // number of nodes 
        float *data,              // volume data in unstructured mesh
        float *grad,              // gradient
        QuadNodeD *bvh,
        float *disp,
        float *invdet,
        float *tf,
        float2 trans,             // range transformation 
        float3 rayO,              // ray origin 
        float3 rayD,              // ray direction
        float stepsize, 
        float tnear, float tfar)
{
  float4 src;
  // float3 N, L = make_float3(1, 0, 0), V = rayD; 
  // float3 Ka = make_float3(0.04), 
  //        Kd = make_float3(0.3), 
  //        Ks = make_float3(0.2); 
  // const float delta = 0.5f / dsz.x;   // for shading 
  float3 pos;
  float value;
  float t = tnear;

  while (t < tfar) {
    pos = rayO + rayD*t;

    const int nid = interpolateXGC(value, bvh, pos, nPhi, nNodes, data, invdet);
    // const int nid = interpolateXGC2(value, bvh, pos, nPhi, nNodes, data, disp, invdet);
    if (nid >= 0) {
      src = value2color(value, tf, trans);

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
      
      src.w = 1.f - pow(1.f - src.w, stepsize*4); // alpha correction  

      dst.x += (1.0 - dst.w) * src.x * src.w;
      dst.y += (1.0 - dst.w) * src.y * src.w;
      dst.z += (1.0 - dst.w) * src.z * src.w;
      dst.w += (1.0 - dst.w) * src.w;
    }
    
    // if (dst.w > 0.98) return; // early ray termination
    t += stepsize; 
  }
  // dst.x = 1; dst.y = 0; dst.z = 0; dst.w = 1;
}

template <int SHADING>
__device__ __host__ static inline void raycasting(
        float4 &dst,              // destination color
        int nPhi,                 // number of planes
        int nNodes,               // number of nodes 
        float *data,              // volume data in unstructured mesh
        float *grad,
        QuadNodeD *bvh,
        float *disp,
        float *invdet,
        float *tf, 
        float2 trans,             // range transformation 
        float3 rayO,              // ray origin 
        float3 rayD,              // ray direction
        float stepsize)           // stepsize
{
  const QuadNodeD &root = bvh[0];
  const float innerRadius = root.Ax, outerRadius = root.Bx;
  const float z0 = root.Ay, z1 = root.By; 
  float tnear0=-FLT_MAX, tfar0=FLT_MAX, tnear1=-FLT_MAX, tfar1=FLT_MAX;
  bool b0 = intersectCylinder(rayO, rayD, tnear0, tfar0, outerRadius, z0, z1), 
       b1 = intersectCylinder(rayO, rayD, tnear1, tfar1, innerRadius, z0, z1);
  
#if 1
  if (b0 && (!b1))
    rc<SHADING>(dst, nPhi, nNodes, data, grad, bvh, disp, invdet, tf, trans, rayO, rayD, stepsize, tnear0, tfar0);
  else if (b0 && b1) {
    rc<SHADING>(dst, nPhi, nNodes, data, grad, bvh, disp, invdet, tf, trans, rayO, rayD, stepsize, tnear0, tnear1);
    rc<SHADING>(dst, nPhi, nNodes, data, grad, bvh, disp, invdet, tf, trans, rayO, rayD, stepsize, tfar1, tfar0);
  }
#else
  if (b0) {
    rc<SHADING>(dst, nPhi, nNodes, data, bvh, disp, tf, trans, rayO, rayD, stepsize, tnear0, tfar0);
  }
#endif
}

#if WITH_CUDA
__global__ static void test_point_locator_kernel(
    int *output, 
    float x, float y,
    QuadNodeD *bvh,
    float *invdet)
{
  float3 lambda;
  *output = QuadNodeD_locatePoint(bvh, x, y, lambda, invdet);
}
#endif

__device__ __host__ inline bool setup_ray(
    int *viewport, 
    float *invmvp, 
    uint x, uint y,
    float3 &rayO, float3 &rayD)
{
  float coord[4], obj0[4], obj1[4]; 
  coord[0] = (x-viewport[0])*2.f / viewport[2] - 1.f; 
  coord[1] = (y-viewport[1])*2.f / viewport[3] - 1.f; 
  coord[2] = -1.0; 
  coord[3] = 1.0;

  mulmatvec(invmvp, coord, obj0); 
  coord[2] = 1.0; 
  mulmatvec(invmvp, coord, obj1); 
  if (obj0[3] == 0.f || obj1[3] == 0.f) return false; 

  for (int i=0; i<3; i++)
      obj0[i] /= obj0[3], obj1[i] /= obj1[3]; 

  rayO = make_float3(obj0[0], obj0[1], obj0[2]);
  rayD = normalize(make_float3(obj1[0]-obj0[0], obj1[1]-obj0[1], obj1[2]-obj0[2]));
  return true;
}

#if WITH_CUDA
template <int SHADING>
__global__ static void raycasting_kernel(
        unsigned char *output_rgba8,
        int *viewport, 
        float *invmvp,
        int nPhi, 
        int nNodes, 
        float *data, 
        float *grad,
        QuadNodeD *bvh,
        float *disp,
        float *invdet,
        float *tf,
        float2 trans, 
        float stepsize)
{
  uint x = blockIdx.x*blockDim.x + threadIdx.x;
  uint y = blockIdx.y*blockDim.y + threadIdx.y;

  if (x >= viewport[2] || y>= viewport[3]) return;
 
  float3 rayO, rayD;
  setup_ray(viewport, invmvp, x, y, rayO, rayD);

  float4 dst = make_float4(0.f); 
  raycasting<SHADING>(dst, nPhi, nNodes, data, grad, bvh, disp, invdet, tf, trans, rayO, rayD, stepsize);

  output_rgba8[(y*viewport[2]+x)*4+0] = dst.x * 255;
  output_rgba8[(y*viewport[2]+x)*4+1] = dst.y * 255;
  output_rgba8[(y*viewport[2]+x)*4+2] = dst.z * 255;
  output_rgba8[(y*viewport[2]+x)*4+3] = dst.w * 255;

  // if (y == 300) printf("%f\n", dst.w);

#if 0
  // GL_ONE_MINUS_DST_ALPHA, GL_ONE
  float w0 = 1-output[(y*viewport[2]+x)*4+3]; //, w1 = 1; make the compiler happy :)
  output[(y*viewport[2]+x)*4+0] += w0* dst.x;
  output[(y*viewport[2]+x)*4+1] += w0* dst.y;
  output[(y*viewport[2]+x)*4+2] += w0* dst.z;
  output[(y*viewport[2]+x)*4+3] += w0* dst.w;
#endif
}
#endif

template <int SHADING>
static void raycasting_cpu(
        unsigned char *output_rgba8,
        int *viewport, 
        float *invmvp,
        int nPhi, 
        int nNodes, 
        float *data, 
        float *grad,
        QuadNodeD *bvh,
        float *disp,
        float *invdet,
        float *tf,
        float2 trans, 
        float stepsize)
{
  fprintf(stderr, "[volren] CPU rendering, width=%d, height=%d\n", viewport[2], viewport[3]);
#pragma omp parallel for collapse(2)
  for (uint x = 0; x < viewport[2]; x ++) {
    for (uint y = 0; y < viewport[3]; y ++) {
      float3 rayO, rayD;
      setup_ray(viewport, invmvp, x, y, rayO, rayD);

      float4 dst = make_float4(0.f); 
      raycasting<SHADING>(dst, nPhi, nNodes, data, grad, bvh, disp, invdet, tf, trans, rayO, rayD, stepsize);

      output_rgba8[(y*viewport[2]+x)*4+0] = clamp(dst.x, 0.f, 1.f) * 255;
      output_rgba8[(y*viewport[2]+x)*4+1] = clamp(dst.y, 0.f, 1.f) * 255;
      output_rgba8[(y*viewport[2]+x)*4+2] = clamp(dst.z, 0.f, 1.f) * 255;
      output_rgba8[(y*viewport[2]+x)*4+3] = clamp(dst.w, 0.f, 1.f) * 255;
    }
  }
}

/////////////////////////////
extern "C" {

void rc_test_point_locator(ctx_rc *ctx, float x, float y)
{
#if WITH_CUDA
  test_point_locator_kernel<<<1, 1>>>(
      (int*)ctx->d_output_rgba8, 
      x, y, ctx->d_bvh, ctx->d_invdet);

  int nid;
  cudaMemcpy(&nid, ctx->d_output_rgba8, sizeof(int), cudaMemcpyDeviceToHost);
  fprintf(stderr, "[rc_test_point_locator] x={%f, %f}, nid=%d\n", x, y, nid);
#endif
}

void rc_render(ctx_rc *ctx)
{
#if WITH_CUDA
  const dim3 blockSize(16, 16); 
  const dim3 gridSize = dim3(iDivUp(ctx->viewport[2], blockSize.x), iDivUp(ctx->viewport[3], blockSize.y));

  cudaMemcpy(ctx->d_viewport, ctx->viewport, sizeof(int)*4, cudaMemcpyHostToDevice);
  // cudaMemcpyToSymbol(c_viewport, ctx->viewport, sizeof(int)*4);
  cudaMemcpy(ctx->d_invmvp, ctx->invmvp, sizeof(float)*16, cudaMemcpyHostToDevice);
  // cudaMemcpyToSymbol(c_invmvp, ctx->invmvp, sizeof(float)*16);
 
  raycasting_kernel<0><<<gridSize, blockSize>>>(
          ctx->d_output_rgba8,
          ctx->d_viewport, 
          ctx->d_invmvp,
          ctx->nPhi, 
          ctx->nNodes,
          ctx->d_data, 
          ctx->d_grad,
          ctx->d_bvh,
          ctx->d_disp,
          ctx->d_invdet,
          ctx->d_tf,
          make_float2(ctx->trans[0], ctx->trans[1]), 
          ctx->stepsize);

  cudaDeviceSynchronize();
  checkLastCudaError("[rc_render]");
#else
  rc_render_cpu(ctx);
#endif
}

void rc_render_cpu(ctx_rc *ctx)
{
  raycasting_cpu<0>(
          (unsigned char*)ctx->h_output,
          ctx->viewport, 
          ctx->invmvp,
          ctx->nPhi, 
          ctx->nNodes,
          ctx->h_data, 
          ctx->h_grad,
          ctx->h_bvh,
          ctx->h_disp,
          ctx->h_invdet,
          ctx->h_tf,
          make_float2(ctx->trans[0], ctx->trans[1]), 
          ctx->stepsize);
}

void rc_bind_bvh(ctx_rc *ctx, int nQuadNodes, QuadNodeD *bvh)
{
  ctx->h_bvh = bvh;
#if WITH_CUDA
  if (ctx->d_bvh != NULL)
    cudaFree(ctx->d_bvh);

  cudaMalloc((void**)&ctx->d_bvh, sizeof(QuadNodeD)*nQuadNodes);
  cudaMemcpy(ctx->d_bvh, bvh, sizeof(QuadNodeD)*nQuadNodes, cudaMemcpyHostToDevice);
#endif
}

void rc_set_default_tf(ctx_rc *ctx)
{
  float r[3] = {0.7215686274509804, 0.1803921568627451, 0.1803921568627451}, 
        b[3] = {0.2, 0.4, 0.8};

  float *tf = ctx->h_tf;
  for (int i=0; i<size_tf; i++) {
    float x = (float)i / (size_tf-1);
    tf[i*4] = x*b[0] + (1-x)*r[0];
    tf[i*4+1] = x*b[1] + (1-x)*r[1];
    tf[i*4+2] = x*b[2] + (1-x)*r[2];
    tf[i*4+3] = fmin(0.999f, (x-0.5)*(x-0.5)*40);
  }
  rc_set_tf(ctx, tf);
}

void rc_set_tf(ctx_rc *ctx, float *tf) 
{
  if (tf != ctx->h_tf)
    memcpy(ctx->h_tf, tf, sizeof(float)*size_tf*4);

#if WITH_CUDA
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();

  cudaMemcpy(ctx->d_tf, tf, sizeof(float)*size_tf*4, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( ctx->d_tfArray, 0, 0, tf, size_tf*4, cudaMemcpyHostToDevice ); 

  texTransferFunc.normalized = true; 
  texTransferFunc.filterMode = cudaFilterModeLinear; 
  texTransferFunc.addressMode[0] = cudaAddressModeClamp; 
  cudaBindTextureToArray(texTransferFunc, ctx->d_tfArray, channelDesc); 

  checkLastCudaError("[rc_set_tf]");
#endif
}

void rc_bind_disp(ctx_rc *ctx, int nNodes, float *disp)
{
  ctx->h_disp = disp;
#if WITH_CUDA
  if (ctx->d_disp == NULL)
    cudaMalloc((void**)&ctx->d_disp, sizeof(float)*nNodes*2);
  cudaMemcpy(ctx->d_disp, disp, sizeof(float)*nNodes*2, cudaMemcpyHostToDevice);
#endif
}

void rc_bind_invdet(ctx_rc *ctx, int nTriangles, float *invdet)
{
  ctx->h_invdet = invdet;
#if WITH_CUDA
  if (ctx->d_invdet == NULL)
    cudaMalloc((void**)&ctx->d_invdet, sizeof(float)*nTriangles);
  cudaMemcpy(ctx->d_invdet, invdet, sizeof(float)*nTriangles, cudaMemcpyHostToDevice);
#endif
}

void rc_bind_data(ctx_rc *ctx, int nNodes, int nTriangles, int nPhi, float *data, float *grad)
{
  ctx->h_data = data;
  ctx->h_grad = grad;
  ctx->nNodes = nNodes;
  ctx->nPhi = nPhi;
#if WITH_CUDA
  if (ctx->d_data == NULL)
    cudaMalloc((void**)&ctx->d_data, sizeof(float)*nNodes*nPhi);
  cudaMemcpy(ctx->d_data, data, sizeof(float)*nNodes*nPhi, cudaMemcpyHostToDevice);
  if (ctx->d_grad == NULL)
    cudaMalloc((void**)&ctx->d_grad, sizeof(float)*nNodes*nTriangles*2);
  cudaMemcpy(ctx->d_grad, grad, sizeof(float)*nNodes*nTriangles*2, cudaMemcpyHostToDevice);
#endif
}

void rc_create_ctx(ctx_rc **ctx)
{
  *ctx = (ctx_rc*)malloc(sizeof(ctx_rc));
  memset(*ctx, 0, sizeof(ctx_rc));

  (*ctx)->h_tf = (float*)malloc(sizeof(float)*size_tf*4);

  const size_t max_npx = 4096*4096;
  (*ctx)->h_output = malloc(4*max_npx);

#if WITH_CUDA
  cudaSetDevice(0);
  cudaMalloc((void**)&((*ctx)->d_output_rgba8), 4*max_npx); 
  
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
  cudaMallocArray( &(*ctx)->d_tfArray, &channelDesc, size_tf*4, 1 ); 

  cudaMalloc((void**)&((*ctx)->d_tf), sizeof(float)*size_tf*4);
  cudaMalloc((void**)&((*ctx)->d_viewport), sizeof(int)*4);
  cudaMalloc((void**)&((*ctx)->d_invmvp), sizeof(float)*16);
  
  checkLastCudaError("[rc_init]");
#endif
}

void rc_destroy_ctx(ctx_rc **ctx)
{
  // TODO: free any resources
#if WITH_CUDA
  cudaFree((*ctx)->d_output_rgba8);
  cudaFree((*ctx)->d_disp);
  cudaFree((*ctx)->d_invdet);
#endif
  free((*ctx)->h_output);
  free((*ctx)->h_tf);
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
  memset(ctx->h_output, 0, ctx->viewport[2]*ctx->viewport[3]*4);
#if WITH_CUDA
  cudaMemset(ctx->d_output_rgba8, 0, 4*ctx->viewport[2]*ctx->viewport[3]);
#endif
}

void rc_copy_output_to_host(ctx_rc *ctx)
{
#if WITH_CUDA
  cudaMemcpy(ctx->h_output, ctx->d_output_rgba8, 4*ctx->viewport[2]*ctx->viewport[3], cudaMemcpyDeviceToHost);
#if 0
  unsigned char* img = (unsigned char*)ctx->h_output;
  const int j = 0;
  for (int i=0; i<720; i++) {
    fprintf(stderr, "%f\n", img[(j*720+i)*4+3]);
  }
#endif
#endif
}

} // extern "C" 
/////////////////
