#define __host__
#define __device__

#include <cmath>
#include <cstdlib>
#include <cstring>

typedef unsigned int uint;

struct dim3 {
  uint x, y, z;
  dim3(int xx, int yy) : x(xx), y(yy) {}
};

typedef struct {
  int x, y;
} int2;

typedef struct {
  int x, y, z;
} int3;

typedef struct {
  int x, y, z, w;
} int4;

typedef struct {
  uint x, y;
} uint2;

typedef struct {
  uint x, y, z;
} uint3;

typedef struct {
  uint x, y, z, w;
} uint4;

typedef struct {
  float x, y;
} float2;

typedef struct {
  float x, y, z;
} float3;

typedef struct {
  float x, y, z, w;
} float4;

inline static float2 make_float2(float x, float y) {
  float2 f; 
  f.x = x;
  f.y = y;
  return f;
}

inline static float3 make_float3(float x, float y, float z) {
  float3 f; 
  f.x = x;
  f.y = y;
  f.z = z;
  return f;
}

inline static float4 make_float4(float x, float y, float z, float w) {
  float4 f; 
  f.x = x;
  f.y = y;
  f.z = z;
  f.w = w;
  return f;
}

inline static int2 make_int2(int x, int y) {
  int2 f;
  f.x = x;
  f.y = y;
  return f;
}

inline static int3 make_int3(int x, int y, int z) {
  int3 f;
  f.x = x;
  f.y = y;
  f.z = z;
  return f;
}

inline static int4 make_int4(int x, int y, int z, int w) {
  int4 f;
  f.x = x;
  f.y = y;
  f.z = z;
  f.w = w;
  return f;
}

inline static uint2 make_uint2(uint x, uint y) {
  uint2 f;
  f.x = x;
  f.y = y;
  return f;
}

inline static uint3 make_uint3(uint x, uint y, uint z) {
  uint3 f;
  f.x = x;
  f.y = y;
  f.z = z;
  return f;
}

inline static uint4 make_uint4(uint x, uint y, uint z, uint w) {
  uint4 f;
  f.x = x;
  f.y = y;
  f.z = z;
  f.w = w;
  return f;
}


#include "volren.cu"
