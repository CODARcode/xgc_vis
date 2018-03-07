#ifndef __VOLRENENGINE_H
#define __VOLRENENGINE_H

#include <queue>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <thread>
#include "volren/png_utils.hpp"
#include "volren/bvh.h"
#include "volren/volren.cuh"
#include "volren/volrenEngine.h"
#include "io/xgcMesh.h"
#include "io/xgcData.h"

enum {
  VOLREN_RENDER = 0,
  VOLREN_EXIT = 1
};

enum {
  VOLREN_FORMAT_RGBA8 = 0,
  VOLREN_FORMAT_PNG = 1
};

struct VolrenTask {
  int tag = VOLREN_RENDER;
  int format = VOLREN_FORMAT_PNG;
  int viewport[4];
  double invmvpd[16];
  float *tf = NULL;
  png_mem_buffer png;
  unsigned char *fb = NULL;
  std::condition_variable cond;
  std::mutex mutex;

  bool enable_angle = false;
  float start_angle = 0, end_angle = M_PI;

  float psi_start=0, psi_end=0.2;

  bool enable_shading = true;
  float Ks = 0.2f, Kd = 0.3f, Ka = 0.04f;
  float light_direction[3] = {-1, 0, 0};

  ~VolrenTask();

  static VolrenTask* createVolrenTaskFromString(const std::string& s);
};

struct VolrenEngine {
  VolrenEngine();
  ~VolrenEngine();

  void start(XGCMesh& m, XGCData& d);
  void start_(XGCMesh& m, XGCData& d);
  void stop();
  void enqueueAndWait(VolrenTask *task);
  VolrenTask* enqueueAndWait(const std::string& s);

  bool started() const {return thread != NULL;}
  std::thread *thread = NULL;
  std::queue<VolrenTask*> volrenTaskQueue;
  std::mutex mutex_volrenTaskQueue;
  std::condition_variable cond_volrenTaskQueue;
};

#endif
