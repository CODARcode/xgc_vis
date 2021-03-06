#include "volren/volrenEngine.h"
#include "volren/kdbvh.h"
#include <json.hpp>
#include <execinfo.h>

using json = nlohmann::json;

VolrenTask::~VolrenTask()
{
  if (tf != NULL) free(tf);
  if (fb != NULL) free(fb);
  if (png.buffer != NULL) free(png.buffer);
}

VolrenTask* VolrenEngine::createTaskFromString(const std::string& query)
{
  const int defaulat_viewport[4] = {0, 0, 720, 382};
  // const int defaulat_viewport[4] = {0, 0, 1440, 764};
  const double defaulat_invmvpd[16] = {1.1604, 0.0367814, -0.496243, 0, -0.157143, 0.563898, -0.325661, 0, -9.79775, -16.6755, -24.1467, -4.95, 9.20395, 15.6649, 22.6832, 5.05};
    
  // create task
  VolrenTask *task = new VolrenTask;
  task->str = query;
 
  // parse query
  try {
    json j = json::parse(query);
    fprintf(stderr, "[volren,rank=%d] json parameter: %s\n", rank, j.dump().c_str());

    if (!j["tag"].is_null()) {
      if (j["tag"] == "exit") task->tag = VOLREN_EXIT;
    }

    if (j["width"].is_null() || j["height"].is_null())
      memcpy(task->viewport, defaulat_viewport, sizeof(int)*4);
    else {
      task->viewport[0] = 0;
      task->viewport[1] = 0;
      task->viewport[2] = j["width"];
      task->viewport[3] = j["height"];
    }

    if (j["matrix"].is_null())
      memcpy(task->invmvpd, defaulat_invmvpd, sizeof(double)*16);
    else {
      for (int i=0; i<16; i++) 
        task->invmvpd[i] = j["matrix"][i];
    }

    if (!j["tf"].is_null()) {
      json jtf = j["tf"];
      task->tf = (float*)malloc(sizeof(float)*1024);
      for (int i=0; i<1024; i++) 
        task->tf[i] = jtf.at(i).get<float>() / 255.f;
    }

    if (!j["enableAngle"].is_null()) task->enable_angle = j["enableAngle"].get<bool>();
    if (!j["startAngle"].is_null()) task->start_angle = j["startAngle"].get<float>();
    if (!j["endAngle"].is_null()) task->end_angle = j["endAngle"].get<float>();

    if (!j["enableShading"].is_null()) task->enable_shading = j["enableShading"].get<bool>();
    if (!j["Ka"].is_null()) task->Ka = j["Ka"].get<float>();
    if (!j["Kd"].is_null()) task->Kd = j["Kd"].get<float>();
    if (!j["Ks"].is_null()) task->Ks = j["Ks"].get<float>();
    if (!j["lightingDirectionX"].is_null()) task->light_direction[0] = j["lightingDirectionX"].get<float>();
    if (!j["lightingDirectionY"].is_null()) task->light_direction[1] = j["lightingDirectionY"].get<float>();
    if (!j["lightingDirectionZ"].is_null()) task->light_direction[2] = j["lightingDirectionZ"].get<float>();

    if (j["psiStart"].is_null()) task->psi_start = 0.f; 
    else task->psi_start = j["psiStart"].get<float>();
    
    if (j["psiEnd"].is_null()) task->psi_end = 0.f; 
    else task->psi_end = j["psiEnd"].get<float>();

  } catch (...) {
    fprintf(stderr, "[volren,rank=%d] json parse failed, using defaulat parameters.\n", rank);
    memcpy(task->viewport, defaulat_viewport, sizeof(int)*4);
    memcpy(task->invmvpd, defaulat_invmvpd, sizeof(double)*16);
  }

  return task;
}

VolrenEngine::VolrenEngine()
{
}

VolrenEngine::~VolrenEngine()
{
  stop();
}


static void
print_trace (void)
{
  void *array[10];
  size_t size;
  char **strings;
  size_t i;

  size = backtrace (array, 10);
  strings = backtrace_symbols (array, size);

  printf ("Obtained %zd stack frames.\n", size);

  for (i = 0; i < size; i++)
    printf ("%s\n", strings[i]);

  free (strings);
}

void VolrenEngine::stop()
{
  // fprintf(stderr, "[volren,rank=%d] stopping volren engine!\n", rank);
  // print_trace();

  if (started()) {
    if (rank == 0) {
      fprintf(stderr, "[volren,rank=%d] submitting a exit job!\n", rank);
      json j;
      j["tag"] = "exit";
      VolrenTask *t = enqueueAndWait(j.dump());
      delete t;
    }
    fprintf(stderr, "[volren,rank=%d] joining volren thread\n", rank);
    thread->join();
    delete thread;
    thread = NULL;
  }
}

void VolrenEngine::start(MPI_Comm comm, XGCMesh& m, XGCData& d)
{
  thread = new std::thread(&VolrenEngine::start_, this, comm, std::ref(m), std::ref(d));
}

void VolrenEngine::start_(MPI_Comm comm, XGCMesh& m, XGCData& d)
{
  MPI_Comm_size(comm, &np);
  MPI_Comm_rank(comm, &rank);

  fprintf(stderr, "[volren,rank=%d] building BVH...\n", rank);
  std::vector<BVHNodeD> bvh = buildBVHGPU(m.nNodes, m.nTriangles, m.coords, m.conn);
  // std::vector<BVHNodeD> bvh = buildKDBVHGPU(m);

  fprintf(stderr, "[volren,rank=%d] initialize volren...\n", rank);
  ctx_rc *rc;
  rc_create_ctx(&rc);
  rc_set_range(rc, -100.f, 100.f); // TODO
  // rc_set_range(rc, -50.f, 50.f); // TODO
  rc_set_default_tf(rc);
  // rc_set_stepsize(rc, 0.001);
  rc_set_stepsize(rc, 0.005);
  rc_bind_bvh(rc, bvh.size(), (BVHNodeD*)bvh.data());
  // rc_test_point_locator(rc, 2.3f, -0.4f);
 
  // mesh 
  rc_bind_neighbors(rc, m.nTriangles, &m.neighbors[0]);
  rc_bind_psi(rc, m.nNodes, &m.psif[0], m.psi_min, m.psi_max);
  rc_bind_disp(rc, m.nNodes, &m.dispf[0]);
  rc_bind_invdet(rc, m.nTriangles, &m.invdetf[0]);

  // dpot
  rc_bind_data(rc, m.nNodes, m.nTriangles, m.nPhi, m.iPhi, &d.dpotf[0], &d.graddpotf[0]);

  while (1) { // volren loop
    VolrenTask *task = NULL;
    int str_len;

    // dequeue task
    if (rank == 0) {
      std::unique_lock<std::mutex> mlock(mutex_volrenTaskQueue);
      while (volrenTaskQueue.empty()) {
        cond_volrenTaskQueue.wait(mlock);
      }
      task = volrenTaskQueue.front();
      volrenTaskQueue.pop();

      if (np > 1) { // distributed rendering
        str_len = task->str.size();
        MPI_Bcast(&str_len, 1, MPI_INT, 0, comm);
        MPI_Bcast((char*)task->str.data(), task->str.size(), MPI_CHAR, 0, comm);
      }
    } else {
      MPI_Bcast(&str_len, 1, MPI_INT, 0, comm);
      std::string str;
      str.resize(str_len);
      MPI_Bcast((char*)str.data(), str_len, MPI_CHAR, 0, comm);
      task = createTaskFromString(str);
    }
    
    if (task->tag == VOLREN_RENDER) 
    {
      std::unique_lock<std::mutex> mlock(task->mutex);
     
      typedef std::chrono::high_resolution_clock clock;

      fprintf(stderr, "[volren,rank=%d] rendering...\n", rank);
      auto t0 = clock::now();
      rc_set_viewport(rc, 0, 0, task->viewport[2], task->viewport[3]);
      rc_set_invmvpd(rc, task->invmvpd);
      // rc_set_psi_range(rc, false, 0, 0.2); // TODO: argument
      rc_set_psi_range(rc, true, task->psi_start, task->psi_end); // TODO: argument
      rc_set_angle_range(rc, task->enable_angle, task->start_angle, task->end_angle);
      rc_set_shading(rc, task->enable_shading, task->Ka, task->Kd, task->Ks, 
          task->light_direction[0], task->light_direction[1], task->light_direction[2]);
      rc_set_slice_highlight_ratio(rc, false, 0.999); // TODO: argument
      if (task->tf) 
        rc_set_tf(rc, task->tf);
      // rc_clear_output(rc);
      rc_render(rc);
      // rc_render_cpu(rc);
      auto t1 = clock::now();
      float tt0 = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
      fprintf(stderr, "[volren,rank=%d] volren time: %f ns\n", rank, tt0);
      
      auto t2 = clock::now();
      rc_copy_output_to_host(rc);
      auto t3 = clock::now();
      float tt2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count();
      fprintf(stderr, "[volren,rank=%d] volren download time: %f ns\n", rank, tt2);
      
      // fprintf(stderr, "[volren] converting to png...\n");
      auto t4 = clock::now();
      if (task->format == VOLREN_FORMAT_PNG) {
#if WITH_PNG
        task->png = save_png(task->viewport[2], task->viewport[3], 8, 
            PNG_COLOR_TYPE_RGBA, (unsigned char*)rc->h_output, 4*task->viewport[2], PNG_TRANSFORM_IDENTITY);
#endif
      } else if (task->format == VOLREN_FORMAT_RGBA8) {
        task->fb = (unsigned char*)malloc(task->viewport[2]*task->viewport[3]*4);
        memcpy(task->fb, rc->h_output, task->viewport[2]*task->viewport[3]*4);
      }
      auto t5 = clock::now();
      float tt4 = std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count();
      task->cond.notify_one();
      
      fprintf(stderr, "[volren,rank=%d] png compression time: %f ns, size=%zu\n", rank, tt4, task->png.size);
    } else if (task->tag == VOLREN_EXIT) {
      rc_destroy_ctx(&rc);
      task->cond.notify_one();
      fprintf(stderr, "[volren,rank=%d] exiting...\n", rank);
      return;
    }
  }

#if 0 // for testing
  const int npx = 1024*768;
  unsigned char *fb = (unsigned char*)malloc(npx*3);
  for (int i=0; i<npx; i++) {
    fb[i*3] = 255; 
    fb[i*3+1] = 0;
    fb[i*3+2] = 0;
  }
  volren_png_buffer = save_png(1024, 768, 8, PNG_COLOR_TYPE_RGB, fb, 3*1024, PNG_TRANSFORM_IDENTITY);
  free(fb);

  // FILE *fp = fopen("test1.png", "wb");
  // fwrite(png.buffer, 1, png.size, fp);
  // free(png.buffer);
  // fclose(fp);
#endif
}

VolrenTask* VolrenEngine::enqueueAndWait(const std::string& s)
{
  VolrenTask *task = NULL; 

  // enqueue
  {
    std::unique_lock<std::mutex> mlock(mutex_volrenTaskQueue);
    task = createTaskFromString(s);
    volrenTaskQueue.push(task);
    mlock.unlock();
    cond_volrenTaskQueue.notify_one();
  }

  // wait
  {
    std::unique_lock<std::mutex> mlock(task->mutex);
    task->cond.wait(mlock);
  }

  return task;
}
