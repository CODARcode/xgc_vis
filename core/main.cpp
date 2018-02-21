#include <mpi.h>
#include <adios.h>
#include <adios_read.h>
#include <adios_error.h>
#include <glob.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <boost/program_options.hpp>
#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>
#include <concurrentqueue.h>
#include "core/bp_utils.hpp"
#include "core/xgcBlobExtractor.h"

#ifdef VOLREN
#include "volren/bvh.h"
#include "volren/volren.cuh"
#endif

const std::string pointsName = "points";
const std::string numPointsName = "numPoints";
const std::string cellsName = "cells";
const std::string numCellsName = "numCells";
const std::string cellShape = "triangle"; // FIXME
const std::string dpotName = "dpot";
const std::string psiName = "psi";
const std::string labelsName = "labels";
const std::string centering = "point";

const std::string groupName = "xgc_blobs", meshName = "xgc_mesh2D";
  
typedef websocketpp::server<websocketpp::config::asio> server;
typedef server::message_ptr message_ptr;

// using ConcurrentQueue = moodycamel::ConcurrentQueue;
using json = nlohmann::json;
using websocketpp::lib::placeholders::_1;
using websocketpp::lib::placeholders::_2;
using websocketpp::lib::bind;

server wss;

XGCBlobExtractor *ex = NULL;
std::mutex mutex_ex;

float *framebuf = (float*)malloc(sizeof(float)*4*2048*2048);

#if 0
struct VolrenTask {
  int viewport[4];
  double invmvpd[16];
  float *fb = NULL; // framebuffer
};

struct SendTask {

};
#endif

std::set<size_t> loadSkipList(const std::string& filename)
{
  std::ifstream ifs(filename, std::ifstream::in);
  json j;
  ifs >> j;
  ifs.close();

  std::set<size_t> skip;
  for (json::iterator it = j.begin();  it != j.end();  it ++) {
    size_t t = *it;
    skip.insert(t);
    // fprintf(stderr, "will skip time step %lu\n", t);
  }
  return skip;
}

void onMessage(server* s, websocketpp::connection_hdl hdl, message_ptr msg) {
  std::cout << "onMessage called with hdl: " << hdl.lock().get()
            << " and message: " << msg->get_payload()
            << std::endl;

  // check for a special command to instruct the server to stop listening so
  // it can be cleanly exited.

  json incoming = json::parse(msg->get_payload());
  json outgoing;
  bool binary = false;
  char *buffer = NULL;
  size_t buffer_length;
 
  if (incoming["type"] == "requestMesh") {
    fprintf(stderr, "requesting mesh!\n");
    outgoing["type"] = "mesh";
    outgoing["data"] = ex->jsonfyMesh();
  } else if (incoming["type"] == "requestSingleSliceRawData") {
    binary = true;
    const int nNodes = ex->getNNodes();
    const double *array0 = ex->getData();
    buffer_length = sizeof(int) + nNodes*sizeof(float);
    buffer = (char*)malloc(buffer_length);
    int *type = (int*)buffer;
    float *array = (float*)(buffer + sizeof(int));
    *type = 10;
    for (int i=0; i<nNodes; i++) 
      array[i] = array0[i];
  } else if (incoming["type"] == "requestMultipleSliceRawData") {
    binary = true;
    const int nNodes = ex->getNNodes();
    const int nPhi = ex->getNPhi();
    const double *array0 = ex->getData();
    buffer_length = sizeof(int) + nNodes*nPhi*sizeof(float);
    buffer = (char*)malloc(buffer_length);
    int *type = (int*)buffer;
    float *array = (float*)(buffer + sizeof(int));
    *type = 11;
    for (int i=0; i<nNodes*nPhi; i++) 
      array[i] = array0[i];
  } else if (incoming["type"] == "requestData") {
    // int client_current_time_index = incoming["client_current_time_index"];
    outgoing["type"] = "data";

    mutex_ex.lock();
    const size_t timestep = ex->getTimestep();
    const int nnodes = ex->getNNodes();
    const double *ptr = ex->getData(); 
    std::vector<double> dpot(ptr, ptr + nnodes);
    // std::vector<int> labels = ex->getFlattenedLabels(0); // ex->getLabels(0);
    std::vector<int> labels = ex->getLabels(0);
    mutex_ex.unlock();
   
    outgoing["timestep"] = timestep;
    outgoing["data"] = dpot;
    outgoing["labels"] = labels;
  } else if (incoming["type"] == "requestVolren") {
    binary = true;
    buffer_length = 1024*768*4*sizeof(float);
    buffer = (char*)malloc(buffer_length);
    for (int i=0; i<1024*768; i++) 
      buffer[i*4] = 1;
    // memcpy(buffer, framebuf, buffer_length);
  } else if (incoming["type"] == "stopServer") {
    s->stop_listening();
    return;
  }

  try {
    // s->send(hdl, msg->get_payload(), msg->get_opcode()); // echo
    if (binary) {
      s->send(hdl, buffer, buffer_length, websocketpp::frame::opcode::binary);
      free(buffer);
    }
    else
      s->send(hdl, outgoing.dump(), websocketpp::frame::opcode::text);
  } catch (const websocketpp::lib::error_code& e) {
    std::cout << "Sending failed failed because: " << e
              << "(" << e.message() << ")" << std::endl;
  }
}

void startWebsocketServer(int port)
{
  try {
    // Set logging settings
    wss.set_access_channels(websocketpp::log::alevel::all);
    wss.clear_access_channels(websocketpp::log::alevel::frame_payload);

    // Initialize Asio
    wss.init_asio();

    // Register our message handler
    wss.set_message_handler(bind(&onMessage, &wss, ::_1, ::_2));

    // Listen on port 9002
    wss.listen(port);

    // Start the server accept loop
    wss.start_accept();
#if 0 
    char hostname[1024];
    gethostname(hostname, 1024);
    fprintf(stderr, "=================WEBSOCKET================\n");
    // fprintf(stderr, "hostname=%s\n", hostname);
    fprintf(stderr, "In order to view analysis results, you need to follow two steps:\n");
    fprintf(stderr, " 1. Create a SSH tunnel on your machine:\n\n");
    fprintf(stderr, "  $ ssh -L %d:%s:%d titan-internal.ccs.ornl.gov -N\n\n", port, hostname, port);
    fprintf(stderr, " On PASSWORD, you need to enter your PIN and 6-digits numbers on your token.\n\n");
    fprintf(stderr, " 2. Open your web browser (we recommend Google Chrome or Safari), and open the following webpage:\n\n   http://www.mcs.anl.gov/~hguo/xgc\n\n");
    fprintf(stderr, "If you have any technical difficulties, please contact Hanqi Guo (hguo@anl.gov) directly.\n");
    fprintf(stderr, "==========================================\n");
#endif
    // Start the ASIO io_service run loop
    wss.run();
  } catch (websocketpp::exception const & e) {
    std::cerr << e.what() << std::endl;
    exit(EXIT_FAILURE);
  } catch (...) {
    std::cerr << "other exception" << std::endl;
  }
}

void startVolren(int nNodes, int nTriangles, double *coords, int *conn)
{
#ifdef VOLREN
  std::vector<QuadNodeD> bvh = buildBVHGPU(nNodes, nTriangles, coords, conn);

  ctx_rc *rc;
  rc_create_ctx(&rc);
  rc_set_stepsize(rc, 0.5);
  rc_set_viewport(rc, 0, 0, 1024, 768);

  rc_clear_output(rc);
  // rc_set_invmvpd(rc, imvmvpd);
  fprintf(stderr, "rendering...\n");
  rc_render(rc);
  for (int i=0; i<1024*768; i++) 
    framebuf[i*4] = 1;
  rc_dump_output(rc, framebuf);
  fprintf(stderr, "done.\n");
#endif
}

void writeUnstructredMeshDataFile(int timestep, MPI_Comm comm, int64_t groupHandle, const std::string& fileName, const std::string& writeMethod, const std::string& writeMethodParams,
    int nNodes, int nTriangles, double *coords, int *conn_, double *dpot, double *psi, int *labels)
{
  int64_t fileHandle = -1;
  if (writeMethod == "POSIX" || writeMethod == "MPI")
    adios_open(&fileHandle, groupName.c_str(), fileName.c_str(), "w", comm);
  else 
    adios_open(&fileHandle, groupName.c_str(), fileName.c_str(), (timestep == 1 ? "w" : "a"), comm);

  // fprintf(stderr, "groupHandle=%lld, fileHandle=%lld\n", groupHandle, fileHandle);

  adios_define_mesh_unstructured(
      (char*)pointsName.c_str(), 
      (char*)cellsName.c_str(), 
      (char*)numCellsName.c_str(), 
      (char*)cellShape.c_str(), 
      (char*)numPointsName.c_str(), 
      (char*)"2",
      groupHandle, 
      meshName.c_str());
 
  // if (timestep == 1) {
  if (1) {
    // points
    adios_define_var(groupHandle, numPointsName.c_str(), "", adios_integer, 0, 0, 0);
    adios_write(fileHandle, numPointsName.c_str(), &nNodes);

    int64_t ptId = adios_define_var(
        groupHandle, 
        pointsName.c_str(), 
        "", 
        adios_double, 
        std::string(numPointsName + ",2").c_str(), 
        std::string(numPointsName + ",2").c_str(), 
        "0,0");
    adios_write_byid(fileHandle, ptId, &coords[0]);

    // cells
    adios_define_var(groupHandle, numCellsName.c_str(), "", adios_integer, 0, 0, 0);
    adios_write(fileHandle, numCellsName.c_str(), &nTriangles);

    const int ptsInCell = 3;
    std::string cellDim = std::to_string(nTriangles) + "," + std::to_string(ptsInCell);
    int64_t cellId = adios_define_var(
        groupHandle, 
        cellsName.c_str(), 
        "", 
        adios_integer, 
        cellDim.c_str(),
        cellDim.c_str(), 
        "0,0");
    adios_write_byid(fileHandle, cellId, &conn_[0]);
    
    // psi
    if (psi != NULL) {
      adios_define_var_mesh(groupHandle, psiName.c_str(), meshName.c_str());
      adios_define_var_centering(groupHandle, psiName.c_str(), centering.c_str());
      int64_t psiId = adios_define_var(
          groupHandle, 
          psiName.c_str(), 
          "",
          adios_double, 
          std::to_string(nNodes).c_str(), 
          std::to_string(nNodes).c_str(),
          "0");
      adios_write_byid(fileHandle, psiId, psi);
    }
  }

  // dpot
  adios_define_var_mesh(groupHandle, dpotName.c_str(), meshName.c_str());
  adios_define_var_centering(groupHandle, dpotName.c_str(), centering.c_str());
  int64_t dpotId = adios_define_var(
      groupHandle, 
      dpotName.c_str(), 
      "",
      adios_double, 
      std::to_string(nNodes).c_str(), 
      std::to_string(nNodes).c_str(),
      "0");
  adios_write_byid(fileHandle, dpotId, dpot);
  
  // labels
  if (labels != NULL) {
    double *d_labels = (double*)malloc(sizeof(double)*nNodes);
    for (int i=0; i<nNodes; i++) 
      d_labels[i] = static_cast<double>(labels[i]);

    adios_define_var_mesh(groupHandle, labelsName.c_str(), meshName.c_str());
    adios_define_var_centering(groupHandle, labelsName.c_str(), centering.c_str());
    int64_t labelsId = adios_define_var(
        groupHandle, 
        labelsName.c_str(), 
        "",
        adios_double, // adios_integer, 
        std::to_string(nNodes).c_str(), 
        std::to_string(nNodes).c_str(),
        "0");
    adios_write_byid(fileHandle, labelsId, d_labels);
    free(d_labels);
  }

  adios_close(fileHandle);
  // adios_finalize(0);
}


int main(int argc, char **argv)
{
  std::thread *ws_thread = NULL;
  std::thread *volren_thread = NULL;

  MPI_Init(&argc, &argv);
  
  using namespace boost::program_options;
  options_description opts(argv[0]);
  opts.add_options()
    ("mesh,m", value<std::string>(), "mesh_file")
    ("input,i", value<std::string>(), "input_file")
    ("pattern", value<std::string>(), "input file pattern, e.g. 'xgc.3d.*.bp'.  only works for BP.")
    ("output,o", value<std::string>()->default_value(""), "output_file")
    ("output_prefix", value<std::string>()->default_value(""), "output file prefix (e.g. 'features' will generate features.%05d.bp)")
    ("output_branches,o", value<std::string>()->default_value(""), "output_branches")
    ("output_prefix_branches,o", value<std::string>()->default_value(""), "output_prefix_branches")
    ("read_method,r", value<std::string>()->default_value("BP"), "read_method (BP|DATASPACES|DIMES|FLEXPATH)")
    ("write_method,w", value<std::string>()->default_value("BIN"), "write_method (POSIX|MPI|DIMES|BIN), bin for raw binary")
    ("write_method_params", value<std::string>()->default_value(""), "write_method_params")
    ("skip", value<std::string>(), "skip timesteps that are specified in a json file")
    ("server,s", "enable websocket server")
    ("volren,v", "enable volren (-s required)")
    ("port,p", value<int>()->default_value(9002), "websocket server port")
    ("help,h", "display this information");
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opts)
      .style(command_line_style::unix_style | command_line_style::allow_long_disguise)
      .run(), vm);
  notify(vm);

  if (!vm.count("mesh") || vm.count("h") || !(vm.count("input") || vm.count("pattern"))) {
    std::cout << opts << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  std::vector<std::string> input_filename_list;
  if (vm.count("pattern")) {
    const std::string pattern = vm["pattern"].as<std::string>();
    glob_t results; 
    glob(pattern.c_str(), 0, NULL, &results); 
    for (int i=0; i<results.gl_pathc; i++)
      input_filename_list.push_back(results.gl_pathv[i]); 
    globfree(&results);

    if (input_filename_list.empty()) {
      fprintf(stderr, "FATAL: cannot find any data file matching '%s'\n", pattern.c_str());
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  std::set<size_t> skipped_timesteps;
  if (vm.count("skip")) {
    const std::string filename = vm["skip"].as<std::string>();
    skipped_timesteps = loadSkipList(filename);
  }

  const std::string filename_mesh = vm["mesh"].as<std::string>();
  const std::string output_prefix = vm["output_prefix"].as<std::string>();
  const std::string output_prefix_branches = vm["output_prefix_branches"].as<std::string>();
  const std::string read_method_str = vm["read_method"].as<std::string>();
  const std::string write_method_str = vm["write_method"].as<std::string>();
  const std::string write_method_params_str = vm["write_method_params"].as<std::string>();
  std::string filename_output = vm["output"].as<std::string>();
  std::string filename_output_branches = vm["output_branches"].as<std::string>();
  std::string filename_input;
  bool single_input = false;
  if (vm.count("input")) {
    filename_input = vm["input"].as<std::string>();
    single_input = true;
  }

  bool write_binary = (vm["write_method"].as<std::string>() == "BIN");
  bool volren = vm.count("volren") && vm.count("server");

  fprintf(stderr, "==========================================\n");
  fprintf(stderr, "filename_mesh=%s\n", filename_mesh.c_str());
  if (single_input) {
    fprintf(stderr, "filename_input=%s\n", filename_input.c_str());
  } else {
    fprintf(stderr, "filename_input_pattern=%s\n", vm["pattern"].as<std::string>().c_str());
    for (int i=0; i<input_filename_list.size(); i++) 
      fprintf(stderr, "filename_input_%d=%s\n", i, input_filename_list[i].c_str());
  }
  if (filename_output.size() > 0) 
    fprintf(stderr, "filename_output=%s\n", filename_output.c_str());
  else if (output_prefix.size() > 0)
    fprintf(stderr, "output_prefix=%s\n", output_prefix.c_str());
  else
    fprintf(stderr, "filename_output=(null)\n");
  fprintf(stderr, "read_method=%s\n", read_method_str.c_str());
  fprintf(stderr, "write_method=%s\n", write_method_str.c_str());
  fprintf(stderr, "write_method_params=%s\n", write_method_params_str.c_str());
  fprintf(stderr, "==========================================\n");

  ADIOS_READ_METHOD read_method;
  if (read_method_str == "BP") read_method = ADIOS_READ_METHOD_BP;
  else if (read_method_str == "DATASPACES") read_method = ADIOS_READ_METHOD_DATASPACES;
  else if (read_method_str == "DIMES") read_method = ADIOS_READ_METHOD_DIMES;
  else if (read_method_str == "FLEXPATH") read_method = ADIOS_READ_METHOD_FLEXPATH;
  else {
    fprintf(stderr, "unsupported read method: %s\n", read_method_str.c_str());
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  /// read mesh
  adios_read_init_method(read_method, MPI_COMM_WORLD, "");
  ADIOS_FILE *meshFP = NULL;
  while (1) {
    meshFP = adios_read_open_file(filename_mesh.c_str(), ADIOS_READ_METHOD_BP, MPI_COMM_WORLD); // always use ADIOS_READ_METHOD_BP for mesh
    if (meshFP == NULL) {
      fprintf(stderr, "failed to open mesh: %s, will retry in 5 seconds.\n", filename_mesh.c_str()); 
      sleep(5);
    } else break;
  }
  adios_read_bp_reset_dimension_order(meshFP, 0);
  
  int nNodes, nTriangles;

  fprintf(stderr, "reading mesh...\n");
  double *coords;
  int *conn, *nextNode;
  double *psi;
  readTriangularMesh(meshFP, nNodes, nTriangles, &coords, &conn, &nextNode);
  readScalars<double>(meshFP, "psi", &psi);
  // adios_close(*meshFP);

  ex = new XGCBlobExtractor(nNodes, nTriangles, coords, conn);
  
  // starting server
  if (vm.count("server")) {
    ws_thread = new std::thread(startWebsocketServer, vm["port"].as<int>());
  }

  // starting volren
  if (volren) {
    volren_thread = new std::thread(startVolren, nNodes, nTriangles, coords, conn);
  }

  // read data
  fprintf(stderr, "opening data stream...\n");
  ADIOS_FILE *varFP;
  if (single_input > 0) {
    if (read_method == ADIOS_READ_METHOD_BP)
      varFP = adios_read_open_file(filename_input.c_str(), read_method, MPI_COMM_WORLD);
    else 
      varFP = adios_read_open (filename_input.c_str(), read_method, MPI_COMM_WORLD, ADIOS_LOCKMODE_ALL, -1.0);
    fprintf(stderr, "varFP=%p, errno=%d\n", varFP, err_end_of_stream);
  }
  adios_read_bp_reset_dimension_order(varFP, 0);
    
  double *dpot = NULL;
  size_t current_time_index = 0; // only for multiple inputs.

  // output
  int64_t groupHandle = -1;
  if (!write_binary) {
    adios_init_noxml(MPI_COMM_WORLD);
    adios_declare_group(&groupHandle, groupName.c_str(), "", adios_stat_default);
    adios_select_method(groupHandle, write_method_str.c_str(), "", "");
    adios_define_schema_version(groupHandle, (char*)"1.1");
    adios_define_mesh_timevarying("no", groupHandle, meshName.c_str());
    adios_delete_vardefs(groupHandle);
  }

  while (1) {
    if (single_input) {
      if (adios_errno == err_end_of_stream) break;
      current_time_index ++;
      fprintf(stderr, "reading data, time_index=%lu\n", current_time_index);
    } else { // multiple input;
      if (current_time_index >= input_filename_list.size()) {
        fprintf(stderr, "all done.\n");
        break;
      }
      const std::string &current_input_filename = input_filename_list[ current_time_index ++ ];
      fprintf(stderr, "reading data from %s\n", current_input_filename.c_str());
      varFP = adios_read_open_file(current_input_filename.c_str(), read_method, MPI_COMM_WORLD);
    }

    int nPhi = 1;
    readValueInt(varFP, "nphi", &nPhi);
    // readScalars<double>(varFP, "dpot", &dpot);

    ADIOS_VARINFO *avi = adios_inq_var(varFP, "dpot");
    assert(avi != NULL);

    adios_inq_var_stat(varFP, avi, 0, 0);
    adios_inq_var_blockinfo(varFP, avi);
    adios_inq_var_meshinfo(varFP, avi);

    // uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {static_cast<uint64_t>(nNodes), static_cast<uint64_t>(nPhi), 1, 1};
    uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {static_cast<uint64_t>(nPhi), static_cast<uint64_t>(nNodes), 1, 1};
    ADIOS_SELECTION *sel = adios_selection_boundingbox(avi->ndim, st, sz);
    // fprintf(stderr, "%d, {%d, %d, %d, %d}\n", avi->ndim, avi->dims[0], avi->dims[1], avi->dims[2], avi->dims[3]);

    assert(sel->type == ADIOS_SELECTION_BOUNDINGBOX);

    if (dpot == NULL) 
      dpot = (double*)malloc(sizeof(double)*nPhi*nNodes);

    adios_schedule_read_byid(varFP, sel, avi->varid, 0, 1, dpot);
    adios_perform_reads(varFP, 1);
    adios_selection_delete(sel);

    if (skipped_timesteps.find(current_time_index) != skipped_timesteps.end()) {
      fprintf(stderr, "skipping timestep %lu.\n", current_time_index);
      adios_advance_step(varFP, 0, 1.0);
      continue;
    }

    fprintf(stderr, "starting analysis..\n");

    mutex_ex.lock();
    ex->setData(current_time_index, nPhi, dpot);
    // ex->setPersistenceThreshold(persistence_threshold);
    // ex->buildContourTree3D();
    // std::map<ctBranch*, size_t> branchSet = ex->buildContourTree2D(0);
    ex->buildContourTree2DAll();
    mutex_ex.unlock();

    // write labels
    int *labels = ex->getLabels(0).data();
    if (output_prefix.length() > 0) {
      std::stringstream ss_filename;
      ss_filename << output_prefix << "." << std::setfill('0') << std::setw(5) << current_time_index 
        << (write_binary ? ".bin" : ".bp");
      filename_output = ss_filename.str();
    }
      
    if (filename_output.length() > 0) {
      fprintf(stderr, "writing results for timestep %zu to %s\n", 
          current_time_index, filename_output.c_str());
      if (write_binary)
        ex->dumpLabels(filename_output);
      else 
        writeUnstructredMeshDataFile(current_time_index, MPI_COMM_WORLD, groupHandle, filename_output, write_method_str, write_method_params_str,
            nNodes, nTriangles, coords, conn, dpot, psi, labels);
    }
   
    // write branches
    if (output_prefix_branches.length() > 0) {
      std::stringstream ss_filename;
      ss_filename << output_prefix_branches << "." << std::setfill('0') << std::setw(5) << current_time_index << ".json";
      filename_output_branches = ss_filename.str();
    }

    if (filename_output_branches.length() > 0) {
      fprintf(stderr, "writing branch decompositions for timestep %zu to %s\n", 
          current_time_index, filename_output_branches.c_str()); 
      // ex->dumpBranches(filename_output_branches, branchSet);
    }
    
    fprintf(stderr, "done.\n");

    // sleep(6);
    if (single_input) {
      if (read_method == ADIOS_READ_METHOD_BP) break;
      else adios_advance_step(varFP, 0, 1.0);
    }
  }

#if 0
  if (filename_output.length() > 0) {
    adios_close(fileHandleOutput);
    adios_finalize(0);
  }
#endif
  
  if (ws_thread) {
    fprintf(stderr, "shutting down wss...\n");
    ws_thread->join();
  }

  if (volren_thread) {
    fprintf(stderr, "shutting down volren...\n");
    volren_thread->join();
  }
   
  // adios_close(*varFP);

  if (dpot != NULL) free(dpot);
  delete ex;
  ex = NULL;
  free(coords);
  free(conn);
  free(nextNode);
  adios_finalize(0);

  fprintf(stderr, "exiting...\n");
  MPI_Finalize();
  return 0;
}
