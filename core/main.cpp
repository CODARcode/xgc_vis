#include <mpi.h>
#include <adios.h>
#include <adios_read.h>
#include <adios_error.h>
#include <glob.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>
#include "core/bp_utils.hpp"
#include "core/xgcBlobExtractor.h"
  
const std::string pointsName = "points";
const std::string numPointsName = "numPoints";
const std::string cellsName = "cells";
const std::string numCellsName = "numCells";
const std::string cellShape = "triangle"; // FIXME
const std::string dpotName = "dpot";
const std::string psiName = "psi";
const std::string labelsName = "labels";
const std::string centering = "point";

  
typedef websocketpp::server<websocketpp::config::asio> server;
typedef server::message_ptr message_ptr;

using json = nlohmann::json;
using websocketpp::lib::placeholders::_1;
using websocketpp::lib::placeholders::_2;
using websocketpp::lib::bind;

server wsserver;

XGCBlobExtractor *ex = NULL;
std::mutex mutex_ex;

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
 
  if (incoming["type"] == "requestMesh") {
    fprintf(stderr, "requesting mesh!\n");
    outgoing["type"] = "mesh";
    outgoing["data"] = ex->jsonfyMesh();
  } else if (incoming["type"] == "requestData") {
    // int client_current_time_index = incoming["client_current_time_index"];

    outgoing["type"] = "data";

    mutex_ex.lock();
    const size_t timestep = ex->getTimestep();
    const int nnodes = ex->getNNodes();
    const double *ptr = ex->getData(); 
    std::vector<double> dpot(ptr, ptr + nnodes);
    std::vector<int> labels = ex->getLabels(0);
    mutex_ex.unlock();
   
    outgoing["timestep"] = timestep;
    outgoing["data"] = dpot;
    outgoing["labels"] = labels;
  } else if (incoming["type"] == "stopServer") {
    s->stop_listening();
    return;
  }

  try {
    // s->send(hdl, msg->get_payload(), msg->get_opcode()); // echo
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
    wsserver.set_access_channels(websocketpp::log::alevel::all);
    wsserver.clear_access_channels(websocketpp::log::alevel::frame_payload);

    // Initialize Asio
    wsserver.init_asio();

    // Register our message handler
    wsserver.set_message_handler(bind(&onMessage, &wsserver, ::_1, ::_2));

    // Listen on port 9002
    wsserver.listen(port);

    // Start the server accept loop
    wsserver.start_accept();
 
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

    // Start the ASIO io_service run loop
    wsserver.run();
  } catch (websocketpp::exception const & e) {
    std::cout << e.what() << std::endl;
  } catch (...) {
    std::cout << "other exception" << std::endl;
  }
}

void openUnstructuredMeshDataFile(MPI_Comm comm, const std::string& fileName, const std::string& writeMethod, const std::string& writeMethodParams, 
    int nNodes, int nTriangles, double *coords, int *conn, double *psi, int64_t &fileHandle, int64_t &dpotId, int64_t &labelsId)
{
  const std::string groupName = "xgc_blobs", meshName="xgc_mesh2D"; // , fileName="test.bp";
  int64_t groupHandle = -1;

  adios_init_noxml(comm);
  adios_declare_group(&groupHandle, groupName.c_str(), "", adios_stat_default);
  adios_select_method(groupHandle, writeMethod.c_str(), writeMethodParams.c_str(), "");
  adios_define_schema_version(groupHandle, (char*)"1.1");
  adios_define_mesh_timevarying("no", groupHandle, meshName.c_str());

  adios_delete_vardefs(groupHandle);
  adios_open(&fileHandle, groupName.c_str(), fileName.c_str(), "w", comm);
  
  adios_define_mesh_unstructured(
      (char*)pointsName.c_str(), 
      (char*)cellsName.c_str(), 
      (char*)numCellsName.c_str(), 
      (char*)cellShape.c_str(), 
      (char*)numPointsName.c_str(), 
      (char*)"2",
      groupHandle, 
      meshName.c_str());
  
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
  adios_write_byid(fileHandle, cellId, &conn[0]);
  
  // dpot
  adios_define_var_mesh(groupHandle, dpotName.c_str(), meshName.c_str());
  adios_define_var_centering(groupHandle, dpotName.c_str(), centering.c_str());
  dpotId = adios_define_var(
      groupHandle, 
      dpotName.c_str(), 
      "",
      adios_double, 
      std::to_string(nNodes).c_str(), 
      std::to_string(nNodes).c_str(),
      "0");
  
  // psi
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
 
  // labels
  adios_define_var_mesh(groupHandle, labelsName.c_str(), meshName.c_str());
  adios_define_var_centering(groupHandle, labelsName.c_str(), centering.c_str());
  labelsId = adios_define_var(
      groupHandle, 
      labelsName.c_str(), 
      "",
      adios_double, // adios_integer, 
      std::to_string(nNodes).c_str(), 
      std::to_string(nNodes).c_str(),
      "0");
  // adios_write_byid(fileHandle, labelsId, d_labels);
}

void writeUnstructredMeshData(int64_t fileHandle, int64_t dpotId, int64_t labelsId, int nNodes, int nTriangles, double *dpot, int *labels)
{
  // dpot
  adios_write_byid(fileHandle, dpotId, dpot);

  // labels
  if (labels != NULL) {
    double *d_labels = (double*)malloc(sizeof(double)*nNodes);
    for (int i=0; i<nNodes; i++) 
      d_labels[i] = static_cast<double>(labels[i]);

    adios_write_byid(fileHandle, labelsId, d_labels);
    free(d_labels);
  }

  // adios_close(fileHandle);
  // adios_finalize(0);
}

void writeUnstructredMeshDataFile_Legacy(int timestep, MPI_Comm comm, const std::string& fileName, const std::string& writeMethod, const std::string& writeMethodParams,
    int nNodes, int nTriangles, double *coords, int *conn_, double *dpot, double *psi, int *labels)
{
  const std::string groupName = "xgc_blobs", meshName="xgc_mesh2D"; // , fileName="test.bp";
  int64_t groupHandle = -1, fileHandle = -1;

  adios_init_noxml(comm);
  adios_declare_group(&groupHandle, groupName.c_str(), "", adios_stat_default);
  adios_select_method(groupHandle, writeMethod.c_str(), "", "");
  adios_define_schema_version(groupHandle, (char*)"1.1");
  adios_define_mesh_timevarying("no", groupHandle, meshName.c_str());

  adios_delete_vardefs(groupHandle);
  // adios_open(&fileHandle, groupName.c_str(), fileName.c_str(), "w", comm);
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
  adios_finalize(0);
}


int main(int argc, char **argv)
{
  std::thread *ws_thread = NULL;

  MPI_Init(&argc, &argv);
  
  using namespace boost::program_options;
  options_description opts(argv[0]);
  opts.add_options()
    ("mesh,m", value<std::string>(), "mesh_file")
    ("input,i", value<std::string>(), "input_file")
    ("pattern", value<std::string>(), "input file pattern, e.g. 'xgc.3d.*.bp'.  only works for BP.")
    ("output,o", value<std::string>()->default_value(""), "output_file")
    ("read_method,r", value<std::string>()->default_value("BP"), "read_method (BP|DATASPACES|DIMES|FLEXPATH)")
    ("write_method,w", value<std::string>()->default_value("MPI"), "write_method (POSIX|MPI|DIMES)")
    ("write_method_params", value<std::string>()->default_value(""), "write_method_params")
    ("skip", value<std::string>(), "skip timesteps that are specified in a json file")
    ("server,s", "enable websocket server")
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
  const std::string filename_output = vm["output"].as<std::string>();
  const std::string read_method_str = vm["read_method"].as<std::string>();
  const std::string write_method_str = vm["write_method"].as<std::string>();
  const std::string write_method_params_str = vm["write_method_params"].as<std::string>();
  std::string filename_input;
  bool single_input = false;
  if (vm.count("input")) {
    filename_input = vm["input"].as<std::string>();
    single_input = true;
  }

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
  int *conn;
  double *psi; 
  readTriangularMesh(meshFP, nNodes, nTriangles, &coords, &conn);
  readScalars<double>(meshFP, "psi", &psi);
  // adios_close(*meshFP);

#if 0
  int64_t fileHandleOutput, dpotIdOutput, labelsIdOutput;
  if (filename_output.size() > 0) 
    openUnstructuredMeshDataFile(MPI_COMM_WORLD, filename_output, write_method_str, write_method_params_str,
        nNodes, nTriangles, coords, conn, psi, fileHandleOutput, dpotIdOutput, labelsIdOutput);
#endif
  
  // starting server
  ex = new XGCBlobExtractor(nNodes, nTriangles, coords, conn);
  
  if (vm.count("server")) {
    ws_thread = new std::thread(startWebsocketServer, vm["port"].as<int>());
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
  // adios_read_bp_reset_dimension_order(varFP, 0);
    
  double *dpot = NULL;
  size_t current_time_index = 0; // only for multiple inputs.

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

    const int nPhi=1;
    // readValueInt(varFP, "nphi", &nPhi);
    // readScalars<double>(varFP, "dpot", &dpot);

    ADIOS_VARINFO *avi = adios_inq_var(varFP, "dpot");
    assert(avi != NULL);

    adios_inq_var_stat(varFP, avi, 0, 0);
    adios_inq_var_blockinfo(varFP, avi);
    adios_inq_var_meshinfo(varFP, avi);

    // uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {nPhi, static_cast<uint64_t>(nNodes), 1, 1};
    uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {static_cast<uint64_t>(nNodes), nPhi, 1, 1};
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
  
#if 0
    std::stringstream ssfilename;
    ssfilename << "original-" << current_time_index << ".bp";
    fprintf(stderr, "writing original data for test.\n");
    writeUnstructredMeshData(MPI_COMM_WORLD, ssfilename.str().c_str(), 
        write_method_str, nNodes, nTriangles, coords, conn, dpot, NULL, NULL); // psi, labels);
    fprintf(stderr, "original data written.\n");
    // FIXME
    if (read_method == ADIOS_READ_METHOD_BP) continue;
    else adios_advance_step(varFP, 0, 1.0);
    continue;
#endif

    mutex_ex.lock();
    ex->setData(current_time_index, nPhi, dpot);
    // ex->setPersistenceThreshold(persistence_threshold);
    // ex->buildContourTree3D();
    ex->buildContourTree2D(0);
    mutex_ex.unlock();

    int *labels = ex->getLabels(0).data();

    if (filename_output.length() > 0) {
      fprintf(stderr, "writing results for timestep %d\n", current_time_index); 
#if 0
      writeUnstructredMeshData(fileHandleOutput, dpotIdOutput, labelsIdOutput, 
          nNodes, nTriangles, dpot, labels); 
#endif
#if 1 // legacy single files
      // std::stringstream ss_filename;
      // ss_filename << filename_output << "." << std::setfill('0') << std::setw(5) << current_time_index << ".bp";
      writeUnstructredMeshDataFile_Legacy(current_time_index, MPI_COMM_WORLD, filename_output, write_method_str, write_method_params_str,
          nNodes, nTriangles, coords, conn, dpot, psi, labels);
#endif
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
    fprintf(stderr, "shutting down wsserver...\n");
    ws_thread->join();
  }
   
  // adios_close(*varFP);

  if (dpot != NULL) free(dpot);
  delete ex;
  ex = NULL;
  free(coords);
  free(conn);

  fprintf(stderr, "exiting...\n");
  MPI_Finalize();
  return 0;
}
