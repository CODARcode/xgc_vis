#include <mpi.h>
#include <adios.h>
#include <adios_read.h>
#include <adios_error.h>
#include <cassert>
#include <iostream>
#include <boost/program_options.hpp>
#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>
#include "core/bp_utils.hpp"
#include "core/xgcBlobExtractor.h"
  
typedef websocketpp::server<websocketpp::config::asio> server;
typedef server::message_ptr message_ptr;

using json = nlohmann::json;
using websocketpp::lib::placeholders::_1;
using websocketpp::lib::placeholders::_2;
using websocketpp::lib::bind;

server wsserver;

XGCBlobExtractor *ex = NULL;
std::mutex mutex_ex;

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
    outgoing["type"] = "data";
    mutex_ex.lock();
    std::vector<double> dpot(ex->getData(), ex->getData() + ex->getNNodes());
    std::vector<int> labels = ex->getLabels(0);
    mutex_ex.unlock();
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

void writeUnstructredMeshData(MPI_Comm comm, const std::string& fileName, const std::string& writeMethod, // POSIX, MPI, or others
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
  adios_open(&fileHandle, groupName.c_str(), fileName.c_str(), "w", comm);

  // fprintf(stderr, "groupHandle=%lld, fileHandle=%lld\n", groupHandle, fileHandle);

  const std::string pointsName = "points";
  const std::string numPointsName = "numPoints";
  const std::string cellsName = "cells";
  const std::string numCellsName = "numCells";
  const std::string cellShape = "triangle"; // FIXME
  const std::string dpotName = "dpot";
  const std::string psiName = "psi";
  const std::string labelsName = "labels";
  const std::string centering = "point";

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
  adios_write_byid(fileHandle, cellId, &conn_[0]);

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
  int64_t labelsId = adios_define_var(
      groupHandle, 
      labelsName.c_str(), 
      "",
      adios_integer, 
      std::to_string(nNodes).c_str(), 
      std::to_string(nNodes).c_str(),
      "0");
  adios_write_byid(fileHandle, labelsId, labels);

  adios_close(fileHandle);
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
    ("output,o", value<std::string>()->default_value(""), "output_file")
    ("read_method,r", value<std::string>()->default_value("BP"), "read_method (BP|DATASPACES|DIMES|FLEXPATH)")
    ("write_method,w", value<std::string>()->default_value("MPI"), "write_method (POSIX|MPI)")
    ("ws,s", "enable websocket server")
    ("port,p", value<int>()->default_value(9002), "websocket server port")
    ("h", "display this information");
  
  positional_options_description posdesc;
  posdesc.add("mesh", 1);
  posdesc.add("input", 1);
  // posdesc.add("output", 1);
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opts)
      .positional(posdesc)
      .style(command_line_style::unix_style | command_line_style::allow_long_disguise)
      .run(), vm);
  notify(vm);

  if (!vm.count("mesh") || !vm.count("input") || vm.count("h")) {
    std::cout << opts << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  const std::string filename_mesh = vm["mesh"].as<std::string>();
  const std::string filename_input = vm["input"].as<std::string>();
  const std::string filename_output = vm["output"].as<std::string>();
  const std::string read_method_str = vm["read_method"].as<std::string>();
  const std::string write_method_str = vm["write_method"].as<std::string>();

  fprintf(stderr, "==========================================\n");
  fprintf(stderr, "filename_mesh=%s\n", filename_mesh.c_str());
  fprintf(stderr, "filename_input=%s\n", filename_input.c_str());
  fprintf(stderr, "filename_output=%s\n", filename_output.c_str());
  fprintf(stderr, "read_method=%s\n", read_method_str.c_str());
  fprintf(stderr, "write_method=%s\n", write_method_str.c_str());
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

  adios_read_init_method(read_method, MPI_COMM_WORLD, "");

  /// read mesh
  ADIOS_FILE *meshFP = NULL;
  while (1) {
    meshFP = adios_read_open_file(filename_mesh.c_str(), ADIOS_READ_METHOD_BP, MPI_COMM_WORLD); // always use ADIOS_READ_METHOD_BP for mesh
    if (meshFP == NULL) {
      fprintf(stderr, "failed to open mesh: %s, will retry in 1 second.\n", filename_mesh.c_str()); 
      sleep(1);
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
  
  // starting server
  ex = new XGCBlobExtractor(nNodes, nTriangles, coords, conn);
  
  if (vm.count("ws")) {
    ws_thread = new std::thread(startWebsocketServer, vm["port"].as<int>());
  }

  // read data
  ADIOS_FILE *varFP;
  if (read_method == ADIOS_READ_METHOD_BP)
    varFP = adios_read_open_file(filename_input.c_str(), read_method, MPI_COMM_WORLD);
  else 
    varFP = adios_read_open (filename_input.c_str(), read_method, MPI_COMM_WORLD, ADIOS_LOCKMODE_ALL, -1.0);
  adios_read_bp_reset_dimension_order(varFP, 0);

  while (adios_errno != err_end_of_stream) {
    fprintf(stderr, "reading data...\n");

    int nPhi;
    readValueInt(varFP, "nphi", &nPhi);

    double *dpot;
    readScalars<double>(varFP, "dpot", &dpot);

    fprintf(stderr, "starting analysis..\n");

    mutex_ex.lock();
    ex->setData(nPhi, dpot);
    // ex->setPersistenceThreshold(persistence_threshold);
    // ex->buildContourTree3D();
    ex->buildContourTree2D(0);
    mutex_ex.unlock();

    int *labels = ex->getLabels(0).data();

    if (filename_output.length() > 0) {
      fprintf(stderr, "dumping results..\n"); 
      writeUnstructredMeshData(MPI_COMM_WORLD, filename_output, write_method_str, nNodes, nTriangles, coords, conn, dpot, psi, labels);
#if 0
      ex->dumpMesh("xgc.mesh.json"); // only need to dump once
      ex->dumpBranchDecompositions("xgc.branches.json");
      ex->dumpLabels("xgc.labels.bin");
#endif
    }
    
    free(dpot);
    fprintf(stderr, "done.\n");
    if (read_method == ADIOS_READ_METHOD_BP) break;
    else 
      adios_advance_step(varFP, 0, 1.0);
  }
  
  if (ws_thread) {
    fprintf(stderr, "shutting down wsserver...\n");
    ws_thread->join();
  }

  delete ex;
  ex = NULL;
  free(coords);
  free(conn);

  fprintf(stderr, "exiting...\n");
  MPI_Finalize();
  return 0;
}
