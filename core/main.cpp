#include "def.h"
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
#include <chrono>
#include <condition_variable>
#include <boost/program_options.hpp>
#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>
// #include <concurrentqueue.h>
#include "io/bp_utils.hpp"
#include "io/xgcMesh.h"
#include "io/xgcData.h"
#include "core/xgcBlobExtractor.h"
#include "volren/volrenEngine.h"
// #include "volren/bvh.h"
// #include "volren/volren.cuh"

typedef websocketpp::server<websocketpp::config::asio> server;
typedef server::message_ptr message_ptr;
server wss;

// using ConcurrentQueue = moodycamel::ConcurrentQueue;
using json = nlohmann::json;

XGCMesh xgcMesh;
XGCData xgcData;

VolrenEngine volrenEngine;

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

void onHttp(server *s, websocketpp::connection_hdl hdl) 
{
  server::connection_ptr con = s->get_con_from_hdl(hdl);
  con->append_header("Access-Control-Allow-Origin", "*");
  con->append_header("Access-Control-Allow-Headers", "*");
  
  std::string query = con->get_resource();
  // fprintf(stderr, "query=%s\n", query.c_str());

  if (query == "/requestMeshInfo") {
    con->set_body(xgcMesh.jsonfyMeshInfo().dump());
    con->set_status(websocketpp::http::status_code::ok);
  } else if (query == "/requestMesh") { 
    con->set_body(xgcMesh.jsonfyMesh().dump());
    con->set_status(websocketpp::http::status_code::ok);
  } else if (query == "/requestDataInfo") {
    con->set_body(xgcData.jsonfyDataInfo(xgcMesh).dump());
    con->set_status(websocketpp::http::status_code::ok);
  } else if (query == "/requestSingleSliceData") {
    con->set_body(xgcData.jsonfySingleSliceData(xgcMesh).dump());
    con->set_status(websocketpp::http::status_code::ok);
  } else if (query == "/requestData") {
    con->set_body(xgcData.jsonfyData(xgcMesh).dump());
    con->set_status(websocketpp::http::status_code::ok);
  } else if (query.find("/requestSampleAlongPsiContour?") == 0) {
    std::string str = query.substr(query.find("?")+1);
    // fprintf(stderr, "Len=%d, Strval=%s\n", strval.size(), strval.c_str());
    double isoval = std::stod(str);

    json j = xgcData.sampleAlongPsiContour(xgcMesh, isoval);
    con->set_body(j.dump());
    con->set_status(websocketpp::http::status_code::ok);
  } else if (query == "/requestSampleAlongPsiContourPolar") {
    json j = xgcData.sampleAlongPsiContourPolar(xgcMesh, 0.2);
    con->set_body(j.dump());
    con->set_status(websocketpp::http::status_code::ok);
  } else if (query == "/exitServer") {
    con->close(0, "exit");
    wss.stop_listening();
    volrenEngine.stop();
  } else if (query == "/testPost") {
    std::string posts = con->get_request_body();
    fprintf(stderr, "posts=%s\n", posts.c_str());
    con->set_body("hello world😱😱😱");
    con->set_status(websocketpp::http::status_code::ok);
  } else if (query == "/requestVolren") {
    con->append_header("Content-Type", "image/png");
    
    std::string posts = con->get_request_body();
    VolrenTask *task = volrenEngine.enqueueAndWait(posts);

    // response
    std::string response(task->png.buffer, task->png.size);
    con->set_body(response);
    con->set_status(websocketpp::http::status_code::ok);
    // con->defer_http_response();

    // clean
    delete task;
  } else {
    std::string response = "<html><body>404 not found</body></html>";
    con->set_body(response);
    con->set_status(websocketpp::http::status_code::not_found);
  }
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
#if 0
    binary = true;
    const int npx = viewport[2]*viewport[3];
    buffer_length = sizeof(int) + npx*3;
    buffer = (char*)malloc(buffer_length);
    int *type = (int*)(buffer);
    *type = 13;
    unsigned char *fb = (unsigned char*)(buffer+sizeof(int));
    memset(fb, 0, npx*3);
    for (int i=0; i<npx; i++) 
      fb[i*3] = 255;
    // memcpy(buffer, framebuf, buffer_length);
#endif
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
  using websocketpp::lib::placeholders::_1;
  using websocketpp::lib::placeholders::_2;
  using websocketpp::lib::bind;
  
  try {
    // Set logging settings
    wss.set_access_channels(websocketpp::log::alevel::all);
    wss.clear_access_channels(websocketpp::log::alevel::frame_payload);
    wss.set_open_handshake_timeout(0); // disable timer

    // Initialize Asio
    wss.init_asio();

    // Register our message handler
    wss.set_message_handler(bind(&onMessage, &wss, _1, _2));
    wss.set_http_handler(bind(&onHttp, &wss, _1));

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
    exit(EXIT_FAILURE);
  }
}

void sigint_handler(int)
{
  volrenEngine.stop();
  wss.stop_listening();
  exit(0);
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
  xgcMesh.readMeshFromADIOS(filename_mesh, ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);

  ex = new XGCBlobExtractor(xgcMesh.nNodes, xgcMesh.nTriangles, xgcMesh.coords, xgcMesh.conn);
  
  // starting server
  if (vm.count("server")) {
    ws_thread = new std::thread(startWebsocketServer, vm["port"].as<int>());
    signal(SIGINT, sigint_handler);
  }

  // read data
  fprintf(stderr, "opening data stream...\n");
  ADIOS_FILE *varFP;
  if (single_input > 0) {
    if (read_method == ADIOS_READ_METHOD_BP)
      varFP = adios_read_open_file(filename_input.c_str(), read_method, MPI_COMM_WORLD);
    else 
      varFP = adios_read_open (filename_input.c_str(), read_method, MPI_COMM_WORLD, ADIOS_LOCKMODE_ALL, -1.0);
    // fprintf(stderr, "varFP=%p, errno=%d\n", varFP, err_end_of_stream);
  }
  adios_read_bp_reset_dimension_order(varFP, 0);
    
  size_t current_time_index = 0; // only for multiple inputs.

  // output
  int64_t groupHandle = -1;
  if (!write_binary) {
    const std::string groupName = "xgc_blobs", meshName = "xgc_mesh2D";
    
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

    xgcMesh.nPhi = 1;
    readValueInt(varFP, "nphi", &xgcMesh.nPhi);
    // readScalars<double>(varFP, "dpot", &dpot);

    ADIOS_VARINFO *avi = adios_inq_var(varFP, "dpot");
    assert(avi != NULL);

    adios_inq_var_stat(varFP, avi, 0, 0);
    adios_inq_var_blockinfo(varFP, avi);
    adios_inq_var_meshinfo(varFP, avi);

    // uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {static_cast<uint64_t>(nNodes), static_cast<uint64_t>(nPhi), 1, 1};
    uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {static_cast<uint64_t>(xgcMesh.nPhi), static_cast<uint64_t>(xgcMesh.nNodes), 1, 1};
    ADIOS_SELECTION *sel = adios_selection_boundingbox(avi->ndim, st, sz);
    // fprintf(stderr, "%d, {%d, %d, %d, %d}\n", avi->ndim, avi->dims[0], avi->dims[1], avi->dims[2], avi->dims[3]);

    assert(sel->type == ADIOS_SELECTION_BOUNDINGBOX);

    if (xgcData.dpot == NULL)
      xgcData.dpot = (double*)malloc(sizeof(double)*xgcMesh.nPhi*xgcMesh.nNodes);

    adios_schedule_read_byid(varFP, sel, avi->varid, 0, 1, xgcData.dpot);
    adios_perform_reads(varFP, 1);
    adios_selection_delete(sel);

    xgcData.deriveSinglePrecisionDpot(xgcMesh);
    xgcData.deriveGradient(xgcMesh);


    if (skipped_timesteps.find(current_time_index) != skipped_timesteps.end()) {
      fprintf(stderr, "skipping timestep %lu.\n", current_time_index);
      adios_advance_step(varFP, 0, 1.0);
      continue;
    }

    fprintf(stderr, "starting analysis..\n");
  
    if (volren) {
      volrenEngine.start(xgcMesh, xgcData); 
      volrenEngine.enqueueAndWait("");
    }

    mutex_ex.lock();
    ex->setData(current_time_index, xgcMesh.nPhi, xgcData.dpot);
    // ex->setPersistenceThreshold(persistence_threshold);
    // ex->buildContourTree3D();
    // std::map<ctBranch*, size_t> branchSet = ex->buildContourTree2D(0);
    // ex->buildContourTree2DAll();
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
            xgcMesh.nNodes, xgcMesh.nTriangles, xgcMesh.coords, xgcMesh.conn, xgcData.dpot, xgcMesh.psi, labels);
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
    fprintf(stderr, "waiting for wss to exit...\n");
    ws_thread->join();
  }

  // adios_close(*varFP);
  delete ex;
  ex = NULL;

  adios_finalize(0);

  fprintf(stderr, "exiting...\n");
  MPI_Finalize();
  return 0;
}
