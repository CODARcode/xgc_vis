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
#include "io/xgcDataReader.h"
#include "core/xgcBlobExtractor.h"
#include "core/xgcLevelSetAnalysis.h"
#include "volren/volrenEngine.h"
// #include "volren/bvh.h"
// #include "volren/volren.cuh"

#include <ftk/storage/storage.h>
#if WITH_ROCKSDB
#include <ftk/storage/rocksdbStorage.h>
#endif

#if WITH_QT
#include <QApplication>
#include "gui/widget.h"
#endif

typedef websocketpp::server<websocketpp::config::asio> server;
typedef server::message_ptr message_ptr;
server wss;

// using ConcurrentQueue = moodycamel::ConcurrentQueue;
using json = nlohmann::json;

XGCMesh xgcMesh;
XGCData xgcData;

VolrenEngine *volrenEngine = NULL;

XGCBlobExtractor *ex = NULL;
std::mutex mutex_ex;


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
    volrenEngine->stop();
  } else if (query == "/testPost") {
    std::string posts = con->get_request_body();
    fprintf(stderr, "posts=%s\n", posts.c_str());
    con->set_body("hello world😱😱😱");
    con->set_status(websocketpp::http::status_code::ok);
  } else if (query == "/requestVolren") {
    con->append_header("Content-Type", "image/png");
    
    std::string posts = con->get_request_body();
    VolrenTask *task = volrenEngine->enqueueAndWait(posts);

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
    // exit(EXIT_FAILURE);
  } catch (...) {
    std::cerr << "[wss] other exception" << std::endl;
    // exit(EXIT_FAILURE);
  }
}

void sigint_handler(int)
{
  volrenEngine->stop();
  wss.stop_listening();
  exit(0);
}

int main(int argc, char **argv)
{
  std::thread *ws_thread = NULL;

  // int np=1, rank=0;
  int np, rank;
  int required = MPI_THREAD_MULTIPLE, provided;
  // MPI_Init(&argc, &argv);
  MPI_Init_thread(&argc, &argv, required, &provided);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // assert(required == provided);

  using namespace boost::program_options;
  options_description opts(argv[0]);
  opts.add_options()
    ("mesh,m", value<std::string>(), "mesh_file")
    ("input_type,t", value<std::string>(), "ADIOS_STREAM|ADIOS_FILES|H5_FILES")
    ("input,i", value<std::string>(), "input")
    ("output,o", value<std::string>()->default_value(""), "output_file")
    ("db", value<std::string>()->default_value(""), "database name")
    ("output_prefix", value<std::string>()->default_value(""), "output file prefix (e.g. 'features' will generate features.%05d.bp)")
    ("output_branches,o", value<std::string>()->default_value(""), "output_branches")
    ("output_prefix_branches,o", value<std::string>()->default_value(""), "output_prefix_branches")
    ("read_method,r", value<std::string>()->default_value("BP"), "read_method (BP|DATASPACES|DIMES|FLEXPATH)")
    ("write_method,w", value<std::string>()->default_value("BIN"), "write_method (POSIX|MPI|DIMES|BIN), bin for raw binary")
    ("write_method_params", value<std::string>()->default_value(""), "write_method_params")
    ("server,s", "enable websocket server")
    ("volren,v", "enable volren (-s required)")
    ("gui,g", "start gui")
    ("port,p", value<int>()->default_value(9002), "websocket server port")
    ("help,h", "display this information");
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opts)
      .style(command_line_style::unix_style | command_line_style::allow_long_disguise)
      .run(), vm);
  notify(vm);

  if (vm.count("h") || !vm.count("mesh") || !(vm.count("input") || !vm.count("input_type"))) {
    if (rank == 0) std::cerr << opts << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  const std::string input_type_str = vm["input_type"].as<std::string>();
  const std::string input = vm["input"].as<std::string>();
  const std::string filename_mesh = vm["mesh"].as<std::string>();
  const std::string db_name = vm["db"].as<std::string>();
  const std::string output_prefix = vm["output_prefix"].as<std::string>();
  const std::string output_prefix_branches = vm["output_prefix_branches"].as<std::string>();
  const std::string read_method_str = vm["read_method"].as<std::string>();
  const std::string write_method_str = vm["write_method"].as<std::string>();
  const std::string write_method_params_str = vm["write_method_params"].as<std::string>();
  std::string filename_output = vm["output"].as<std::string>();
  std::string filename_output_branches = vm["output_branches"].as<std::string>();

  bool write_binary = (vm["write_method"].as<std::string>() == "BIN");
  bool volren = vm.count("volren") && vm.count("server");
  bool gui = vm.count("gui");

  if (rank == 0) {
    fprintf(stderr, "==========================================\n");
    fprintf(stderr, "filename_mesh=%s\n", filename_mesh.c_str());
    fprintf(stderr, "input_type=%s\n", input_type_str.c_str());
    fprintf(stderr, "input=%s\n", input.c_str());
    if (db_name.size() > 0)
      fprintf(stderr, "db_name=%s\n", db_name.c_str());
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
  }

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
  
  if (filename_mesh.rfind(".bp") != std::string::npos)
    xgcMesh.readMeshFromADIOS(filename_mesh, ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  else 
    xgcMesh.readMeshFromH5(filename_mesh);

  ex = new XGCBlobExtractor(xgcMesh.nNodes, xgcMesh.nTriangles, xgcMesh.coords, xgcMesh.conn);
  
  // starting server
  if (rank == 0 && vm.count("server")) {
    ws_thread = new std::thread(startWebsocketServer, vm["port"].as<int>());
    // signal(SIGINT, sigint_handler);
  }

  // read data
  fprintf(stderr, "opening input data\n");
  XGCDataReader *reader = XGCDataReaderFactory::newXGCDataReader(input_type_str, read_method_str, MPI_COMM_WORLD);
  if (reader == NULL) MPI_Abort(MPI_COMM_WORLD, 1);
  reader->open(input);

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

  std::shared_ptr<ftk::Storage> db;
  if (rank == 0 && db_name.size() > 0) {
#if WITH_ROCKSDB
    db = std::make_shared<ftk::RocksDBStorage>(ftk::RocksDBStorage());
#endif
    db->open(db_name);
  }

  if (gui) { // TODO: async
    QApplication app(argc, argv);
    QGLFormat fmt = QGLFormat::defaultFormat();
    fmt.setSampleBuffers(true);
    fmt.setSamples(16); 
    QGLFormat::setDefaultFormat(fmt); 
    
    CGLWidget *widget = new CGLWidget(xgcMesh, xgcData, *reader);
    widget->show(); 

    app.exec();
  } else {
    typedef std::vector<std::set<size_t> > Components;
    std::shared_ptr<Components> lastComponents;
    int lastTimestep = -1;

    ftk::Graph<> g;

    while (1) {
      int currentTimestep = reader->getCurrentTimestep();

      fprintf(stderr, "reading timestep %d\n", reader->getCurrentTimestep());
      reader->read(xgcMesh, xgcData);
      xgcData.deriveSinglePrecisionDpot(xgcMesh);
      xgcData.deriveGradient(xgcMesh);

      fprintf(stderr, "[rank=%d] starting analysis..\n", rank);

      // XGCLevelSetAnalysis::thresholdingByPercentageOfTotalEnergy(xgcMesh, xgcData, 0.6);
      std::shared_ptr<Components> components = 
        std::make_shared<Components>(XGCLevelSetAnalysis::extractSuperLevelSet2D(xgcMesh, xgcData, 220));

      if (db) db->put_obj("cc" + std::to_string(currentTimestep), *components);
   
      if (lastComponents) {
        fprintf(stderr, "associating...\n");
        auto mat = ftk::trackConnectedComponents(*lastComponents, *components);
        for (auto p : mat) {
          g.addNode(lastTimestep, p.first);
          g.addNode(currentTimestep, p.second);
          g.addEdge(lastTimestep, p.first, currentTimestep, p.second);
          // fprintf(stderr, "%zu -- %zu\n", pair.first, pair.second);
        }

        if (db) db->put_obj("mat" + std::to_string(lastTimestep) + "." + std::to_string(currentTimestep), mat);
      }

      lastTimestep = currentTimestep;
      lastComponents = components; 

      if (volren) {
        volrenEngine = new VolrenEngine();
        volrenEngine->start(MPI_COMM_WORLD, xgcMesh, xgcData);
        sleep(1);
        if (rank == 0) volrenEngine->enqueueAndWait("");
      }

      mutex_ex.lock();
      ex->setData(reader->getCurrentTimestep(), xgcMesh.nPhi, xgcData.dpot);
      // ex->setPersistenceThreshold(persistence_threshold);
      // ex->buildContourTree3D();
      // std::map<ctBranch*, size_t> branchSet = ex->buildContourTree2D(0);
      // ex->buildContourTree2DAll();
      mutex_ex.unlock();

      // write labels
      int *labels = ex->getLabels(0).data();
      if (output_prefix.length() > 0) {
        std::stringstream ss_filename;
        ss_filename << output_prefix << "." << std::setfill('0') << std::setw(5) << reader->getCurrentTimestep() 
          << (write_binary ? ".bin" : ".bp");
        filename_output = ss_filename.str();
      }
        
      if (filename_output.length() > 0) {
        fprintf(stderr, "writing results for timestep %d to %s\n", 
            reader->getCurrentTimestep(), filename_output.c_str());
        if (write_binary)
          ex->dumpLabels(filename_output);
        else 
          writeUnstructredMeshDataFile(reader->getCurrentTimestep(), MPI_COMM_WORLD, groupHandle, filename_output, write_method_str, write_method_params_str,
              xgcMesh.nNodes, xgcMesh.nTriangles, xgcMesh.coords, xgcMesh.conn, xgcData.dpot, xgcMesh.psi, labels);
      }
     
      // write branches
      if (output_prefix_branches.length() > 0) {
        std::stringstream ss_filename;
        ss_filename << output_prefix_branches << "." << std::setfill('0') << std::setw(5) << reader->getCurrentTimestep() << ".json";
        filename_output_branches = ss_filename.str();
      }

      if (filename_output_branches.length() > 0) {
        fprintf(stderr, "writing branch decompositions for timestep %d to %s\n", 
            reader->getCurrentTimestep(), filename_output_branches.c_str()); 
        // ex->dumpBranches(filename_output_branches, branchSet);
      }
      
      fprintf(stderr, "[rank=%d] done.\n", rank);

      if (reader->advanceTimestep() < 0) break; 
    } // while

    g.relabel();
    g.detectEvents();
    g.generateDotFile("mydot");
  }

#if 0
  if (filename_output.length() > 0) {
    adios_close(fileHandleOutput);
    adios_finalize(0);
  }
#endif
 
  if (ws_thread != NULL) {
    fprintf(stderr, "[rank=%d] waiting for wss to exit...\n", rank);
    ws_thread->join();
    delete ws_thread;
    ws_thread = NULL;
  }

  if (volrenEngine != NULL) {
    volrenEngine->stop();
    delete volrenEngine;
  }

  // adios_close(*varFP);
  delete ex;
  ex = NULL;

  delete reader;

  adios_finalize(0);

  if (db) db->close();

  fprintf(stderr, "[rank=%d] exiting...\n", rank);
  MPI_Finalize();
  return 0;
}
