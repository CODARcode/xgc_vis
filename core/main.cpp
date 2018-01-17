#include <mpi.h>
#include <adios.h>
#include <adios_read.h>
#include <adios_error.h>
#include <cassert>
#include <iostream>
#include <boost/program_options.hpp>
#include "core/bp_utils.hpp"
#include "core/xgcBlobExtractor.h"

void writeUnstructredMeshData(MPI_Comm comm, int nNodes, int nTriangles, double *coords, int *conn_, double *dpot, double *psi, int *labels)
{
  const std::string groupName = "xgc_blobs", meshName="xgc_mesh2D", fileName="test.bp";
  int64_t groupHandle = -1, fileHandle = -1;

  adios_init_noxml(comm);
  adios_declare_group(&groupHandle, groupName.c_str(), "", adios_stat_default);
  // adios_select_method(groupHandle, "MPI", "", "");
  adios_select_method(groupHandle, "POSIX", "", "");
  adios_define_schema_version(groupHandle, (char*)"1.1");
  adios_define_mesh_timevarying("no", groupHandle, meshName.c_str());

  adios_delete_vardefs(groupHandle);
  adios_open(&fileHandle, groupName.c_str(), fileName.c_str(), "w", comm);

  fprintf(stderr, "groupHandle=%lld, fileHandle=%lld\n", groupHandle, fileHandle);

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
  MPI_Init(&argc, &argv);
  
  using namespace boost::program_options;
  options_description opts(argv[0]);
  opts.add_options()
    ("mesh,m", value<std::string>(), "mesh_file")
    ("input,i", value<std::string>(), "input_file")
    ("output,o", value<std::string>(), "output_file")
    ("read_method,r", value<std::string>()->default_value("BP"), "read_method (BP|DATASPACES|DIMES|FLEXPATH)")
    ("write_method,w", value<std::string>()->default_value("BP"), "write_method (BP)")
    ("persistence_threshold,p", value<double>()->default_value(0), "persistence_threshold")
    ("h", "display this information");
  
  positional_options_description posdesc;
  posdesc.add("mesh", 1);
  posdesc.add("input", 1);
  posdesc.add("output", 1);
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opts)
      .positional(posdesc)
      .style(command_line_style::unix_style | command_line_style::allow_long_disguise)
      .run(), vm);
  notify(vm);

  if (!vm.count("mesh") || !vm.count("input") || !vm.count("output") || vm.count("h")) {
    std::cout << opts << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  const std::string filename_mesh = vm["mesh"].as<std::string>();
  const std::string filename_input = vm["input"].as<std::string>();
  const std::string filename_output = vm["output"].as<std::string>();
  const std::string read_method_str = vm["read_method"].as<std::string>();
  const std::string write_method_str = vm["write_method"].as<std::string>();
  const double persistence_threshold = vm["persistence_threshold"].as<double>();

  fprintf(stderr, "==========================================\n");
  fprintf(stderr, "filename_mesh=%s\n", filename_mesh.c_str());
  fprintf(stderr, "filename_input=%s\n", filename_input.c_str());
  fprintf(stderr, "filename_output=%s\n", filename_output.c_str());
  fprintf(stderr, "read_method=%s\n", read_method_str.c_str());
  fprintf(stderr, "write_method=%s\n", write_method_str.c_str());
  // fprintf(stderr, "persistence_threshold=%.03e\n", persistence_threshold);
  fprintf(stderr, "==========================================\n");

  int read_method;

  if (read_method_str == "BP") read_method = ADIOS_READ_METHOD_BP;
  else if (read_method_str == "DATASPACES") read_method = ADIOS_READ_METHOD_DATASPACES;
  else if (read_method_str == "DIMES") read_method = ADIOS_READ_METHOD_DIMES;
  else if (read_method_str == "FLEXPATH") read_method = ADIOS_READ_METHOD_FLEXPATH;
  else {
    fprintf(stderr, "unsupported read method: %s\n", read_method_str.c_str());
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  /// start
  ADIOS_FILE *meshFP = adios_read_open_file(filename_mesh.c_str(), ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  ADIOS_FILE *varFP = adios_read_open_file(filename_input.c_str(), ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
    
  adios_read_bp_reset_dimension_order(meshFP, 0);
  adios_read_bp_reset_dimension_order(varFP, 0);

  int nNodes, nTriangles, nPhi;
  readValueInt(varFP, "nphi", &nPhi);

  fprintf(stderr, "reading mesh...\n");
  double *coords;
  int *conn;
  double *psi; 
  readTriangularMesh(meshFP, nNodes, nTriangles, &coords, &conn);
  readScalars<double>(meshFP, "psi", &psi);
  // fprintf(stderr, "nNodes=%d, nTriangles=%d, nPhi=%d\n", 
  //     nNodes, nTriangles, nPhi);

  fprintf(stderr, "reading data...\n");
  double *dpot;
  readScalars<double>(varFP, "dpot", &dpot);

  // writeUnstructredMeshData(MPI_COMM_WORLD, nNodes, nTriangles, coords, conn, dpot);

  fprintf(stderr, "starting analysis..\n");
  XGCBlobExtractor *extractor = new XGCBlobExtractor(nNodes, nTriangles, nPhi, coords, conn);
  extractor->setData(dpot);
  // extractor->setPersistenceThreshold(persistence_threshold);
  // extractor->buildContourTree3D();
  extractor->buildContourTree2D(0);

  int *labels = extractor->getLabels(0).data();

  fprintf(stderr, "dumping results..\n"); // TODO
  writeUnstructredMeshData(MPI_COMM_WORLD, nNodes, nTriangles, coords, conn, dpot, psi, labels);
#if 0
  extractor->dumpMesh("xgc.mesh.json"); // only need to dump once
  extractor->dumpBranchDecompositions("xgc.branches.json");
  extractor->dumpLabels("xgc.labels.bin");
#endif

  delete extractor;
  free(dpot);
  free(coords);
  free(conn);

  MPI_Finalize();
  return 0;
}
