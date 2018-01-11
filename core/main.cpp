#include <mpi.h>
#include <adios.h>
#include <adios_read.h>
#include <cassert>
#include <iostream>
#include <boost/program_options.hpp>
#include "core/bp_utils.hpp"
#include "core/xgcBlobExtractor.h"

int main(int argc, char **argv)
{
  using namespace boost::program_options;
  options_description opts(argv[0]);
  opts.add_options()
    ("mesh,m", value<std::string>(), "mesh_file")
    ("input,i", value<std::string>(), "input_file")
    ("output,o", value<std::string>(), "output_file")
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
    return 1;
  }

  const std::string filename_mesh = vm["mesh"].as<std::string>();
  const std::string filename_input = vm["input"].as<std::string>();
  const std::string filename_output = vm["output"].as<std::string>();

  fprintf(stderr, "filename_mesh=%s\n", filename_mesh.c_str());
  fprintf(stderr, "filename_input=%s\n", filename_input.c_str());
  fprintf(stderr, "filename_output=%s\n", filename_output.c_str());


  /// start
  MPI_Init(&argc, &argv);

  ADIOS_FILE *meshFP = adios_read_open_file(filename_mesh.c_str(), ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  ADIOS_FILE *varFP = adios_read_open_file(filename_input.c_str(), ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
    
  adios_read_bp_reset_dimension_order(meshFP, 0);
  adios_read_bp_reset_dimension_order(varFP, 0);

  int nNodes, nTriangles, nPhi;
  readValueInt(varFP, "nphi", &nPhi);

  fprintf(stderr, "reading mesh...\n");
  double *coords; 
  int *conn;
  readTriangularMesh(meshFP, nNodes, nTriangles, &coords, &conn);
  // fprintf(stderr, "nNodes=%d, nTriangles=%d, nPhi=%d\n", 
  //     nNodes, nTriangles, nPhi);

  fprintf(stderr, "reading data...\n");
  double *dpot;
  readScalars<double>(varFP, "dpot", &dpot);

  fprintf(stderr, "starting analysis..\n");
  XGCBlobExtractor *extractor = new XGCBlobExtractor(nNodes, nTriangles, nPhi, coords, conn);
  extractor->setData(dpot);
  extractor->buildContourTree3D();

  fprintf(stderr, "dumping results..\n"); // TODO
  extractor->dumpMesh("xgc.mesh.json"); // only need to dump once
  extractor->dumpBranchDecompositions("xgc.branches.json");
  extractor->dumpLabels("xgc.labels.bin");

  delete extractor;
  free(dpot);
  free(coords);
  free(conn);

  MPI_Finalize();
  return 0;
}
