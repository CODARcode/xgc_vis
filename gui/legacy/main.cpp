#include <mpi.h>
#include <adios.h>
#include <adios_read.h>
#include <cassert>

#include <QApplication>
#include "widget.h"
#include "core/bp_utils.hpp"

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  ADIOS_FILE *meshFP = adios_read_open_file("xgc.mesh.bp", ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  ADIOS_FILE *varFP = adios_read_open_file("xgc.3d.bp", ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
    
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

  fprintf(stderr, "starting GUI...\n");
  QApplication app(argc, argv);
  QGLFormat fmt = QGLFormat::defaultFormat();
  fmt.setSampleBuffers(true);
  fmt.setSamples(16); 
  QGLFormat::setDefaultFormat(fmt); 
  
  CGLWidget *widget = new CGLWidget;
  widget->setTriangularMesh(nNodes, nTriangles, nPhi, coords, conn);
  widget->setData(dpot);
  widget->show();
  app.exec();

  free(coords);
  free(conn);

  MPI_Finalize();
}
