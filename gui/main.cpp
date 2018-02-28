#include <QApplication>
#include "widget.h"

int main(int argc, char **argv)
{
  XGCMesh m;
  m.readMeshFromADIOS(argv[1], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);

  ADIOS_FILE *varFP = adios_read_open_file(argv[2], ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  XGCData d;
  d.readDpotFromADIOS(m, varFP);

  m.buildNeighbors();
  m.buildNodeGraph();
  
  
  QApplication app(argc, argv);
  QGLFormat fmt = QGLFormat::defaultFormat();
  fmt.setSampleBuffers(true);
  fmt.setSamples(16); 
  QGLFormat::setDefaultFormat(fmt); 
  
  CGLWidget *widget = new CGLWidget(m, d);
#if 0
  widget->loadMeshFromJsonFile("xgc.mesh.json");
  widget->loadBranchesFromJsonFile("xgc.branches.json");
  widget->loadLabels("xgc.labels.bin"); // TODO: load labels from ADIOS
#endif
  widget->show();
  return app.exec();
}
