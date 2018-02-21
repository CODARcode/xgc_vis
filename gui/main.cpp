#include <QApplication>
#include "widget.h"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  QGLFormat fmt = QGLFormat::defaultFormat();
  fmt.setSampleBuffers(true);
  fmt.setSamples(16); 
  QGLFormat::setDefaultFormat(fmt); 
  
  CGLWidget *widget = new CGLWidget;
#if 0
  widget->loadMeshFromJsonFile("xgc.mesh.json");
  widget->loadBranchesFromJsonFile("xgc.branches.json");
  widget->loadLabels("xgc.labels.bin"); // TODO: load labels from ADIOS
#endif
  widget->show();
  return app.exec();
}
