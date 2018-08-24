#include <QApplication>
#include <iostream>
#include "widget.h"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  QGLFormat fmt = QGLFormat::defaultFormat();
  fmt.setSampleBuffers(true);
  fmt.setSamples(16); 
  QGLFormat::setDefaultFormat(fmt); 
  
  CGLWidget *widget = new CGLWidget;
  widget->loadData("/Users/hguo/workspace/data/scalar-data");
  widget->show(); 

  app.exec();

  return 0;
}
