#include <mpi.h>
#include <adios.h>
#include <adios_read.h>
#include <cassert>
#include <cfloat>
#include <chrono>
#include <queue>
#include <stack>
#include <QApplication>
#include "widget.h"
#include "bvh.h"
#include "volren.cuh"

extern "C"
{
extern void adios_read_bp_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran);
}

bool readValueInt(ADIOS_FILE *fp, const char *nm, int *val)
{
  ADIOS_VARINFO *avi = adios_inq_var(fp, nm);
  if (avi == NULL || avi->type != adios_integer) return false;
  *val = *((int*)avi->value);
  return true;
}


template <typename T> int readScalars(ADIOS_FILE *fp, const char *var, T **arr)
{
  ADIOS_VARINFO *avi = adios_inq_var(fp, var);
  if (avi == NULL) return false;

  adios_inq_var_stat(fp, avi, 0, 0);
  adios_inq_var_blockinfo(fp, avi);
  adios_inq_var_meshinfo(fp, avi);
    
  int nt = 1;
  uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {0, 0, 0, 0};
  
  for (int i = 0; i < avi->ndim; i++) {
    st[i] = 0;
    sz[i] = avi->dims[i];
    nt = nt * sz[i];
  }
  // fprintf(stderr, "%d, %d, %d, %d\n", sz[0], sz[1], sz[2], sz[3]);
  
  ADIOS_SELECTION *sel = adios_selection_boundingbox(avi->ndim, st, sz);
  assert(sel->type == ADIOS_SELECTION_BOUNDINGBOX);

  *arr = (T*)malloc(sizeof(double)*nt);

  adios_schedule_read_byid(fp, sel, avi->varid, 0, 1, *arr);
  int retval = adios_perform_reads(fp, 1);
  
  adios_selection_delete(sel);
  return avi->ndim;
}

bool readTriangularMesh(ADIOS_FILE *fp, int &nNodes, int &nTriangles, double **coords, int **conn, int **nextNode) 
{
  // read number of nodes and triangles
  readValueInt(fp, "n_n", &nNodes);
  readValueInt(fp, "n_t", &nTriangles);
  
  // read coordinates
  readScalars<double>(fp, "/coordinates/values", coords);

  // read triangles
  // readScalars<int>(fp, "/cell_set[0]/node_connect_list", conn);
  readScalars<int>(fp, "nd_connect_list", conn);

  // read nextNode
  readScalars<int>(fp, "nextnode", nextNode);
  
  fprintf(stderr, "mesh loaded: nNodes=%d, nTriangles=%d\n", 
      nNodes, nTriangles);

  return true;
}

void writeVTK(const int nNodes, const int nTriangles, double *coords, int *conn) 
{

}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  ADIOS_FILE *meshFP = adios_read_open_file("xgc.mesh.bp", ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  ADIOS_FILE *varFP = adios_read_open_file("xgc.3d.bp", ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
    
  adios_read_bp_reset_dimension_order(meshFP, 0);
  adios_read_bp_reset_dimension_order(varFP, 0);

  int nNodes, nTriangles, nPhi = 1;
  readValueInt(varFP, "nphi", &nPhi);

  fprintf(stderr, "reading mesh...\n");
  double *coords; 
  int *conn, *nextNode;
  readTriangularMesh(meshFP, nNodes, nTriangles, &coords, &conn, &nextNode);
  // fprintf(stderr, "nNodes=%d, nTriangles=%d, nPhi=%d\n", 
  //     nNodes, nTriangles, nPhi);

  buildBVHGPU(nNodes, nTriangles, coords, conn);

  // twisting mesh using nextNode
  fprintf(stderr, "%d, %d, %d\n", nextNode[conn[3374*3]], conn[3374*3+1], conn[3374*3+2]);
#if 0
  for (int i=0; i<nNodes*3; i++) {
    // conn[i] = nextNode[nextNode[nextNode[nextNode[nextNode[nextNode[nextNode[nextNode[nextNode[conn[i]]]]]]]]]];
    conn[i] = nextNode[conn[i]];
  }
#endif

#if 0
  const double p[2] = {10, 10}, p0[2] = {0, 0}, p1[2] = {10, 0}, p2[2] = {0, 10};
  fprintf(stderr, "awegawegawegaweg, %d\n", 
      insideTriangle(p, p0, p1, p2));
#endif

  fprintf(stderr, "reading data...\n");
  double *dpot;
  readScalars<double>(varFP, "dpot", &dpot);

  fprintf(stderr, "starting GUI...\n");
  QApplication app(argc, argv);
  QGLFormat fmt = QGLFormat::defaultFormat();
  fmt.setSampleBuffers(true);
  fmt.setSamples(16); 
  QGLFormat::setDefaultFormat(fmt); 
 
  const size_t nLimit = atoi(argv[1]);
  fprintf(stderr, "nLimit=%zu\n", nLimit);

  CGLWidget *widget = new CGLWidget;
  widget->setTriangularMesh(nNodes, nTriangles, nPhi, coords, conn);
  widget->setNLimit(nLimit);
  widget->setData(dpot);
  widget->show();
  app.exec();

  free(coords);
  free(conn);

  MPI_Finalize();
}
