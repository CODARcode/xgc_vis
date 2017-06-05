#include <mpi.h>
#include <adios.h>
#include <adios_read.h>
#include <cassert>

extern "C"
{
extern void adios_read_bp_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran);
}

bool readIntValue(ADIOS_FILE *fp, const char *nm, int *val)
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
  
  ADIOS_SELECTION *sel = adios_selection_boundingbox(avi->ndim, st, sz);
  assert(sel->type == ADIOS_SELECTION_BOUNDINGBOX);

  *arr = (T*)malloc(sizeof(double)*nt);

  adios_schedule_read_byid(fp, sel, avi->varid, 0, 1, *arr);
  int retval = adios_perform_reads(fp, 1);
  
  adios_selection_delete(sel);
  return avi->ndim;
}

bool readTriangularMesh(ADIOS_FILE *fp, int &nNodes, int &nTriangles, double **coords, int **conn) 
{
  // read number of nodes and triangles
  readIntValue(fp, "n_n", &nNodes);
  readIntValue(fp, "n_t", &nTriangles);
  
  // read coordinates
  readScalars<double>(fp, "/coordinates/values", coords);

  // read triangles
  readScalars<int>(fp, "/cell_set[0]/node_connect_list", conn);
  
  fprintf(stderr, "nNodes=%d, nTriangles=%d\n", 
      nNodes, nTriangles);

  return true;
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  ADIOS_FILE *meshFP = adios_read_open_file("xgc.mesh.bp", ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
  ADIOS_FILE *varFP = adios_read_open_file("xgc.3d.bp", ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
    
  adios_read_bp_reset_dimension_order(meshFP, 0);
  adios_read_bp_reset_dimension_order(varFP, 0);

  int nNodes, nTriangles, nPhi;
  readIntValue(varFP, "nphi", &nPhi);

  double *coords; 
  int *conn;
  readTriangularMesh(meshFP, nNodes, nTriangles, &coords, &conn);
  // fprintf(stderr, "nNodes=%d, nTriangles=%d, nPhi=%d\n", 
  //     nNodes, nTriangles, nPhi);

  free(coords);
  free(conn);

  MPI_Finalize();
}
