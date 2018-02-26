#ifndef _BP_UTILS_HPP
#define _BP_UTILS_HPP

#include <cassert>

extern "C"
{
extern void adios_read_bp_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran);
}

static bool readValueInt(ADIOS_FILE *fp, const char *nm, int *val)
{
  ADIOS_VARINFO *avi = adios_inq_var(fp, nm);
  if (avi == NULL || avi->type != adios_integer) return false;
  *val = *((int*)avi->value);
  return true;
}

template <typename T> static int readScalars(ADIOS_FILE *fp, const char *var, T **arr, bool allocate=true)
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

  if (allocate) 
    *arr = (T*)malloc(sizeof(double)*nt);

  adios_schedule_read_byid(fp, sel, avi->varid, 0, 1, *arr);
  int retval = adios_perform_reads(fp, 1);
  
  adios_selection_delete(sel);
  return avi->ndim;
}

static void writeUnstructredMeshDataFile(int timestep, MPI_Comm comm, int64_t groupHandle, const std::string& fileName, const std::string& writeMethod, const std::string& writeMethodParams,
    int nNodes, int nTriangles, double *coords, int *conn_, double *dpot, double *psi, int *labels)
{
  const std::string groupName = "xgc_blobs", meshName = "xgc_mesh2D";
  
  const std::string pointsName = "points";
  const std::string numPointsName = "numPoints";
  const std::string cellsName = "cells";
  const std::string numCellsName = "numCells";
  const std::string cellShape = "triangle"; // FIXME
  const std::string dpotName = "dpot";
  const std::string psiName = "psi";
  const std::string labelsName = "labels";
  const std::string centering = "point";

  int64_t fileHandle = -1;
  if (writeMethod == "POSIX" || writeMethod == "MPI")
    adios_open(&fileHandle, groupName.c_str(), fileName.c_str(), "w", comm);
  else 
    adios_open(&fileHandle, groupName.c_str(), fileName.c_str(), (timestep == 1 ? "w" : "a"), comm);

  // fprintf(stderr, "groupHandle=%lld, fileHandle=%lld\n", groupHandle, fileHandle);

  adios_define_mesh_unstructured(
      (char*)pointsName.c_str(), 
      (char*)cellsName.c_str(), 
      (char*)numCellsName.c_str(), 
      (char*)cellShape.c_str(), 
      (char*)numPointsName.c_str(), 
      (char*)"2",
      groupHandle, 
      meshName.c_str());
 
  // if (timestep == 1) {
  if (1) {
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
    
    // psi
    if (psi != NULL) {
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
    }
  }

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
  
  // labels
  if (labels != NULL) {
    double *d_labels = (double*)malloc(sizeof(double)*nNodes);
    for (int i=0; i<nNodes; i++) 
      d_labels[i] = static_cast<double>(labels[i]);

    adios_define_var_mesh(groupHandle, labelsName.c_str(), meshName.c_str());
    adios_define_var_centering(groupHandle, labelsName.c_str(), centering.c_str());
    int64_t labelsId = adios_define_var(
        groupHandle, 
        labelsName.c_str(), 
        "",
        adios_double, // adios_integer, 
        std::to_string(nNodes).c_str(), 
        std::to_string(nNodes).c_str(),
        "0");
    adios_write_byid(fileHandle, labelsId, d_labels);
    free(d_labels);
  }

  adios_close(fileHandle);
  // adios_finalize(0);
}

#endif
