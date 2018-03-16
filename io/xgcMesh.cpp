#include "def.h"
#include <unistd.h>
#include <mpi.h>
#include <adios_error.h>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <list>
#include <iostream>
#include <map>
#include "common/base64.h"
#include "io/xgcMesh.h"
#include "io/bp_utils.hpp"

#if WITH_H5
#include <hdf5.h>
#endif

using json = nlohmann::json;

#if WITH_H5
void XGCMesh::readMeshFromH5(const std::string& filename)
{
  hid_t h5fid = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  hid_t h5id_nn = H5Dopen2(h5fid, "/n_n", H5P_DEFAULT);
  H5Dread(h5id_nn, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nNodes);
  H5Dclose(h5id_nn);
 
  hid_t h5id_nt = H5Dopen2(h5fid, "/n_t", H5P_DEFAULT);
  H5Dread(h5id_nt, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nTriangles);
  H5Dclose(h5id_nt);

  coords = (double*)malloc(sizeof(double)*nNodes*2);
  hid_t h5id_coords = H5Dopen2(h5fid, "/coordinates/values", H5P_DEFAULT);
  H5Dread(h5id_coords, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords);
  H5Dclose(h5id_coords);

  conn = (int*)malloc(sizeof(int)*nTriangles*3);
  hid_t h5id_conn = H5Dopen2(h5fid, "/cell_set[0]/node_connect_list", H5P_DEFAULT);
  H5Dread(h5id_conn, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, conn);
  H5Dclose(h5id_conn);

  nextNode = (int*)malloc(sizeof(int)*nNodes);
  hid_t h5id_nextnode = H5Dopen2(h5fid, "/nextnode", H5P_DEFAULT);
  H5Dread(h5id_nextnode, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, conn);
  H5Dclose(h5id_nextnode);

  psi = (double*)malloc(sizeof(double)*nNodes);
  hid_t h5id_psi = H5Dopen2(h5fid, "/psi", H5P_DEFAULT);
  H5Dread(h5id_psi, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, conn);
  H5Dclose(h5id_psi);

  H5Fclose(h5fid);

  preprocessMesh();
}
#endif

#if WITH_ADIOS
void XGCMesh::readMeshFromADIOS(const std::string& filename, ADIOS_READ_METHOD readMethod, MPI_Comm comm)
{
  adios_read_init_method(readMethod, comm, "");
  ADIOS_FILE *fp = NULL;
  while (1) {
    fp = adios_read_open_file(filename.c_str(), ADIOS_READ_METHOD_BP, comm); // always use ADIOS_READ_METHOD_BP for mesh
    if (fp == NULL) {
      fprintf(stderr, "failed to open mesh: %s, will retry in 5 seconds.\n", filename.c_str()); 
      sleep(5);
    } else break;
  }
  adios_read_bp_reset_dimension_order(fp, 0);
  
  readValueInt(fp, "n_n", &nNodes);
  readValueInt(fp, "n_t", &nTriangles);
  readScalars<double>(fp, "/coordinates/values", &coords);
  readScalars<int>(fp, "/cell_set[0]/node_connect_list", &conn);
  readScalars<int>(fp, "nextnode", &nextNode);
  readScalars<double>(fp, "psi", &psi);
 
  adios_read_finalize_method (ADIOS_READ_METHOD_BP);
  adios_read_close(fp);
  
  preprocessMesh();
}
#endif

void XGCMesh::preprocessMesh() {
  // psi
  psif = (float*)realloc(psif, sizeof(float)*nNodes);
  for (int i=0; i<nNodes; i++) {
    psif[i] = psi[i];
    if (psi_min > psif[i]) {
      psi_min = psif[i]; 
      psi_min_node = i;
    }
    // psi_min = std::min(psi_min, psif[i]);
    psi_max = std::max(psi_max, psif[i]);
  }
  psi_min_x = coords[psi_min_node*2];
  psi_min_y = coords[psi_min_node*2+1];
  
  // centroids and inverse of determinants
  centroidsf = (float*)realloc(centroidsf, sizeof(float)*nTriangles*2);
  invdetf = (float*)realloc(invdetf, sizeof(float)*nTriangles);
  for (int i=0; i<nTriangles; i++) {
    const int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    float x0 = coords[i0*2], x1 = coords[i1*2], x2 = coords[i2*2],
          y0 = coords[i0*2+1], y1 = coords[i1*2+1], y2 = coords[i2*2+1];
    float det = (y1-y2)*(x0-x2) + (x2-x1)*(y0-y2);
    invdetf[i] = 1.f / det;
    centroidsf[i*2] = (x0 + x1 + x2) / 3;
    centroidsf[i*2+1] = (y0 + y1 + y2) / 3;
  }
  
  // displacements
  dispf = (float*)realloc(dispf, sizeof(float)*nNodes*2);
  for (int i=0; i<nNodes; i++) {
    int j = nextNode[i];
    dispf[i*2] = coords[j*2] - coords[i*2];
    dispf[i*2+1] = coords[j*2+1] - coords[i*2+1];
    // fprintf(stderr, "%d, %d, %f, %f\n", i, j, dispf[i*2], dispf[i*2+1]);
  }

  // bounding box and centroid
  coords_min_x = coords_min_y = DBL_MAX;
  coords_max_x = coords_max_y = -DBL_MAX;
  coords_centroid_x = coords_centroid_y = 0;

  for (int i=0; i<nNodes; i++) {
    double x = coords[i*2], y = coords[i*2+1];
    coords_min_x = std::min(coords_min_x, x);
    coords_min_y = std::min(coords_min_y, y);
    coords_max_x = std::max(coords_max_x, x);
    coords_max_y = std::max(coords_max_y, y);
    coords_centroid_x += x; 
    coords_centroid_y += y;  
  }
  coords_centroid_x /= nNodes;
  coords_centroid_y /= nNodes;

  // others
  buildNeighbors();
  buildNodeGraph();
}

void XGCMesh::buildNeighbors()
{
  typedef std::tuple<int, int> edgeType;
  std::map<edgeType, int> edgeTriangleMap;

  auto makeEdge = [](int i0, int i1) {return std::make_tuple(i0, i1);};

  auto addEdge = [&makeEdge, &edgeTriangleMap](int triangleId, int i0, int i1) {
    // fprintf(stderr, "triangleId=%d, adding edge {%d, %d}\n", triangleId, i0, i1);
    edgeTriangleMap[ makeEdge(i0, i1) ] = triangleId;
  };

  for (int i=0; i<nTriangles; i++) {
    const int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    addEdge(i, i0, i1);
    addEdge(i, i1, i2);
    addEdge(i, i2, i0);
  }

  neighbors = (int*)realloc(neighbors, nTriangles*3*sizeof(int));

  for (int i=0; i<nTriangles; i++) {
    const int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];
    const edgeType edges[3] = {makeEdge(i1, i0), makeEdge(i2, i1), makeEdge(i0, i2)};

    for (int j=0; j<3; j++) {
      const edgeType& e = edges[j];
      if (edgeTriangleMap.find(e) != edgeTriangleMap.end()) neighbors[i*3+j] = edgeTriangleMap[e];
      else neighbors[i*3+j] = -1;
    }

    // fprintf(stderr, "triangleId=%d, neighbors={%d, %d, %d}\n", i, 
    //     neighbors[i*3], neighbors[i*3+1], neighbors[i*3+2]);
  }
}

void XGCMesh::buildNodeGraph()
{
  nodeGraph.clear();
  nodeGraph.resize(nNodes);
  for (int i=0; i<nTriangles; i++) {
    int i0 = conn[i*3], i1 = conn[i*3+1], i2 = conn[i*3+2];

    nodeGraph[i0].insert(i1);
    nodeGraph[i0].insert(i2);
    nodeGraph[i1].insert(i0);
    nodeGraph[i1].insert(i2);
    nodeGraph[i2].insert(i0);
    nodeGraph[i2].insert(i1);
  }
}


#if 0
std::vector<double> XGCMesh::testMarchingTriangles(double *scalar, double isoval)
{
  auto findZero = [scalar, isoval](int i0, int i1, double &alpha) {
    double f0 = scalar[i0], f1 = scalar[i1];
    alpha = (isoval - f0) / (f1 - f0);
    bool b = alpha >= 0 && alpha < 1;
    if (!b) alpha = std::nan("");
    return b;
  };

  std::vector<double> zeroPoints;

  auto findIntersection = [this, &findZero, &zeroPoints, scalar, isoval](int triangleId, int i[3]) {
    double a[3];
    bool b[3] = {findZero(i[0], i[1], a[0]), findZero(i[1], i[2], a[1]), findZero(i[2], i[0], a[2])};
    for (int j=0; j<3; j++) {
      if (b[j]) {
        int n0 = i[j], n1 = i[(j+1)%3];
        double alpha = a[j];
        double X = (1-alpha) * coords[n0*2] + alpha * coords[n1*2], 
               Y = (1-alpha) * coords[n0*2+1] + alpha * coords[n1*2+1];
        zeroPoints.push_back(X);
        zeroPoints.push_back(Y);
      }
    }
    if (b[0] || b[1] || b[2])
      fprintf(stderr, "triangleId=%d, alpha={%f, %f, %f}, neighbors={%d, %d, %d}\n", 
          triangleId, a[0], a[1], a[2], neighbors[triangleId*3], neighbors[triangleId*3+1], neighbors[triangleId*3+2]);
  };
  
  for (int i=0; i<nTriangles; i++)
    findIntersection(i, conn+i*3);

  return zeroPoints;
}

std::list<std::list<double> > XGCMesh::marchingTriangles(double *scalar, double isoval)
{
  auto findZero = [scalar, isoval](int i0, int i1, double &alpha) {
    double f0 = scalar[i0], f1 = scalar[i1];
    alpha = (isoval - f0) / (f1 - f0);
    bool b = alpha >= 0 && alpha < 1;
    if (!b) alpha = std::nan("");
    return b;
  };

  typedef std::tuple<double, double, double> TriangleIntersections;
  std::map<int, TriangleIntersections> candidateTriangles;

  auto findIntersection = [this, &findZero, &candidateTriangles, scalar, isoval](int triangleId, int i[3]) {
    double a[3];
    bool b[3] = {findZero(i[0], i[1], a[0]), findZero(i[1], i[2], a[1]), findZero(i[2], i[0], a[2])};
    if (b[0] || b[1] || b[2]) {
      candidateTriangles[triangleId] = std::make_tuple(a[0], a[1], a[2]);
      // fprintf(stderr, "triangleId=%d, alpha={%f, %f, %f}, neighbors={%d, %d, %d}\n", 
      //     triangleId, a[0], a[1], a[2], neighbors[triangleId*3], neighbors[triangleId*3+1], neighbors[triangleId*3+2]);
    }
  };

  for (int i=0; i<nTriangles; i++)
    findIntersection(i, conn+i*3);

  // traversal
  auto traceContourTrianglesOnSingleDirection = [this, &candidateTriangles](std::list<double>& contour, int seed, bool forward) {
    // fprintf(stderr, "seed=%d, direction=%d\n", seed, forward);
    int current = seed;
    while (candidateTriangles.find(current) != candidateTriangles.end()) {
      std::tuple<double, double, double> intersections = candidateTriangles[current];
      
      candidateTriangles.erase(current);
      
      // fprintf(stderr, "seed=%d, adding %d, neighbors={%d, %d, %d}, intersections={%f, %f, %f}\n", seed, current, 
      //     neighbors[current*3], neighbors[current*3+1], neighbors[current*3+2], 
      //     std::get<0>(intersections), std::get<1>(intersections), std::get<2>(intersections));
   
      // if (forward) contourTriangles.push_front(current);
      // else contourTriangles.push_back(current);
  
      int edge;
      double alpha;
      if ((!std::isnan(std::get<0>(intersections))) && candidateTriangles.find(neighbors[current*3]) != candidateTriangles.end()) {
        edge = 0;
        alpha = std::get<0>(intersections);
      } else if ((!std::isnan(std::get<1>(intersections))) && candidateTriangles.find(neighbors[current*3+1]) != candidateTriangles.end()) {
        edge = 1;
        alpha = std::get<1>(intersections);
      } else if ((!std::isnan(std::get<2>(intersections))) && candidateTriangles.find(neighbors[current*3+2]) != candidateTriangles.end()) {
        edge = 2;
        alpha = std::get<2>(intersections);
      } else break;
      
      int n0 = conn[current*3+edge], n1 = conn[current*3+(edge+1)%3];
      double X = (1-alpha) * coords[n0*2] + alpha * coords[n1*2], 
             Y = (1-alpha) * coords[n0*2+1] + alpha * coords[n1*2+1];
      // double X = alpha * coords[n0*2] + (1-alpha) * coords[n1*2], 
      //        Y = alpha * coords[n0*2+1] + (1-alpha) * coords[n1*2+1];
      if (forward) {contour.push_back(X); contour.push_back(Y);}
      else {contour.push_front(Y); contour.push_front(X);}
      // fprintf(stderr, "{%f, %f}\n", X, Y); //  TODO

      current = neighbors[current*3+edge];
    }
  };

  std::list<std::list<double> > contours;
  while (!candidateTriangles.empty()) {
    int seed = candidateTriangles.begin()->first;
    std::list<double> contour;
    traceContourTrianglesOnSingleDirection(contour, seed, true);
    traceContourTrianglesOnSingleDirection(contour, seed, false);
    contours.emplace_back(contour);
  }

  return contours;
}

std::vector<double> XGCMesh::sampleScalarsAlongPsiContour(double *scalar, int nSamples, double isoval)
{
  std::vector<double> results;
  std::list<double> contour_ = marchingTriangles(psi, isoval).front();
  std::vector<double> contour(std::begin(contour_), std::end(contour_));

  // for (int i=0; i<contour.size()/2; i++) {
  //   fprintf(stderr, "%f, %f\n", contour[i*2], contour[i*2+1]);
  // }

  return contour;
}
#endif


XGCMesh::~XGCMesh() {
  free(conn);
  free(psi);
  free(nextNode);
  free(coords);
  free(dispf);
  free(invdetf);
  free(neighbors);
}

json XGCMesh::jsonfyMeshInfo() const {
  json j;
  
  j["nPhi"] = nPhi;
  j["iPhi"] = iPhi;
  j["nNodes"] = nNodes;
  j["nTriangles"] = nTriangles;
  
  j["psi_min"] = psi_min;
  j["psi_max"] = psi_max;
  j["coords_min_x"] = coords_min_x;
  j["coords_min_y"] = coords_min_y;
  j["coords_max_x"] = coords_max_x;
  j["coords_max_y"] = coords_max_y;
  j["coords_centroid_x"] = coords_centroid_x;
  j["coords_centroid_y"] = coords_centroid_y;

  return j;
}

json XGCMesh::jsonfyMesh() const { 
  json j = jsonfyMeshInfo();

  j["coords"] = base64_encode((unsigned char*)coords, sizeof(double)*nNodes*2);
  // j["coords"] = std::vector<double>(coords, coords+nNodes*2);
  j["conn"] = base64_encode((unsigned char*)conn, sizeof(int)*nNodes*3);
  // j["conn"] = std::vector<int>(conn, conn+nNodes*3);
 
  return j;
}
