#include <stdio.h>
#include <mpi.h>
#include <adios.h>
#include <adios_read.h>

#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkLongArray.h>
#include <vtkGenericCell.h>
#include <vtkDataSetWriter.h>
#include <vtkPointData.h>

#include <string>
#include <iostream>

extern "C"
{
extern void adios_read_bp_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran);
}


using namespace std;

static vtkDataArray *
AllocateTypedArray(ADIOS_VARINFO *avi)
{
  vtkDataArray *array = NULL;
  switch (avi->type)
  {
  case adios_unsigned_integer:
    array = vtkUnsignedIntArray::New(); 
    break;
  case adios_integer:
    array = vtkIntArray::New(); 
    break;
  case adios_unsigned_long:
    array = vtkUnsignedLongArray::New(); 
    break;
  case adios_long:
    array = vtkLongArray::New(); 
    break;
  case adios_real:
    array = vtkFloatArray::New(); 
    break;
  case adios_double:
    array = vtkDoubleArray::New(); 
    break;
  default:
    cout<<"Unknown Type: "<<avi->type<<endl;
    break;

  /*
  case adios_unsigned_byte:
    array = vtkCharArray::New();
    break;
  case adios_byte:
    array = vtkUnsignedCharArray::New();
    break;
  case adios_string:
    array = vtkCharArray::New();
    break;
  case adios_unsigned_short:
    array = vtkUnsignedShortArray::New();
    break;
  case adios_short:
    array = vtkShortArray::New();
    break;
  case adios_unsigned_integer:
    array = vtkUnsignedIntArray::New(); 
    break;
  case adios_integer:
    array = vtkIntArray::New(); 
    break;
  case adios_unsigned_long:
    array = vtkUnsignedLongArray::New(); 
    break;
  case adios_long:
    array = vtkLongArray::New(); 
    break;
  case adios_real:
    array = vtkFloatArray::New(); 
    break;
  case adios_double:
    array = vtkDoubleArray::New(); 
    break;
  case adios_complex:
    array = vtkFloatArray::New();
    break;
  case adios_double_complex:
    array = vtkDoubleArray::New();
    break;
      
  case adios_long_double: // 16 bytes
  default:
    std::string str = "Inavlid variable type";
    EXCEPTION1(InvalidVariableException, str);
    break;
  */
  }
  
  return array;
}


static vtkDataArray *
AllocateScalarArray(ADIOS_VARINFO *avi, ADIOS_SELECTION *sel)
{
  vtkDataArray *array = AllocateTypedArray(avi);
  int nt = 0;
  if (sel->type == ADIOS_SELECTION_BOUNDINGBOX)
  {
    nt = 1;
    for (int i = 0; i < sel->u.bb.ndim; i++)
      nt *= sel->u.bb.count[i];
  }
  if (avi->type == adios_complex || avi->type == adios_double_complex)
    array->SetNumberOfComponents(2);
  else
    array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(nt);
  fprintf(stderr, "%d, %d, %d\n", sel->type, ADIOS_SELECTION_BOUNDINGBOX, nt);
  return array;
}

static bool
ReadScalar(ADIOS_FILE *fp, const char *nm, int &val)
{
  ADIOS_VARINFO *avi = adios_inq_var(fp, nm);
  if (avi == NULL)
    return false;
  if (avi->type != adios_integer)
    return false;
  val = *((int*)avi->value);
  return true;
}

static bool
ReadScalarData(ADIOS_FILE *fp, const char *var, vtkDataArray **arr)
{
  ADIOS_VARINFO *avi = adios_inq_var(fp, var);
  if (avi == NULL)
    return false;
  adios_inq_var_stat(fp, avi, 0, 0);
  adios_inq_var_blockinfo(fp, avi);
  adios_inq_var_meshinfo(fp, avi);
  
  uint64_t start[4] = {0,0,0,0}, count[4] = {0,0,0,0};

  for (int i = 0; i < avi->ndim; i++)
    count[i] = avi->dims[i];
  
  ADIOS_SELECTION *sel = adios_selection_boundingbox(avi->ndim, start, count);
  *arr = AllocateScalarArray(avi, sel);
  if (*arr)
  {
    adios_schedule_read_byid(fp, sel, avi->varid, 0, 1, (*arr)->GetVoidPointer(0));
    int retval = adios_perform_reads(fp, 1);
  }
  adios_selection_delete(sel);
  return (*arr != NULL);
}



static vtkDataSet *
ReadMesh(ADIOS_FILE *fp, int numNodes, int numTris, int numPhi)
{
  int numPlanes = numPhi+1;
  numPlanes = numPhi;
  double dPhi = 2.0*M_PI/(double)(numPlanes-1);

  cout<<"NNodes = "<<numNodes<<endl;
  cout<<"NTris = "<<numTris<<endl;
  cout<<"NPlns = "<<numPlanes<<endl;    
  
  vtkDataArray *coords = NULL;
  if (!ReadScalarData(fp, "/coordinates/values", &coords))
    return NULL;
      
  vtkDataArray *conn = NULL;
  if (!ReadScalarData(fp, "/cell_set[0]/node_connect_list", &conn))
    return NULL;

  vtkUnstructuredGrid *grid = vtkUnstructuredGrid::New();

  //Create coordinates
  vtkPoints *pts = vtkPoints::New();
  pts->SetNumberOfPoints(numNodes*numPlanes);
  grid->SetPoints(pts);

  for (int i = 0; i < numPlanes; i++)
  {
    double phi = - (double)i * dPhi;
    for (int p = 0; p < numNodes; p++)
    {
      double R = coords->GetTuple1(p*2 +0);
      double Z = coords->GetTuple1(p*2 +1);
      pts->SetPoint(p+i*numNodes, R,phi,Z);
    }
  }

  for (int p=0; p<numNodes; p++) {
    double R = coords->GetTuple1(p*2 +0);
    double Z = coords->GetTuple1(p*2 +1);
    
    /// fprintf(stderr, "%d: (%f, %f)\n", p, R, Z);
  }
  
  coords->Delete();
  pts->Delete();

#if 0
  //Create wedges.
  int *connPtr = (int *)(conn->GetVoidPointer(0));
  vtkIdType wedge[6];
  for (int i = 0; i < numPlanes-1; i++)
  {
    for (int p = 0; p < numTris*3; p+=3)
    {
      int off = i*(numNodes);
      int p0 = connPtr[p+0];
      int p1 = connPtr[p+1];
      int p2 = connPtr[p+2];
      
      wedge[0] = p0 + off;
      wedge[1] = p1 + off;
      wedge[2] = p2 + off;

      off = (i+1)*(numNodes);
      wedge[3] = p0 + off;
      wedge[4] = p1 + off;
      wedge[5] = p2 + off;
      grid->InsertNextCell(VTK_WEDGE, 6, wedge);
    }
  }
  conn->Delete();
#endif
  
  //Create triangles.
  int *connPtr = (int *)(conn->GetVoidPointer(0));
  vtkIdType triangle[6]; 
  // for (int i = 0; i < numPlanes-1; i++)
  for (int i = 0; i < 1; i++)
  {
    for (int p = 0; p < numTris*3; p+=3)
    {
      int off = i*(numNodes);
      int p0 = connPtr[p+0];
      int p1 = connPtr[p+1];
      int p2 = connPtr[p+2];
      
      triangle[0] = p0 + off;
      triangle[1] = p1 + off;
      triangle[2] = p2 + off;

      /// fprintf(stderr, "%d, %d, %d\n", p0, p1, p2);

      grid->InsertNextCell(VTK_TRIANGLE, 3, triangle);
    }
  }
  conn->Delete();

  return grid;
}

vtkDataArray *
NextNodeMap(vtkDataArray *arr, ADIOS_FILE *meshFP, int numPhi, int numNodes, bool isXYZ)
{
  vtkDoubleArray *arrNN = vtkDoubleArray::New();
  int numPlanes = (isXYZ ? numPhi : numPhi+1);

  int nc = arr->GetNumberOfComponents();
  
  arrNN->SetNumberOfComponents(nc);
  arrNN->SetNumberOfTuples(numPlanes * numNodes);

  vtkDataArray *nextNode = NULL;
  ReadScalarData(meshFP, "nextnode", &nextNode);
  int *nnPtr = (int *)(nextNode->GetVoidPointer(0));    

  int idx = 0;
  for (int i = 0; i < numPhi; i++)
    for (int p = 0; p < numNodes; p++)
    {
      int off = i*(numNodes);
      int idx2 = (i==0 ? p : nnPtr[p]+off);
      idx2 = idx;
      
      for (int j = 0; j < nc; j++) {
          arrNN->SetComponent(idx, j, arr->GetComponent(idx2,j));
          // fprintf(stderr, "%d, %f\n", idx, arr->GetComponent(idx2,j));
      }
      idx++;
    }
  //For RPZ mesh, copy plane 1 into plane N
  if (!isXYZ)
  {
    for (int p = 0; p < numNodes; p++, idx++)
      for (int j = 0; j < nc; j++)
        arrNN->SetComponent(idx, j, arrNN->GetComponent(p,j));                
  }
  nextNode->Delete();

  return arrNN;
}

vtkDataArray *
ReadVar(ADIOS_FILE *fp, ADIOS_FILE *mfp, const char *nm, int numNodes, int numTris, int numPhi)
{
  vtkDataArray *var = NULL;
  if (!ReadScalarData(fp, nm, &var))
      return NULL;

  vtkDataArray *varNN = NextNodeMap(var, mfp, numPhi, numNodes, false);
  var->Delete();
  
  varNN->SetName(nm);
  return varNN;
}

int main (int argc, char ** argv) 
{
  int rank;
  MPI_Comm comm =  MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  
  string meshFname = "xgc.mesh.bp", varFname = "xgc.3d.bp";

  ADIOS_FILE *meshFP = adios_read_open_file(meshFname.c_str(), ADIOS_READ_METHOD_BP, comm);
  ADIOS_FILE *varFP = adios_read_open_file(varFname.c_str(), ADIOS_READ_METHOD_BP, comm);

  //Issue with how XGC writes data....
  adios_read_bp_reset_dimension_order(meshFP, 0);
  adios_read_bp_reset_dimension_order(varFP, 0);


  int numNodes, numTris, numPhi;
  ReadScalar(meshFP, "n_n", numNodes);
  ReadScalar(meshFP, "n_t", numTris);
  ReadScalar(varFP, "nphi", numPhi);
  
  vtkDataSet *grid = ReadMesh(meshFP, numNodes, numTris, numPhi);
  vtkDataArray *var = ReadVar(varFP, meshFP, "dpot", numNodes, numTris, numPhi);
  grid->GetPointData()->SetScalars(var);

  vtkDataSetWriter *wrt = vtkDataSetWriter::New();
  wrt->SetFileTypeToBinary();
  wrt->SetFileName("xgc.vtk");
  wrt->SetInputData(grid);
  wrt->Write();
  wrt->Delete();
  
  MPI_Finalize();
  return 0;
}

