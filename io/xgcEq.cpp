#include "xgcEq.h"
#include <fstream>
#include <cassert>

#if WITH_NC
#include <netcdf.h>
#endif

void XGCEq::writeToNetCDF(const std::string& filename) 
{
#if WITH_NC
  int ncid;
  int dimid_psi, dimid_z, dimid_r;
  int varid_psigrid, varid_I, varid_psi_rz, varid_rgrid, varid_zgrid;

  NC_SAFE_CALL( nc_create(filename.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &ncid) );
  
  // NC_SAFE_CALL( nc_def_dim(ncid, "psi", mpsi, &dimid_psi) );
  NC_SAFE_CALL( nc_def_dim(ncid, "z", mz, &dimid_z) );
  NC_SAFE_CALL( nc_def_dim(ncid, "r", mr, &dimid_r) );

  int dimids_rz[2] = {dimid_r, dimid_z};

  // NC_SAFE_CALL( nc_def_var(ncid, "psigrid", NC_DOUBLE, 1, &dimid_psi, &varid_psigrid) );
  // NC_SAFE_CALL( nc_def_var(ncid, "rgrid", NC_DOUBLE, 1, &dimid_r, &varid_rgrid) );
  // NC_SAFE_CALL( nc_def_var(ncid, "zgrid", NC_DOUBLE, 1, &dimid_z, &varid_zgrid) );
  // NC_SAFE_CALL( nc_def_var(ncid, "I", NC_DOUBLE, 1, &dimid_psi, &varid_I) );
  NC_SAFE_CALL( nc_def_var(ncid, "psi_rz", NC_DOUBLE, 2, dimids_rz, &varid_psi_rz) );
  NC_SAFE_CALL( nc_enddef(ncid) );

  size_t st[2] = {0, 0};
  size_t sz_psi[1] = {mpsi}, sz_z[1] = {mz}, sz_r[1] = {mr}, sz_rz[2] = {mr, mz};

  // NC_SAFE_CALL( nc_put_vara_double(ncid, varid_psigrid, st, sz_psi, psigrid.data()) );
  // NC_SAFE_CALL( nc_put_vara_double(ncid, varid_rgrid, st, sz_r, rgrid.data()) );
  // NC_SAFE_CALL( nc_put_vara_double(ncid, varid_zgrid, st, sz_z, zgrid.data()) );
  // NC_SAFE_CALL( nc_put_vara_double(ncid, varid_I, st, sz_psi, I.data()) );
  NC_SAFE_CALL( nc_put_vara_double(ncid, varid_psi_rz, st, sz_rz, psi_rz.data()) );

  NC_SAFE_CALL( nc_close(ncid) );
#endif
}

void XGCEq::parseFromFile(const std::string& filename)
{
  std::ifstream ifs(filename, std::ifstream::in);
  
  // ignore the first line
  ifs.ignore(1000, '\n');

  ifs >> mr >> mz >> mpsi;
  ifs >> min_r >> max_r >> min_z >> max_z;
  ifs >> axis_r >> axis_z >> axis_b;
  ifs >> x_psi >> x_r >> x_z;
  
  psigrid.resize(mpsi);
  for (int i=0; i<mpsi; i++) 
    ifs >> psigrid[i];

  I.resize(mpsi);
  for (int i=0; i<mpsi; i++) 
    ifs >> I[i];

  psi_rz.resize(mr*mz);
  for (int i=0; i<mr*mz; i++) 
    ifs >> psi_rz[i];

  int end_flag;
  ifs >> end_flag;
  assert(end_flag == -1);

  ifs.close();

  rgrid.resize(mr);
  for (int i=0; i<mr; i++)
    rgrid[i] = min_r + (max_r - min_r) / (mr - 1) * i;

  zgrid.resize(mz);
  for (int i=0; i<mz; i++) 
    zgrid[i] = min_z + (max_z - min_z) / (mz - 1) * i;

  printInfo();
}

void XGCEq::printInfo() const 
{
  fprintf(stderr, "mr=%zu, mz=%zu, mpsi=%zu\n", mr, mz, mpsi);
  fprintf(stderr, "min_r=%f, max_r=%f, min_z=%f, max_z=%f\n", 
      min_r, max_r, min_z, max_z);
  fprintf(stderr, "axis_r=%f, axis_z=%f, axis_b=%f\n", 
      axis_r, axis_z, axis_b);
  fprintf(stderr, "x_psi=%f, x_r=%f, x_z=%f\n", 
      x_psi, x_r, x_z);
}
