#ifndef _MU_CUH
#define _MU_CUH

typedef struct {
  int n_units; // typically 2 or 3, maximum 8
  int tscheme[8];  // ex. [7, 7, 2]
  int acc_tscheme[8]; // ex. [0, 7, 14]
  int sum_tscheme; // ex. 16
} mu_t;

__device__ __host__
inline int multiunits_lid2gid(const mu_t& mu, int unit, int lid) // return gid
{
  return (int)(lid/mu.tscheme[unit]) * mu.sum_tscheme + lid%mu.tscheme[unit] + mu.acc_tscheme[unit];
}

__device__ __host__
int multiunits_gid2lid(const mu_t& mu, int unit, int gid) // return lid
{
  return (int)(gid/mu.sum_tscheme) * mu.tscheme[unit] + (gid%mu.sum_tscheme - mu.acc_tscheme[unit]);
}

#endif
