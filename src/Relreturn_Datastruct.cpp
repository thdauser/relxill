/*
   This file is part of the RELXILL model code.

   RELXILL is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   RELXILL is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.
   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

    Copyright 2022 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "Relreturn_Datastruct.h"
#include "Relreturn_Table.h"

extern "C" {
#include "common.h"
#include "relutility.h"
}

/* using a universal energy grid for calculating a BBODY spectrum
 * this will save us many mallocs and loops
 * Function "get_std_bbody_energy_grid" will set up the grid
 */
#define N_BBODY_ENER 500
#define EMIN_BBODY 0.01
#define EMAX_BBODY 50
double *global_bbody_ener_std = nullptr;


returnSpec2D *getReturnradOutputStructure(const returningFractions *dat,
                                          double **spec_rr_zones,
                                          double **spec_prim_zones,
                                          double *ener,
                                          int n_ener,
                                          int *status) {
  returnSpec2D *returnSpec = new_returnSpec2D(dat->rlo, dat->rhi, dat->nrad, ener, n_ener, status);
  returnSpec->specRet = spec_rr_zones;
  returnSpec->specPri = spec_prim_zones;
  return returnSpec;
}


/* calculate g from gstar */
static double gstar2gfac(double gmin, double gmax, double gstar) {
  return gstar * (gmax - gmin) + gmin;
}

void get_gfac_grid(double* gfac, double gmin, double gmax, int ng) {
  for (int ii=0; ii<ng; ii++){
    gfac[ii] =  gstar2gfac(gmin,gmax, (ii+0.5)/ng);
  }
}

void get_std_bbody_energy_grid(int *n_ener, double **ener, int *status) {
  CHECK_STATUS_VOID(*status);
  if (global_bbody_ener_std == NULL) {
    global_bbody_ener_std = (double *) malloc((N_BBODY_ENER + 1) * sizeof(double));
    CHECK_MALLOC_VOID_STATUS(global_bbody_ener_std, status)
    get_log_grid(global_bbody_ener_std, (N_BBODY_ENER + 1), EMIN_BBODY, EMAX_BBODY);
  }
  (*n_ener) = N_BBODY_ENER;
  (*ener) = global_bbody_ener_std;
}


void sum_2Dspec(double *spec, double **spec_arr, int nener, int nrad, const int *status) {
  CHECK_STATUS_VOID(*status);

  for (int jj = 0; jj < nener; jj++) {
    spec[jj] = 0.0;
    for (int ii = 0; ii < nrad; ii++) {
      spec[jj] += spec_arr[ii][jj];
    }
  }
}

double **new_specZonesArr(int nener_inp, int nrad, int *status) {
  auto **spec_zones = (double **) malloc(sizeof(double *) * nrad);
  CHECK_MALLOC_RET_STATUS(spec_zones, status, spec_zones);
  for (int ii = 0; ii < nrad; ii++) {
    spec_zones[ii] = (double *) malloc(sizeof(double) * nener_inp);
    CHECK_MALLOC_RET_STATUS(spec_zones[ii], status, spec_zones);
  }
  return spec_zones;
}


returnSpec2D *new_returnSpec2D(double *rlo, double *rhi, int nrad, double* ener, int n_ener, int *status) {

  auto rspec = (returnSpec2D *) malloc(sizeof(returnSpec2D));
  CHECK_MALLOC_RET_STATUS(rspec, status, rspec)

  rspec->rlo = rlo;
  rspec->rhi = rhi;
  rspec->nrad = nrad;

  rspec->ener = ener;
  rspec->n_ener = n_ener;

  rspec->specPri = NULL;
  rspec->specRet = NULL;

  return rspec;
}

