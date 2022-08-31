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

#ifndef RELXILL_RELRETURN_C_RELRETURN_GRIDS_H_
#define RELXILL_RELRETURN_C_RELRETURN_GRIDS_H_

#include "Relreturn_Table.h"

typedef struct {

  double *rlo;
  double *rhi;
  int nrad;

  double* ener; // length is n_ener+1
  int n_ener;

  double **specRet;
  double **specPri;

} returnSpec2D;


void get_gfac_grid(double* gfac, double gmin, double gmax, int ng);
void get_std_bbody_energy_grid(int *n_ener, double **ener, int *status);

double **new_specZonesArr(int nener_inp, int nrad, int *status);
void sum_2Dspec(double *spec, double **spec_arr, int nener, int nrad, const int *status);

returnSpec2D *new_returnSpec2D(double *rlo, double *rhi, int nrad, double* ener, int n_ener, int *status);

returnSpec2D *getReturnradOutputStructure(const returningFractions *dat,
                                          double **spec_rr_zones,
                                          double **spec_prim_zones,
                                          double* ener, int n_ener,
                                          int *status);

#endif //RELXILL_RELRETURN_C_RELRETURN_GRIDS_H_
