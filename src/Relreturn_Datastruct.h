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
#include "Types.h"

typedef struct {

  VecD rlo;
  VecD rhi;
  int nrad;

  const double* ener; // length is n_ener+1
  int n_ener;

  /**
   * @brief specRet and specPri are 2D matrix used to store the returning and primary emitted
   * spectra for each radial zone given by the radial grid above
   *
   * The unit is bin integrated (Xspec standard) in photons/bin
   *
   * This double pointer represents an array of arrays, where each sub-array corresponds to
   * the spectral data within a particular radial zone. The number of sub-arrays is equal to
   * the number of radial divisions (`nrad`), and each sub-array's size corresponds to the
   * number of energy bins (`n_ener`).
   */
  double **specRet;
  double **specPri;

} returnSpec2D;


VecD get_gfac_grid(double gmin, double gmax, int ng);
void get_std_bbody_energy_grid(int *n_ener, double **ener, int *status);

double **new_specZonesArr(int nener_inp, int nrad, int *status);
VecU_D new_specZonesArrVec(int nener_inp, int nrad);
void sum_2Dspec(double *spec, double **spec_arr, int nener, int nrad, const int *status);

returnSpec2D *new_returnSpec2D(const VecD& rlo, const VecD& rhi, int nrad, const double* ener, int n_ener, int *status);
void delete_returnSpec2D(returnSpec2D** rspec);

returnSpec2D *getReturnradOutputStructure(const returningFractions *dat,
                                          double **spec_rr_zones,
                                          double **spec_prim_zones,
                                          double* ener, int n_ener,
                                          int *status);

void assert_double_array_isnormal(double* spec, int n);

#endif //RELXILL_RELRETURN_C_RELRETURN_GRIDS_H_
