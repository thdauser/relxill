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

#ifndef XILLSPEC_H_
#define XILLSPEC_H_

extern "C" {
#include "xilltable.h"
#include "relutility.h"
#include "common.h"
}

#define EMIN_XILLVER_NORMALIZATION 0.1
#define EMAX_XILLVER_NORMALIZATION 1000.0
#define EMIN_XILLVER 0.01
#define EMAX_XILLVER EMAX_XILLVER_NORMALIZATION
#define N_ENER_COARSE 500


/** parameters for the convolution **/
#define N_ENER_CONV  4096  // number of bins for the convolution, not that it needs to follow 2^N because of the FFT
#define EMIN_RELXILL_CONV 0.00035  // minimal energy of the convolution (in keV)
#define EMAX_RELXILL_CONV 2000.0 // maximal energy of the convolution (in keV)

EnerGrid *get_relxill_conv_energy_grid();


double norm_factor_semi_infinite_slab(double incl_deg);
void norm_xillver_spec(xillSpec *spec, double incl);



xillSpec *get_xillver_spectra_table(xillTableParam *param, int *status);

xillSpec *get_xillver_spectra(xillParam *param, int *status);


xillSpec *new_xill_spec(int n_incl, int n_ener, int *status);
void free_xill_spec(xillSpec *spec);

void calc_xillver_angdep(double *xill_flux, xillSpec *xill_spec, const double *dist, const int *status);

void get_xillver_angdep_spec(double *o_xill_flux,
                             double *ener,
                             int n_ener,
                             double *rel_dist,
                             xillSpec *xill_spec,
                             int *status);

void getNormalizedXillverSpec(double* xill_flux, double* ener, int n_ener, xillParam* xill_param,
                              double *rel_cosne_dist, int *status);

void get_xillver_fluxcorrection_factors(const xillSpec *xill_spec,
                                        double *fac_fluxcorr,
                                        double *fac_gshift_fluxcorr,
                                        xillTableParam *xill_table_param,
                                        int *status);

void calc_primary_spectrum(double *pl_flux_xill, const double *ener, int n_ener,
                           const xillTableParam *xill_param, int *status,
                           double ener_shift_source_obs = 1.0);

double *calc_normalized_primary_spectrum(const double *ener, int n_ener,
                                         const relParam *rel_param, const xillTableParam *xill_param, int *status);

double calc_xillver_normalization_change(double energy_shift, const xillTableParam *xill_param_0);
double *calc_xillver_normalization_change_source_to_disk(const double *energy_shift,
                                                         int n_zones,
                                                         const xillTableParam *xill_param_0);

double calcNormWrtXillverTableSpec(const double *flux, const double *ener, int n, int *status);
EnerGrid *get_coarse_xillver_energrid();

#endif //XILLSPEC_H_
