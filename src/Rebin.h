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

#ifndef RELXILL__REBIN_H_
#define RELXILL__REBIN_H_

#ifdef __cplusplus
extern "C" {
#endif

/** linear interpolation in 1 dimension **/
double interp_lin_1d(double ifac_r, double rlo, double rhi);
double interp_lin_1d_float(double ifac_r, float rlo, float rhi);

/** log interpolation in 1 dimension **/
double interp_log_1d(double ifac_r, double rlo, double rhi);

/** linear interpolation in 2 dimensions **/
double interp_lin_2d(double ifac1, double ifac2, double r11, double r12, double r21, double r22);
double interp_lin_2d_float(double ifac1, double ifac2, float r11, float r12, float r21, float r22);

/** rebin spectrum to a given energy grid (length of ener is num_flux_bins+1) **/
void _rebin_spectrum(const double* ener, double* flu, int nbins, const double* ener0, const double* flu0, int nbins0);

/** for x0 descending and xn ascending, calculate the mean at xn from y0 **/
void inv_rebin_mean(const double* x0, double* y0, int n0, const double* xn, double* yn, int nn, int* status);

/** for x0 and xn ascending, calculate the mean at xn  (x0 has n0+1 bins and xn has nn+1 bins) **/
void rebin_mean_flux(double* x0, double* y0, int n0, double* xn, double* yn, int nn, int* status);

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus

#include "Types.h"

/** rebin flux0 (on grid ener) onto flux_out (on grid ener*shift) using midpoint interpolation **/
void rebin_mean_flux_shifted(double* ener, int nener, VecD& flux_out, const VecD& flux0, double shift);

#endif

#endif //RELXILL__REBIN_H_