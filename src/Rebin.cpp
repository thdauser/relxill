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

#include "Rebin.h"

extern "C" {
#include "common.h"
}

#include <cmath>
#include <cassert>


/*
 * C functions (interp / rebin)
 */

extern "C" {

double interp_lin_1d(double ifac_r, double rlo, double rhi) {
  return ifac_r * rhi + (1.0 - ifac_r) * rlo;
}

double interp_lin_1d_float(double ifac_r, float rlo, float rhi) {
  return ifac_r * rhi + (1.0 - ifac_r) * rlo;
}

double interp_log_1d(double ifac_r, double rlo, double rhi) {
  return exp(ifac_r * log(rhi) + (1.0 - ifac_r) * log(rlo));
}

double interp_lin_2d(double ifac1, double ifac2, double r11, double r12, double r21, double r22) {
  return (1.0 - ifac1) * (1.0 - ifac2) * r11 +
      (ifac1) * (1.0 - ifac2) * r12 +
      (1.0 - ifac1) * (ifac2) * r21 +
      (ifac1) * (ifac2) * r22;
}

double interp_lin_2d_float(double ifac1, double ifac2, float r11, float r12, float r21, float r22) {
  return (1.0 - ifac1) * (1.0 - ifac2) * r11 +
      (ifac1) * (1.0 - ifac2) * r12 +
      (1.0 - ifac1) * (ifac2) * r21 +
      (ifac1) * (ifac2) * r22;
}

void _rebin_spectrum(const double *ener, double *flu, int nbins, const double *ener0, const double *flu0, int nbins0) {

  int ii;
  int jj;
  int imin = 0;
  int imax = 0;

  for (ii = 0; ii < nbins; ii++) {

    flu[ii] = 0.0;

    /* check of the bin is outside the given energy range */
    if ((ener0[0] <= ener[ii + 1]) && (ener0[nbins0] >= ener[ii])) {

      /* need to make sure we are in the correct bin */
      while (ener0[imin] <= ener[ii] && imin <= nbins0) {
        imin++;
      }
      // need to set it back, as we just crossed to the next bin
      if (imin > 0) {
        imin--;
      }
      while ((ener0[imax] <= ener[ii + 1] && imax < nbins0)) {
        imax++;
      }
      if (imax > 0) {
        imax--;
      }

      double elo = ener[ii];
      double ehi = ener[ii + 1];
      if (elo < ener0[imin]) elo = ener0[imin];
      if (ehi > ener0[imax + 1]) ehi = ener0[imax + 1];

      if (imax == imin) {
        // flu0 needs to be interpolated
        flu[ii] = (ehi - elo) / (ener0[imin + 1] - ener0[imin]) * flu0[imin];
      } else {

        double dmin = (ener0[imin + 1] - elo) / (ener0[imin + 1] - ener0[imin]);
        double dmax = (ehi - ener0[imax]) / (ener0[imax + 1] - ener0[imax]);

        // same here in principle
        flu[ii] += flu0[imin] * dmin + flu0[imax] * dmax;

        for (jj = imin + 1; jj <= imax - 1; jj++) {
          flu[ii] += flu0[jj];
        }

      }

    }

  }
}

// for x0 descending and xn ascending, calculate the mean at xn
void inv_rebin_mean(const double *x0, double *y0, int n0, const double *xn, double *yn, int nn, int *status) {
  if (xn[0] > xn[nn - 1] || x0[nn - 1] > x0[0]) {
    RELXILL_ERROR(" *** grid in wrong order", status);
    return;
  }
  if (xn[0] < x0[n0 - 1] || xn[nn - 1] > x0[0]) {
    RELXILL_ERROR(" *** new grid is larger than the input grid", status);
    return;
  }

  int in = nn - 1; // traverse new array backwards
  int ii;
  for (ii = 0; ii < n0 - 1; ii++) {  // only go to the second to last bin

    while (in >= 0 && x0[ii] > xn[in] && x0[ii + 1] <= xn[in]) {
      double ifac_r = (xn[in] - x0[ii + 1]) / (x0[ii] - x0[ii + 1]);
      yn[in] = interp_lin_1d(ifac_r, y0[ii + 1], y0[ii]);

      in--;
    }
    if (in < 0) {
      break;
    }

  }

}

// for x0 and xn ascending, calculate the mean at xn  (x0 has n0+1 bins and xn has nn+1 bins)
void rebin_mean_flux(double *xn, double *yn, int nn, double *x0, double *y0, int n0, int *status) {
  if (xn[0] > xn[nn] || x0[0] > x0[n0]) {
    RELXILL_ERROR(" *** grid in wrong order", status);
    return;
  }
  /* if (xn[0]<x0[n0-1] || xn[nn-1]>x0[0]){
     RELXILL_ERROR(" *** new grid is larger than the input grid",status);
     return;
   } */

  int ii = 1; // we start at one, as we always need one bin lower than the xn[0] bin
  for (int in = 0; in < nn; in++) {  // only go to the second to last bin

    yn[in] = 0.0; // set it zero by default

    double xn_m = 0.5 * (xn[in] + xn[in + 1]);
    while (0.5 * (x0[ii - 1] + x0[ii]) < xn_m) {
      ii++;
      if (ii >= n0) {
        break;
      }
    }
    ii--;

    if (ii == 0) {  // at the lowest bin, we need to extrapolate
      ii = 1;
    }

    double x0_m_lo = 0.5 * (x0[ii - 1] + x0[ii]);
    double x0_m_hi = 0.5 * (x0[ii] + x0[ii + 1]);

    if (xn_m > x0_m_lo && xn_m <= x0_m_hi) {
      double ifac_r = (xn_m - x0_m_lo) / (x0_m_hi - x0_m_lo);
      yn[in] = interp_lin_1d(ifac_r, y0[ii - 1], y0[ii]);
    }

  }

}

} // extern "C"


/*
 * C++ functions
 */

void rebin_mean_flux_shifted(double* ener, int nener, VecD& flux_out, const VecD& flux0, double shift) {
    // Rebins flux0 (on grid ener[ii]) onto flux_out (on grid ener[ii]*shift).
    // Analogous to rebin_mean_flux: for each output bin, interpolate from the input
    // bin whose unshifted midpoint is closest.
    //
    // Performance: precompute the shifted output midpoint once per output bin,
    // so the while-loop body only computes unshifted input midpoints (no per-iteration multiply).

    int ii = 1;
    for (int in = 0; in < nener; in++) {
        flux_out[in] = 0.0;

        // midpoint of the output bin on the shifted grid (computed once per output bin)
        const double xn_0 = 0.5 * (ener[in] + ener[in + 1]);
        const double xn_m =  xn_0 * shift;

        // advance ii until the unshifted input midpoint at ii >= xn_m
        while (0.5 * (ener[ii - 1] + ener[ii]) < xn_m) {
            ii++;
            if (ii >= nener) { break; }
        }
        ii--;

        if (ii == 0) { ii = 1; }

        const double lo = 0.5 * (ener[ii - 1] + ener[ii]);
        const double hi = 0.5 * (ener[ii]     + ener[ii + 1]);

        if (xn_m > lo && xn_m <= hi) {
            flux_out[in] = interp_lin_1d((xn_m - lo) / (hi - lo), flux0[ii - 1], flux0[ii]);
        }
    }
}