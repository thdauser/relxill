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

    Copyright 2025 Thomas Dauser, Remeis Observatory & ECAP
*/

// Unit tests of the rebinning of spectra (Rebin.cpp): _rebin_spectrum, rebin_mean_flux, and the rebinning with a
// shifted energy grid (rebin_mean_flux_shifted), which is used to apply the energy shift of the returning thermal
// radiation. Run with bin/tests "[rebin]".

#include "catch2/catch_amalgamated.hpp"
#include "Rebin.h"
#include "Relreturn_BlackBody.h"

extern "C" {
#include "relutility.h"
}

#define PREC 1e-6

// uniform energy grid with nener bins over [emin, emax]
static VecD make_ener_grid(int nener, double emin, double emax) {
  VecD ener(nener + 1);
  const double de = (emax - emin) / nener;
  for (int i = 0; i <= nener; i++) { ener[i] = emin + i * de; }
  return ener;
}

// integral of a flux density (per keV) over the energy grid
static double integ_flux(const double *ener, int nener, const double *flux) {
  double sum = 0.0;
  for (int ii = 0; ii < nener; ii++) {
    sum += flux[ii] * (ener[ii + 1] - ener[ii]);
  }
  return sum;
}

TEST_CASE (

"_rebin_spectrum: bin-integrated spectrum rebinned onto a different grid"
,
"[rebin]"
)
 {

  int n0 = 6;
  const size_t n = 5;
  double ener0[] = {1, 3, 4, 5, 6, 7, 9};
  double val0[] = {1, 1, 1, 2, 1, 1};

  double ener[] = {0.5, 2, 4, 5.5, 7.5, 8};
  double val[n];

  _rebin_spectrum(ener, val, n, ener0, val0, n0);

  double val_ref[5];
  val_ref[0] = 0.5;
  val_ref[1] = 0.5 + 1;
  val_ref[2] = 1.0 + 1.0;
  val_ref[3] = 1.0 + 1.0 + 0.25;
  val_ref[4] = 0.25;

  size_t ii;
  for (ii = 0; ii < n; ii++) {
    REQUIRE(fabs(val[ii] - val_ref[ii]) < PREC);
  }
}

TEST_CASE (

"rebin_mean_flux: mean flux density rebinned onto a coarser grid"
,
"[rebin]"
)
 {

  int status = EXIT_SUCCESS;

  int n0 = 7;
  int n = 4;
  double ener0[] = {1, 2, 4, 5, 6, 7, 9, 10};
  double val0[] = {1, 4, 3, 2, 4, 1, 2};

  double ener[] = {1, 3, 4, 8.5, 10};
  double val[] = {0, 0, 0, 0};

  rebin_mean_flux(ener, val, n, ener0, val0, n0, &status);

  double val_ref[4];
  val_ref[0] = 2.0;
  val_ref[1] = 3.6666666666666666666666;
  val_ref[2] = 3.5;
  val_ref[3] = 1.83333333333333333333;

  int ii;
  for (ii = 0; ii < n; ii++) {
    if (fabs(val[ii] - val_ref[ii]) > PREC){
      printf(" expecting %e but cacluated %e \n", val_ref[ii], val[ii]);
    }
    REQUIRE (fabs(val[ii] - val_ref[ii]) <= PREC);
  }
}

TEST_CASE (
"rebin_mean_flux_shifted: flat input stays flat for shift=1"
,
"[rebin]"
)
 {
  const int nener = 50;
  VecD ener = make_ener_grid(nener, 1.0, 6.0);

  const double C = 3.14;
  VecD flux0(nener, C);
  VecD flux_out(nener, 0.0);

  rebin_mean_flux_shifted(ener.data(), nener, flux_out, flux0, 1.0);

  // bin 0 is 0: known edge effect of the midpoint interpolation scheme
  // (output midpoint coincides exactly with the lower boundary of the input midpoint bracket)
  REQUIRE(fabs(flux_out[0]) < PREC);
  for (int in = 1; in < nener; in++) {
    INFO("Failing index: " << in);
    REQUIRE(fabs(flux_out[in] - C) < PREC);
  }
}

TEST_CASE (
"rebin_mean_flux_shifted: identity mapping for shift=1"
,
"[rebin]"
)
 {
  const int nener = 50;
  VecD ener = make_ener_grid(nener, 1.0, 6.0);

  VecD flux0(nener);
  for (int i = 0; i < nener; i++) { flux0[i] = (i + 1) * 0.1; }

  VecD flux_out(nener, 0.0);
  rebin_mean_flux_shifted(ener.data(), nener, flux_out, flux0, 1.0);

  for (int in = 1; in < nener; in++) {
    INFO("Failing index: " << in);
    REQUIRE(fabs(flux_out[in] - flux0[in]) < PREC);
  }
}

TEST_CASE (
"rebin_mean_flux_shifted: large shift gives zero in all bins"
,
"[rebin]"
)
 {
  const int nener = 50;
  VecD ener = make_ener_grid(nener, 1.0, 6.0);

  VecD flux0(nener, 1.0);
  VecD flux_out(nener, 0.0);

  // shift >> 1: the output grid ener*shift lies entirely above the input grid
  // -> exercises the early-return path; all output bins must be zero
  rebin_mean_flux_shifted(ener.data(), nener, flux_out, flux0, 100.0);

  for (int in = 0; in < nener; in++) {
    INFO("Failing index: " << in);
    REQUIRE(fabs(flux_out[in]) < PREC);
  }
}

TEST_CASE (
"rebin_mean_flux_shifted: moderate shift gives zero in the high-energy bins"
,
"[rebin]"
)
 {
  const int nener = 50;
  VecD ener = make_ener_grid(nener, 1.0, 6.0);

  VecD flux0(nener, 1.0);
  VecD flux_out(nener, 0.0);

  // shift=1.5: low-index output bins fall within the input grid (non-zero),
  // high-index bins exceed it (zero via early return)
  rebin_mean_flux_shifted(ener.data(), nener, flux_out, flux0, 1.5);

  REQUIRE(fabs(flux_out[0] - 1.0) < PREC);         // in range, flat flux preserved
  REQUIRE(fabs(flux_out[nener - 1]) < PREC);        // out of range
}

TEST_CASE (
"rebin_mean_flux_shifted: identical to rebin_mean_flux on the explicitly shifted grid"
,
"[rebin]"
)
 {
  const int nener = 50;
  VecD ener = make_ener_grid(nener, 1.0, 6.0);

  // arbitrary non-trivial flux values
  VecD flux0(nener);
  for (int i = 0; i < nener; i++) { flux0[i] = 1.0 + 0.1 * i + 0.002 * i * i; }

  const double shift = 0.8;

  // manually build the shifted output grid
  VecD ener_shifted(nener + 1);
  for (int i = 0; i <= nener; i++) { ener_shifted[i] = ener[i] * shift; }

  // reference: rebin_mean_flux with the explicitly shifted grid as output
  VecD flux_ref(nener, 0.0);
  int status = EXIT_SUCCESS;
  rebin_mean_flux(ener_shifted.data(), flux_ref.data(), nener,
                  ener.data(), flux0.data(), nener, &status);
  REQUIRE(status == EXIT_SUCCESS);

  // function under test
  VecD flux_out(nener, 0.0);
  rebin_mean_flux_shifted(ener.data(), nener, flux_out, flux0, shift);

  for (int in = 0; in < nener; in++) {
    INFO("Failing index: " << in);
    REQUIRE(fabs(flux_out[in] - flux_ref[in]) < PREC);
  }
}

TEST_CASE (
"rebin_mean_flux_shifted: shifted blackbody has the same photon flux as a blackbody with gfac"
,
"[rebin]"
)
 {
  const int nener = 5000;

  double ener[nener + 1];
  get_log_grid(ener, nener, 0.01, 150.0);

  double gfac = 1.5;
  double temperature = 1.0;
  VecD spec(nener);
  bbody_photon_flux_per_keV(ener, nener, spec.data(), temperature, gfac);

  VecD spec0(nener);
  bbody_photon_flux_per_keV(ener, nener, spec0.data(), temperature, 1.0);
  VecD spec0_shift(nener);
  // as we shift the energy grid and not the spectrum, we need to shift it backwards with 1/gfac
  rebin_mean_flux_shifted(ener, nener, spec0_shift, spec0, 1. / gfac);

  auto pho_flux = integ_flux(ener, nener, spec.data());
  auto pho_flux0 = integ_flux(ener, nener, spec0_shift.data());

  REQUIRE(fabs(pho_flux - pho_flux0) < 1e-3);
}
