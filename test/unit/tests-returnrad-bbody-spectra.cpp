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

// Unit tests of the returning thermal radiation (relxillBB, Relreturn_BlackBody.cpp): temperature profile of the
// radial zones, the primary and returning spectra of the zones (returnSpec2D), and the emissivity profile and
// ionization gradient derived from them. Run with bin/tests "[rrad-bbody][spectra]".

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "Relreturn_BlackBody.h"
#include "Relreturn_Table.h"
#include "Relphysics.h"
#include "Rellp.h"

extern "C" {
#include "relutility.h"
}

// spectra of the zones integrated over the (proper) area of the zones
static double *integSpecArea(double **spec, int n_ener, double spin, const returnSpec2D *returnSpec, int *status) {
  auto *areaIntegSpec = (double *) malloc(sizeof(double) * n_ener);
  CHECK_MALLOC_RET_STATUS(areaIntegSpec, status, areaIntegSpec)

  for (int jj = 0; jj < n_ener; jj++) {
    areaIntegSpec[jj] = 0.0;
  }
  for (int ii = 0; ii < returnSpec->nrad; ii++) {
    double prop_area = calc_proper_area_ring(returnSpec->rlo[ii], returnSpec->rhi[ii], spin);
    for (int jj = 0; jj < n_ener; jj++) {
      areaIntegSpec[jj] += spec[ii][jj] * prop_area;
    }
  }
  return areaIntegSpec;
}

static int checkTemperatureProfile(double *temperature, int n) {
  if (temperature[0] < 0.0) {
    printf("\n *** error: lowest bin of temperature profile is %f  and NOT >=0 as expected \n",
           temperature[0]);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

TEST_CASE (
"returning radiation: temperature of the innermost zone is non-negative for a = 0.85-0.9"
,
"[rrad-bbody][spectra]"
)
 {
  int status = EXIT_SUCCESS;

  std::vector<double> spinArray = {0.85, 0.851, 0.86, 0.87, 0.89, 0.9};

  double Rout = 1000;

  returningFractions *dat = nullptr;
  for (double spin: spinArray) {
    double Rin = kerr_rms(spin);
    dat = get_rrad_fractions(spin, Rin, Rout, &status);

    VecD temperature = get_temperature_profile_disk_zones(dat, Rin, 1.0, &status);

    REQUIRE(checkTemperatureProfile(temperature.data(), dat->nrad) == EXIT_SUCCESS);

    delete_returningFractions(&dat);
  }
}

TEST_CASE (
"returning radiation: primary and returning spectra of the zones (written to testrr-rframe-*.fits)"
,
"[rrad-bbody][spectra]"
)
 {
  int status = EXIT_SUCCESS;

  /* create an energy grid */
  int n_ener = 200;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 50.0);

  double Tin = 1.0;
  double spin = 0.998; // simply determines the radial grid here

  double Rin = kerr_rms(spin);
  double Rout = RMAX_RELRET;

  returnSpec2D *returnSpec = spec_returnrad_blackbody(ener, n_ener, Tin, Rin, Rout, spin, &status);
  REQUIRE(returnSpec != nullptr);

  double *photar_areaInteg = integSpecArea(returnSpec->specRet, n_ener, spin, returnSpec, &status);
  double *photar0_areaInteg = integSpecArea(returnSpec->specPri, n_ener, spin, returnSpec, &status);

  fits_rr_write_2Dspec("!testrr-rframe-rr-bbody.fits", returnSpec->specRet, ener, n_ener,
                       returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, &status);
  fits_rr_write_2Dspec("!testrr-rframe-prim-bbody.fits", returnSpec->specPri, ener, n_ener,
                       returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, &status);

  free(photar_areaInteg);
  free(photar0_areaInteg);

  delete_returnSpec2D(&returnSpec);
}

TEST_CASE (
"returning radiation: emissivity profile is positive and decreasing with radius"
,
"[rrad-bbody][spectra]"
)
 {
  int status = EXIT_SUCCESS;

  double spin = 0.998;
  double Rin = kerr_rms(spin);
  double Rout = RMAX_RELRET;
  double Tin = 2.0;

  int n_ener = 500;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 50.0);

  returnSpec2D *returnSpec = spec_returnrad_blackbody(ener, n_ener, Tin, Rin, Rout, spin, &status);
  REQUIRE(status == EXIT_SUCCESS);
  REQUIRE(returnSpec != nullptr);

  emisProfile *emis_profile =
      calculate_emis_profile_from_spectra(returnSpec->specRet,
                                          returnSpec->rlo,
                                          returnSpec->rhi,
                                          returnSpec->nrad,
                                          returnSpec->ener,
                                          returnSpec->n_ener,
                                          &status);
  REQUIRE(status == EXIT_SUCCESS);
  REQUIRE(emis_profile != nullptr);
  REQUIRE(emis_profile->nr == returnSpec->nrad);

  // profile should be in ascending radial order
  REQUIRE(emis_profile->re[0] < emis_profile->re[1]);

  // emissivity should be
  // (1) positive and non-zero for all radial zones
  // (2) decreasing as due to returning radiation
  for (int ii = 0; ii < emis_profile->nr; ii++) {
    REQUIRE(emis_profile->emis[ii] > 0.0);
    if (ii > 0) {
      REQUIRE(emis_profile->emis[ii] <= emis_profile->emis[ii - 1]);
    }
  }

  free_emisProfile(emis_profile);
  delete_returnSpec2D(&returnSpec);
}

TEST_CASE (
"returning radiation: maximum of the ionization gradient equals the logxi parameter"
,
"[rrad-bbody][spectra]"
)
 {
  int status = EXIT_SUCCESS;

  double spin = 0.998;
  double Rin = kerr_rms(spin);
  double Rout = RMAX_RELRET;
  double Tin = 2.0;

  int n_ener = 500;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 50.0);

  returnSpec2D *returnSpec = spec_returnrad_blackbody(ener, n_ener, Tin, Rin, Rout, spin, &status);
  REQUIRE(status == EXIT_SUCCESS);
  REQUIRE(returnSpec != nullptr);

  const double input_logxi = 3.0;
  LocalModel lmod(ModelName::relxillBB);
  lmod.set_par(XPar::logxi, input_logxi);

  auto ion_gradient = calc_ion_gradient_from_returnrad(returnSpec, lmod.get_model_params(), &status);
  REQUIRE(status == EXIT_SUCCESS);

  double lxi_max = ion_gradient.lxi[0];
  for (int ii = 1; ii < ion_gradient.nzones(); ii++) {
    if (ion_gradient.lxi[ii] > lxi_max) {
      lxi_max = ion_gradient.lxi[ii];
    }
  }

  REQUIRE(lxi_max == Catch::Approx(input_logxi).epsilon(1e-2));

  delete_returnSpec2D(&returnSpec);
}
