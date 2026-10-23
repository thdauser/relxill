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
#include "Relreturn_Datastruct.h"

#include <algorithm>
#include <string>

extern "C" {
#include "relutility.h"
}

#if !defined(RELXILL_SOURCE_DIR)
#warning: RELXILL_SOURCE_DIR is not defined
#endif

// ---- reference data of the returning radiation of the zones ---- //
// The returning radiation of the zones (specRet) for kT = 2 keV, the alpha-disk profile, the first order of the
// returning radiation, and a disk from the ISCO to RMAX_RELRET, on a logarithmic grid of 300 bins from 0.01 to 100
// keV. The reference files are written by the hidden test case "[.create-refdata]" (bin/tests "[.create-refdata]"),
// in the format of fits_rr_write_2Dspec (one row per zone; columns rlo, rhi, ener (lower bin edges), spec).
static constexpr int N_ENER_REFDATA_RRAD = 300;
static constexpr double KT_REFDATA_RRAD = 2.0; // default of relxillBB

static std::string get_fname_refdata_rrad(double spin) {
  char fname[100];
  snprintf(fname, sizeof(fname), "%s/test/unit/refdata/refvalues_rrad_bbody_specRet_a%.3f.fits", RELXILL_SOURCE_DIR,
           spin);
  return fname;
}

static void get_energy_grid_refdata_rrad(double *ener) {
  get_log_grid(ener, N_ENER_REFDATA_RRAD + 1, 0.01, 100.0);
}

// returning radiation of the zones for the parameters of the reference data
static returnSpec2D *calc_returnrad_refdata_params(double spin, double *ener, int *status) {
  setenv("RELXILL_RRAD_BBODY_TPROFILE", "1", 1);
  setenv("RELXILL_RRAD_BBODY_ORDER", "1", 1);
  returnSpec2D *return_spec =
      spec_returnrad_blackbody(ener, N_ENER_REFDATA_RRAD, KT_REFDATA_RRAD, kerr_rms(spin), RMAX_RELRET, spin, status);
  unsetenv("RELXILL_RRAD_BBODY_TPROFILE");
  unsetenv("RELXILL_RRAD_BBODY_ORDER");
  return return_spec;
}

// reads the spectra of the zones of a reference file and checks that its zones and its energy grid are the ones of
// the model (return_spec, on the grid ener)
static VecD_2D fits_read_refdata_rrad(const std::string &fname, const returnSpec2D *return_spec, const double *ener,
                                      int *status) {
  const int nrad = return_spec->nrad;
  VecD_2D spec(nrad, VecD(N_ENER_REFDATA_RRAD, 0.0));
  fitsfile *fptr = nullptr;
  if (fits_open_table(&fptr, fname.c_str(), READONLY, status)) {
    relxill_check_fits_error(status);
    return spec;
  }
  long nrows = 0;
  fits_get_num_rows(fptr, &nrows, status);
  if (*status == EXIT_SUCCESS && nrows != nrad) {
    printf(" *** error : reference file %s has %li zones instead of %i\n", fname.c_str(), nrows, nrad);
    *status = EXIT_FAILURE;
  }
  int anynul = 0;
  double nullval = 0.0;
  VecD rlo(nrad, 0.0);
  VecD ener_lo(N_ENER_REFDATA_RRAD, 0.0);
  if (*status == EXIT_SUCCESS) {
    fits_read_col(fptr, TDOUBLE, 1, 1, 1, nrad, &nullval, rlo.data(), &anynul, status);
    fits_read_col(fptr, TDOUBLE, 3, 1, 1, N_ENER_REFDATA_RRAD, &nullval, ener_lo.data(), &anynul, status);
    for (int ii = 0; ii < nrad; ii++) {
      fits_read_col(fptr, TDOUBLE, 4, ii + 1, 1, N_ENER_REFDATA_RRAD, &nullval, spec[ii].data(), &anynul, status);
    }
  }
  relxill_check_fits_error(status);
  int status_close = EXIT_SUCCESS;
  fits_close_file(fptr, &status_close);
  CHECK_STATUS_RET(*status, spec);

  for (int ii = 0; ii < nrad; ii++) {
    if (fabs(rlo[ii] / return_spec->rlo[ii] - 1) > 1e-10) {
      printf(" *** error : zone %i of %s starts at %e instead of %e\n", ii, fname.c_str(), rlo[ii], return_spec->rlo[ii]);
      *status = EXIT_FAILURE;
    }
  }
  for (int jj = 0; jj < N_ENER_REFDATA_RRAD; jj++) {
    if (fabs(ener_lo[jj] / ener[jj] - 1) > 1e-10) {
      printf(" *** error : energy grid of %s differs in bin %i\n", fname.c_str(), jj);
      *status = EXIT_FAILURE;
    }
  }
  return spec;
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

// Regression test of the returning radiation of the zones (specRet) against the reference data (see above), for
// a = 0.998 (relxillBB) and a = 0.5 (interpolation of the returning-radiation table in spin). The primary of the zones
// is tested independently in "disk emission (1b)" and "(3a)", and (3b) tests that specRet is calculated from this
// primary; this test detects any change of the returning radiation itself (returning fractions, interpolation of the
// table, energy shifts). In each zone, the spectrum must agree with the reference to 1e-6 of its maximum.
TEST_CASE (
"returning radiation: returning spectra of the zones agree with the reference data"
,
"[rrad-bbody][spectra][refdata]"
)
 {
  const double spin = GENERATE(0.998, 0.5);
  INFO("a = " << spin << ", reference file " << get_fname_refdata_rrad(spin));

  int status = EXIT_SUCCESS;
  double ener[N_ENER_REFDATA_RRAD + 1];
  get_energy_grid_refdata_rrad(ener);
  returnSpec2D *return_spec = calc_returnrad_refdata_params(spin, ener, &status);
  REQUIRE(status == EXIT_SUCCESS);
  const VecD_2D spec_ref = fits_read_refdata_rrad(get_fname_refdata_rrad(spin), return_spec, ener, &status);
  REQUIRE(status == EXIT_SUCCESS);

  double max_dev = 0.0;
  int izone_max = 0;
  for (int ii = 0; ii < return_spec->nrad; ii++) {
    const double max_ref = *std::max_element(spec_ref[ii].begin(), spec_ref[ii].end());
    REQUIRE(max_ref > 0);
    for (int jj = 0; jj < N_ENER_REFDATA_RRAD; jj++) {
      const double dev = fabs(return_spec->specRet[ii][jj] - spec_ref[ii][jj]) / max_ref;
      if (dev > max_dev) {
        max_dev = dev;
        izone_max = ii;
      }
    }
  }
  INFO("max. deviation from the reference data: " << max_dev << " (zone " << izone_max << ")");
  CHECK(max_dev < 1e-6);

  delete_returnSpec2D(&return_spec);
}

// writes the reference data of the test above (only run if called explicitly: bin/tests "[.create-refdata]")
TEST_CASE (
"returning radiation: create the reference data of the returning spectra of the zones"
,
"[.create-refdata]"
)
 {
  for (double spin: {0.998, 0.5}) {
    int status = EXIT_SUCCESS;
    double ener[N_ENER_REFDATA_RRAD + 1];
    get_energy_grid_refdata_rrad(ener);
    returnSpec2D *return_spec = calc_returnrad_refdata_params(spin, ener, &status);
    const std::string fname = "!" + get_fname_refdata_rrad(spin);
    fits_rr_write_2Dspec(fname.c_str(), return_spec->specRet, ener, N_ENER_REFDATA_RRAD, return_spec->rlo,
                         return_spec->rhi, return_spec->nrad, nullptr, &status);
    REQUIRE(status == EXIT_SUCCESS);
    printf(" written reference data: %s\n", fname.c_str() + 1);
    delete_returnSpec2D(&return_spec);
  }
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
