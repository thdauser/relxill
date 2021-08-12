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

    Copyright 2021 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "catch2/catch_amalgamated.hpp"

#include "common-functions.h"

extern "C" {
#include "relutility.h"
#include "relbase.h"
#include "relphysics.h"
}
#define LIMIT_PREC 1e-6


static double getRatioRelProfileAndRingArea(relline_spec_multizone *rprofile, int ii, double spin) {

  // proper area times the emissivity
  double propAreaRing = calc_proper_area_ring(rprofile->rgrid[ii], rprofile->rgrid[ii + 1], spin);
  double propAreaFullDisk = calc_proper_area_ring(rprofile->rgrid[0], rprofile->rgrid[rprofile->n_zones], spin);

  double ringAreaFraction = propAreaRing / propAreaFullDisk;

  double relFlux = calcSum(rprofile->flux[ii], rprofile->n_ener);

  return relFlux / ringAreaFraction;
}


TEST_CASE(" Relline Table Values ", "[basic]") {

  int status = EXIT_SUCCESS;
  relTable *tab = nullptr;

  const std::string filename_reltable = RELTABLE_FILENAME;
  INFO("loading rel table");
  read_relline_table(filename_reltable.c_str(), &tab, &status);

  // test certain values
  assert(tab != nullptr);
  assert(tab->arr != nullptr);

  double aref_val = 0.9829549;
  REQUIRE(fabs(tab->a[RELTABLE_NA - 2] - aref_val) < LIMIT_PREC);

  double mu0ref_val = 0.09476821;
  REQUIRE(fabs(tab->mu0[1] - mu0ref_val) < LIMIT_PREC);

  const int n = 5;
  const float ref_val[5] = {(float) 985.76074, (float) 0.01127052,
                            (float) 0.01120779, (float) 0.01121218, (float) 0.03022655};
  const float val[5] = {
      tab->arr[0][0]->r[1],
      tab->arr[0][0]->trff1[0][0],
      tab->arr[0][0]->trff1[0][1],
      tab->arr[0][0]->trff1[1][0],
      tab->arr[0][1]->trff1[1][0]
  };
  int ii;
  for (ii = 0; ii < n; ii++) {
    REQUIRE(fabsf(ref_val[ii] - val[ii]) < LIMIT_PREC);
  }

  free_relTable(tab);

}


TEST_CASE( "Lamp Post Table Values read from FITS", "[basic]") {

  /** test the currently implemented relline table
    * [current version used: rel_table_v0.4e]
    **/

  int status = EXIT_SUCCESS;
  lpTable *tab = nullptr;

  // load the table
  INFO(" loading relativistic table");
  std::string lp_filename = LPTABLE_FILENAME;
  read_lp_table(lp_filename.c_str(), &tab, &status);
  REQUIRE(status == EXIT_SUCCESS);

  // test certain values
  assert(tab != nullptr);
  assert(tab->dat != nullptr);
  assert(tab->a != nullptr);

  double aref_val = 0.98374581;
  REQUIRE (fabs(tab->a[LPTABLE_NA - 2] - aref_val) < LIMIT_PREC);

  double href_val = 1.8294963;
    REQUIRE(fabs(tab->dat[1]->h[1] - href_val) < LIMIT_PREC);

  const int n = 3;
  const float ref_val[3] = {(float) 7.618106e-05, (float) 2.6826601, (float) 1.2509402};  // the table value is -1.2509402, but we read all angles as positive values
  const float array_lp_table_values[3] = {
      tab->dat[1]->intens[2][1],
      tab->dat[1]->del[2][1],
      tab->dat[1]->del_inc[2][1]
  };

  int ii;
  for (ii = 0; ii < n; ii++) {
    REQUIRE(fabsf((ref_val[ii] - array_lp_table_values[ii]) / ref_val[ii]) < LIMIT_PREC);
  }

  // free memory
  free_lpTable(tab);

}


TEST_CASE("rebin spectrum", "[basic]") {

  int n0 = 6;
  const size_t n = 5;
  double ener0[] = {1, 3, 4, 5, 6, 7, 9};
  double val0[] = {1, 1, 1, 2, 1, 1};

  double ener[] = {0.5, 2, 4, 5.5, 7.5, 8};
  double val[n];

  rebin_spectrum(ener, val, n, ener0, val0, n0);

  double val_ref[5];
  val_ref[0] = 0.5;
  val_ref[1] = 0.5 + 1;
  val_ref[2] = 1.0 + 1.0;
  val_ref[3] = 1.0 + 1.0 + 0.25;
  val_ref[4] = 0.25;

  size_t ii;
  for (ii = 0; ii < n; ii++) {
    REQUIRE(fabs(val[ii] - val_ref[ii]) < LIMIT_PREC);
  }
}

TEST_CASE(" rebin mean flux ", "[basic]") {

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
    if (fabs(val[ii] - val_ref[ii]) > LIMIT_PREC){
      printf(" expecting %e but cacluated %e \n", val_ref[ii], val[ii]);
    }
    REQUIRE (fabs(val[ii] - val_ref[ii]) <= LIMIT_PREC);
  }
}


TEST_CASE(" 1d Interpolation Routines", "[basic]") {

  double ifac = 0.2;
  double rlo = 1.0;
  double rhi = 2.0;
  double val_1d = 1.2;
  REQUIRE(fabs(interp_lin_1d(ifac, rlo, rhi) - val_1d) < LIMIT_PREC);
}

TEST_CASE(" testing 2d interpolation ", "[basic]"){
  double ifac1 = 0.4;
  double ifac2 = 0.8;
  double val_2d = 2.52;

  double r11 = 1.0;
  double r12 = 2.0;
  double r21 = 2.0;
  double r22 = 4.0;
  INFO(" double 2d interpolation ");
  REQUIRE(fabs(interp_lin_2d(ifac1, ifac2, r11, r12, r21, r22) - val_2d) < LIMIT_PREC);

  auto rf11 = (float) 1.0;
  auto rf12 = (float) 2.0;
  auto rf21 = (float) 2.0;
  auto rf22 = (float) 4.0;
  INFO(" float 2d interpolation ");
  REQUIRE(fabs(interp_lin_2d_float(ifac1, ifac2, rf11, rf12, rf21, rf22) - val_2d) < LIMIT_PREC);
}


TEST_CASE("Relline Normalization Convergence", "[basic]"){
  /* Define criterion: for large radii the normalization of the line should
   * converge towards 1/4 of the disk area fraction of the respective radial ring
   * [definition of rel_table_v0.5a]
   *  - note that this depends on the normalization of the bkn power law emissivity profile
   */

  int status = EXIT_SUCCESS;

  int nzones = 100;
  const int min_number_zones = 10;

  relline_spec_multizone *rel_profile = nullptr;
  relParam *rel_param = nullptr;
  get_RelProfileConstEmisZones(&rel_profile, &rel_param, nzones, &status);

  double beginOuterRadii = 50.;
  int indexBeginOuterRadii = binary_search(rel_profile->rgrid, rel_profile->n_zones, beginOuterRadii);

  REQUIRE(rel_profile->n_zones - indexBeginOuterRadii > min_number_zones);

  double averageRatio = 0.0;
  double contributingFactorOfEachZone = 1.0 / ((double) (rel_profile->n_zones - indexBeginOuterRadii + 1));
  for (int ii = indexBeginOuterRadii; ii < rel_profile->n_zones; ii++) {
    averageRatio += getRatioRelProfileAndRingArea(rel_profile, ii, rel_param->a) * contributingFactorOfEachZone;
  }

  double referenceRatio = 0.5*cos(rel_param->incl);
  /* As we still expect some relat. effects testing for a high precision will lead to
 * errors, while the normalization is correct. Therefore we choose simply a weaker
 * criterion, which can be fulfilled considering GR.
 */
  double PREC = 0.1;
  REQUIRE(fabs(averageRatio - referenceRatio) < PREC);
  INFO("failed comparing the calculated ratio of the relline flux to the expectation");


}

TEST_CASE(" Normlization of the FFT Convolution", "[bbasic]") {

  int status = EXIT_SUCCESS;

  relline_spec_multizone *rel_profile = nullptr;
  double *xill_spec = nullptr;
  init_std_relXill_spec(&rel_profile, &xill_spec, &status);

  assert(rel_profile != nullptr);
  assert(xill_spec != nullptr);

  double relativistic_relline_profile_norm = calcSum(rel_profile->flux[0], rel_profile->n_ener);
  double sumProfileXill =
      calcSumInEnergyBand(xill_spec, rel_profile->n_ener, rel_profile->ener, EMIN_XILLVER, EMAX_XILLVER);

  REQUIRE( fabs(relativistic_relline_profile_norm - 1) < 1e-8);

  // requirements; spec_cache needs to be allocated for the FFT to work
  specCache *spec_cache = init_global_specCache(&status);
  REQUIRE(status == EXIT_SUCCESS);

  REQUIRE( rel_profile->n_zones == 1 );
  int izone = 0;

  double spec_conv_out[rel_profile->n_ener];
  convolveSpectrumFFTNormalized(rel_profile->ener, xill_spec, rel_profile->flux[izone], spec_conv_out,
                                rel_profile->n_ener, 1, 1, izone, spec_cache, &status);

  double sumProfileAfter =
      calcSumInEnergyBand(spec_conv_out, rel_profile->n_ener, rel_profile->ener, EMIN_XILLVER, EMAX_XILLVER);

  // now do the test (standard relline HAS to be normalized to 1)
  REQUIRE( fabs(sumProfileXill - sumProfileAfter) < 1e-8 );

  rel_profile->flux[izone][1000] = 1000.0;
  convolveSpectrumFFTNormalized(rel_profile->ener, xill_spec, rel_profile->flux[izone], spec_conv_out,
                                rel_profile->n_ener, 1, 1, izone, spec_cache, &status);

  double sumProfileAfterWrong =
      calcSumInEnergyBand(spec_conv_out, rel_profile->n_ener, rel_profile->ener, EMIN_XILLVER, EMAX_XILLVER);
  // now do the test (standard relline HAS to be normalized to 1)
  REQUIRE(fabs(sumProfileXill - sumProfileAfterWrong) > 1e-8);
  INFO(" the relativistic convolution is not normalized to 1 \n ");
  INFO(" but still the convoluation is normalized, which should not be the case\n ");

}



