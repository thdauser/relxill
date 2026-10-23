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


// Unit tests of the table of the returning radiation fractions (Relreturn_Table.cpp), which is used by the returning
// radiation of the corona (relxilllp) and of the disk blackbody (relxillBB). Run with bin/tests "[rrad-table]".

#include "catch2/catch_amalgamated.hpp"
#include "Relreturn_Datastruct.h"
#include "Relreturn_Table.h"
#include "Relphysics.h"

extern "C" {
#include "relutility.h"
}

#define PREC 1e-6

static int is_grid_consistent(const double *rlo, const double *rhi, int nrad) {
  for (int ii = 1; ii < nrad; ii++) {
    if (fabs(rhi[ii - 1] - rlo[ii]) > PREC) {
      printf(" *** error: lower boundary rlo[ii]-rhi[ii-1] = %e - %e = %e  for ii=%i \n",
             rlo[ii],
             rhi[ii - 1],
             fabs(rhi[ii - 1] - rlo[ii]),
             ii);
      return 0;
    }
  }
  return 1;
}

static void test_table_radial_grid(tabulatedReturnFractions *dat) {
  REQUIRE(dat->nrad == RETURNRAD_TABLE_NR);
  REQUIRE(dat->ng == RETURNRAD_TABLE_NG);

  REQUIRE(dat->rlo[0] > 1.0);
  REQUIRE(dat->rhi[dat->nrad - 1] <= 1000.1); // don't care about border effects so add the 0.1

  REQUIRE(is_grid_consistent(dat->rlo, dat->rhi, dat->nrad) == 1);
}

static int test_table_fracE_norm(const tabulatedReturnFractions *dat) {
  int status = EXIT_SUCCESS;
  double kSumfrac;
  double kSumfrac_ref = 1.0;

  for (int jj = 0; jj < dat->nrad; jj++) {
    kSumfrac = 0.0;
    for (int ii = 0; ii < dat->nrad; ii++) {
      kSumfrac += dat->frac_e[ii][jj];
    }
    if (fabs(kSumfrac - kSumfrac_ref) > PREC) {
      RELXILL_ERROR("testing the normalization of FRAC_E failed (return radiation table)", &status);
      printf(" expecting a normalization of %e, but found %e\n", kSumfrac_ref, kSumfrac);
      break;
    }
  }
  return status;
}

static int test_table_fracG_norm(const tabulatedReturnFractions *dat) {
  int status = EXIT_SUCCESS;
  double kSumfrac;
  double kSumfrac_ref = 1.0;

  for (int ii = 0; ii < dat->nrad; ii++) {
    for (int jj = 0; jj < dat->nrad; jj++) {
      kSumfrac = 0.0;
      for (int kk = 0; kk < dat->ng; kk++) {
        kSumfrac += dat->frac_g[ii][jj][kk];
      }

      if (fabs(kSumfrac - kSumfrac_ref) > PREC) {
        RELXILL_ERROR("testing the normalization of FRAC_G failed (return radiation table)", &status);
        printf(" [%02i][%02i] expecting a normalization of %e, but found %e\n", ii, jj, kSumfrac_ref, kSumfrac);
        break;
      }
    }
  }
  return status;
}

static void test_table_fractions_normalization(tabulatedReturnFractions *dat) {
  REQUIRE(dat->frac_e[0][0] != 0);
  //  REQUIRE(test_table_fracE_norm(dat) == EXIT_SUCCESS);
  REQUIRE(test_table_fracG_norm(dat) == EXIT_SUCCESS);
}

TEST_CASE (

"returning radiation table: radial grid is consistent and the fractions are normalized"
,
"[rrad-table]"
)
 {

  int status = EXIT_SUCCESS;
  returnTable *tab = get_returnrad_table(&status);
  REQUIRE(tab != NULL);
  REQUIRE(status == EXIT_SUCCESS);

  DYNAMIC_SECTION("testing radial grid") {
    for (int ii = 0; ii < tab->nspin; ii++) {
      test_table_radial_grid(tab->retFrac[ii]);
      REQUIRE(status == EXIT_SUCCESS);
    }
  }

  DYNAMIC_SECTION("testing fractions normalization") {
    for (int ii = 0; ii < tab->nspin; ii++) {
      test_table_fractions_normalization(tab->retFrac[ii]);
    }
  }

}

TEST_CASE (

"returning fractions: number of radial zones decreases if Rin is increased"
,
"[rrad-table]"
)
 {

  int status = EXIT_SUCCESS;

  double spin = 0.9;
  double Rout = 1000;

  double Rin = kerr_rms(spin);
  returningFractions *dat = get_rrad_fractions(spin, Rin, Rout, &status);
  REQUIRE(status==EXIT_SUCCESS);
  int nrad_rms = dat->nrad;

  delete_returningFractions(&dat);

  Rin *= 2;
  dat = get_rrad_fractions(spin, Rin, Rout, &status);
  REQUIRE(status==EXIT_SUCCESS);
  int nrad_rfac2 = dat->nrad;

  delete_returningFractions(&dat);

  REQUIRE(nrad_rms > nrad_rfac2);

}

TEST_CASE (

"returning fractions: increasing Rin slightly only changes the fractions emitted by the innermost zone"
,
"[rrad-table]"
)
 {

  int status = EXIT_SUCCESS;
  double spin = 0.99;
  double rin = kerr_rms(spin);
  double rin_grid = rin * 1.001;
  double rout = 1000;

  returningFractions *rf0 = get_rrad_fractions(spin, rin, rout, &status);
  returningFractions *rf_grid = get_rrad_fractions(spin, rin_grid, rout, &status);

  REQUIRE(rf0 != NULL);
  REQUIRE(rf_grid != NULL);

  INFO("changing Rin slightly has to reduce the returning incident fraction onto the next zone ");
  REQUIRE(abs(rf_grid->tf_r[1][0] - rf0->tf_r[1][0]) > 1e-8); // tf_r[ind_ri][ind_re]

  INFO("but not change for any other emitting radius");
  REQUIRE(abs(rf_grid->tf_r[1][1] - rf0->tf_r[1][1]) < 1e-8);

  delete_returningFractions(&rf0);
  delete_returningFractions(&rf_grid);
}

/* disabled test of the interpolation of the returning fractions (uses get_test_emis_rrad)

static emisProfile* get_test_emis_rrad(double rin, double rout, double spin, int* status)
{
    returningFractions* rf = get_rrad_fractions(spin, rin, rout, status);
    double gamma = 2;

    emisProfile* emis = new_emisProfile(rf->rad.data(), rf->nrad, status);
    get_emis_bkn(emis->emis, emis->re, emis->nr, 3.0, 3.0, emis->re[0]);
    for (int ii = 0; ii < rf->nrad; ii++)
    {
        emis->emis[ii] /= emis->emis[0]; //
    }

    emisProfile* emis_return = calc_rrad_emis_corona(rf, nullptr, emis, gamma, status);

    delete_returningFractions(&rf);
    free_emisProfile(emis);

    return emis_return;
}

TEST_CASE(" Interpolation of Returning Radiation Fractions", "[returnrad]") {

  int status = EXIT_SUCCESS;
  double spin = 0.995;
  int ind_rad_test = 2;
  double rout = 1000;

  returningFractions *rf = get_rrad_fractions(spin, kerr_rms(spin) , rout, &status);
  double rad_grid = rf->rad[ind_rad_test];

  double rad_lo = rad_grid*0.99999;
  double rad_hi  = rad_grid*1.00001;

  emisProfile *emis_lo = get_test_emis_rrad(rad_lo, rout, spin, &status);
  emisProfile *emis_0 = get_test_emis_rrad(rad_grid, rout, spin, &status);
  emisProfile *emis_hi = get_test_emis_rrad(rad_hi, rout, spin, &status);

  double sum_lo = calcSum(emis_lo->emis, emis_lo->nr);
  double sum_0 = calcSum(emis_0->emis, emis_0->nr);
  double sum_hi = calcSum(emis_hi->emis, emis_hi->nr);

  REQUIRE(sum_0 > 0);
  REQUIRE(sum_lo / sum_0 > 1.0001);
  REQUIRE(sum_0 / sum_hi > 1.0001);

  // absolute difference for a tiny change around the grid-point less than 1%
  REQUIRE(sum_lo - sum_0 < 0.01);
  REQUIRE(sum_0 - sum_hi < 0.01);

  // absolute difference for a tiny change around the grid-point less than 1%
  REQUIRE(sum_lo / sum_0 < 1.01);
  REQUIRE(sum_0 / sum_hi < 1.01);

  // fractional change in both direction of the grid point should be fairly similar
  REQUIRE( abs( sum_lo/sum_0 - sum_0/sum_hi) < 1e-4 );

} */
