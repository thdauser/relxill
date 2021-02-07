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
//#include "LocalModel.h"
//#include "XspecSpectrum.h"
//#include "common-functions.h"
extern "C" {
#include "relreturn.h"
#include "relutility.h"
}

#include <vector>

#define PREC 1e-6


static int is_grid_consistent(const double* rlo, const double* rhi, int nrad){
  for (int ii=1; ii<nrad; ii++){
    if (fabs(rhi[ii-1]-rlo[ii])>PREC){
      printf(" *** error: lower boundary rlo[ii]-rhi[ii-1] = %e - %e = %e  for ii=%i \n",
             rlo[ii], rhi[ii-1],
             fabs(rhi[ii-1]-rlo[ii]), ii);
      return 0;
    }
  }
  return 1;
}

static void test_table_radial_grid(returnFracData *dat) {

  REQUIRE(dat->nrad == RETURNRAD_TABLE_NR);
  REQUIRE(dat->ng == RETURNRAD_TABLE_NG);

  REQUIRE(dat->rlo[0] > 1.0);
  REQUIRE(dat->rhi[dat->nrad - 1] <= 1000.1); // don't care about border effects so add the 0.1

  REQUIRE(is_grid_consistent(dat->rlo, dat->rhi, dat->nrad)==1);
}

static int test_table_fracE_norm(const returnFracData *dat) {

  int status = EXIT_SUCCESS;
  double kSumfrac;
  double kSumfrac_ref = 1.0;

  for (int jj=0; jj<dat->nrad; jj++) {

    kSumfrac=0.0;
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

static int test_table_fracG_norm(const returnFracData *dat) {

  int status = EXIT_SUCCESS;
  double kSumfrac;
  double kSumfrac_ref = 1.0;

  for (int ii=0; ii<dat->nrad; ii++) {
    for (int jj=0; jj<dat->nrad; jj++) {

      kSumfrac = 0.0;
      for (int kk = 0; kk < dat->ng; kk++) {
        kSumfrac += dat->frac_g[ii][jj][kk];
      }

      if (fabs(kSumfrac-kSumfrac_ref)>PREC){
        RELXILL_ERROR("testing the normalization of FRAC_G failed (return radiation table)",&status);
        printf (" [%02i][%02i] expecting a normalization of %e, but found %e\n", ii,jj,kSumfrac_ref, kSumfrac);
        break;
      }
    }
  }
  return status;
}


static void test_table_fractions_normalization(returnFracData* dat){
  REQUIRE(dat->frac_e[0][0]!= 0);
  REQUIRE(test_table_fracE_norm(dat) == EXIT_SUCCESS);
  REQUIRE(test_table_fracG_norm(dat) == EXIT_SUCCESS);
}



// ------- //
TEST_CASE(" Returning Radiation Table ", "[table]") {

  int status = EXIT_SUCCESS;
  returnTable *tab = get_returnRadTable(&status);
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

