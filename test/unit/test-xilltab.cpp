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

extern "C" {
#include "relbase.h"
#include "test_relxill.h"
}


TEST_CASE(" creating of new xillver table ", "[xilltab]") {
  int status = EXIT_SUCCESS;

  int dim[] = {5, 6};
  int ndim = 2;

  int ii;
  for (ii = 0; ii < ndim; ii++) {
    xillTable *tab = new_xillTable(dim[ii], &status);
    REQUIRE(status == EXIT_SUCCESS);
    free_xillTable(tab);
  }

}


static void test_tab_num_param(const xillParam *param, int *status, const xillTable *tab) {

  int dim_expect = 5;
  if (is_6dim_table(param->model_type)) {
    dim_expect = 6;
  }
  REQUIRE(tab->num_param == dim_expect);
}


TEST_CASE(" initialize Xillver table ", "[xilltab]") {

  int status = EXIT_SUCCESS;

  xillTable *tab = nullptr;
  xillParam *param = get_std_param_xillver(&status);

  const std::string fname_xilltable = XILLTABLE_FILENAME;
  init_xillver_table(fname_xilltable.c_str(), &tab, &status);

  REQUIRE( status == EXIT_SUCCESS);
  REQUIRE(tab->n_ener > 0);

  int ii;
  for (ii = 0; ii < tab->num_param; ii++) {
    REQUIRE(tab->num_param_vals[ii] > 0);
  }

  REQUIRE(tab->num_param > 0);
  test_tab_num_param(param, &status, tab);

}

static void test_spec_norm(xillSpec *spec) {

  int ii;
  int jj;

  double norm_spectrum = 0.0;

  for (jj = 0; jj < spec->n_incl; jj++) {
    for (ii = 0; ii < spec->n_ener; ii++) {
      REQUIRE(spec->flu[jj][ii] >= 0);
      norm_spectrum += spec->flu[jj][ii];
    }
  }
  REQUIRE(norm_spectrum > 1e-8);

}

static void testAutomaticLoadingTable(xillParam *param) {

  xillTable *tab = nullptr;

  int status = EXIT_SUCCESS;
  const char *fname = get_init_xillver_table(&tab, param, &status);
  REQUIRE(fname != nullptr);
  REQUIRE( status == EXIT_SUCCESS );
  INFO("  -> loading table for model type = %i: " << fname);
  REQUIRE(tab != nullptr);

}

static void test_get_spec(const std::string& specName, xillParam *param) {


  INFO(" loading spec for model %s (type = %i): " << specName);

  int status = EXIT_SUCCESS;

  xillSpec *spec = get_xillver_spectra(param, &status);
  REQUIRE( status == EXIT_SUCCESS);

  REQUIRE(spec->n_ener > 0);
  REQUIRE(spec->n_incl > 0);

  test_spec_norm(spec);
  free_xill_spec(spec);

}

TEST_CASE(" xillver spectra evaluation ", "[xilltab]") {

  xillParam *param;

  int status = EXIT_SUCCESS;

  param = get_std_param_xillver(&status);
  REQUIRE(status == EXIT_SUCCESS);
  test_get_spec("xillver", param);

  param = get_std_param_xillver_co(&status);
  REQUIRE(status == EXIT_SUCCESS);
  test_get_spec("xillverCO", param);

  param = get_std_param_xillver_ns(&status);
  REQUIRE(status == EXIT_SUCCESS);
  test_get_spec("xillverNS", param);

  param = get_std_param_xillver_nthcomp(&status);
  REQUIRE(status == EXIT_SUCCESS);
  test_get_spec("xillverCp", param);

  param = get_std_param_xillver_dens_nthcomp(&status);
  REQUIRE(status == EXIT_SUCCESS);
  test_get_spec("xillverDCp", param);

  if (param != nullptr) {
    param->dens = 17.0;
  }
  test_get_spec("xillverDCp (dens=1e17)", param);

  free(param);

}

TEST_CASE(" automated loading of the xillver tables ", "[xilltab]") {

  int status = EXIT_SUCCESS;
  xillParam *param;


  param = get_std_param_xillver(&status);
  REQUIRE(status == EXIT_SUCCESS);
  testAutomaticLoadingTable(param);
  free(param);

  param = get_std_param_xillver_dens(&status);
  REQUIRE(status == EXIT_SUCCESS);
  testAutomaticLoadingTable(param);
  free(param);

  param = get_std_param_xillver_co(&status);
  REQUIRE(status == EXIT_SUCCESS);
  testAutomaticLoadingTable(param);
  free(param);

  param = get_std_param_xillver_ns(&status);
  REQUIRE(status == EXIT_SUCCESS);
  testAutomaticLoadingTable(param);
  free(param);

  param = get_std_param_xillver_dens_nthcomp(&status);
  REQUIRE(status == EXIT_SUCCESS);
  testAutomaticLoadingTable(param);
  free(param);


}

TEST_CASE(" loading table which does not exist ", "[xilltab]") {

  int status = EXIT_SUCCESS;
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  std::string nonExistingFilename = "no_table_has_this_name_1234.fits";

  int statusFailing = EXIT_SUCCESS;
  xillTable *tab = nullptr;
  init_xillver_table(nonExistingFilename.c_str(), &tab, &statusFailing);

  REQUIRE( statusFailing == EXIT_FAILURE);

}


TEST_CASE(" check for existing table ", "[xilltab]") {

  int status = EXIT_SUCCESS;

  printf(" TEST: checking for table existence ");

  const std::string nonExistingTable = "no_table_has_this_name_1234.fits";
  REQUIRE(checkIfTableExists(nonExistingTable.c_str(), &status) == 0);

  const std::string existingTable = XILLTABLE_FILENAME;
  INFO(" TEST: *** error : trying to open table %s failed "<< existingTable);
  REQUIRE(checkIfTableExists(existingTable.c_str(), &status) == 1);

}

TEST_CASE(" loading alternative xillver table", "[xilltab]") {

  int status = EXIT_SUCCESS;
  const std::string existingTable = XILLTABLE_FILENAME;
  const std::string nonExistingTable = "no_table_has_this_name_1234.fits";

  const char *tableName = getXilltableNameUsingAlternativeIfNotExisting(
      nonExistingTable.c_str(), existingTable.c_str(), &status);

  REQUIRE(checkIfTableExists(tableName, &status) == 1);

}

static void testNormfacBand(xillSpec **spec, double elo, double ehi, double prec) {

  double sum0 = calcSumInEnergyBand(spec[0]->flu[0], spec[0]->n_ener, spec[0]->ener, elo, ehi);
  double sum1 = calcSumInEnergyBand(spec[1]->flu[0], spec[1]->n_ener, spec[0]->ener, elo, ehi);

  double fraction = sum0 / sum1 - 1;

  REQUIRE(fabs(fraction) < prec);

}

TEST_CASE(" renorm xilltable with density and logxi ") {

  int status = EXIT_SUCCESS;
  xillParam *param;

  xillSpec *spec[2];
  param = get_std_param_xillver_dens_nthcomp(&status);
  spec[0] = get_xillver_spectra(param, &status);

  test_get_spec("xillverDCp", param);
  param->dens = 15.5;
  spec[1] = get_xillver_spectra(param, &status);

  testNormfacBand(spec, 30.0, 80.0, 0.01);

}
