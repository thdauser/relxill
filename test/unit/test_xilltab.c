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
#include "relbase.h"
#include "test_relxill.h"

static void print_status_outcome(const int *status) {
  if (*status != EXIT_SUCCESS) {
    printf("  FAILED \n\n");
  } else {
    printf("  SUCCESSFUL \n\n");
  }
}

static void testCreateNewXilltable(int *status) {

  CHECK_STATUS_VOID(*status);

  PRINT_RELXILL_TEST_MSG_DEFAULT();

  int dim[] = {5, 6};
  int ndim = 2;

  int ii;
  for (ii = 0; ii < ndim; ii++) {
    xillTable *tab = new_xillTable(dim[ii], status);

    free_xillTable(tab);
  }

  print_status_outcome(status);

}


static void test_tab_num_param(const xillParam *param, int *status, const xillTable *tab) {
  printf("     - number of parameters: %i ", tab->num_param);
  int dim_expect = 5;
  if (is_6dim_table(param->model_type)) {
    dim_expect = 6;
  }
  printf(" [expected: %i] ", dim_expect);
  if (tab->num_param != dim_expect) {
    *status = EXIT_FAILURE;
  }
}

static void testInitializeXillverTableCorrectly(int *status) {

  PRINT_RELXILL_TEST_MSG_DEFAULT();

  xillTable *tab = NULL;
  xillParam *param = get_std_param_xillver(status);
  init_xillver_table(XILLTABLE_FILENAME, &tab, status);

  if (*status == EXIT_SUCCESS) {
    assert(tab->n_ener > 0);

    int ii;
    for (ii = 0; ii < tab->num_param; ii++) {
      assert(tab->num_param_vals[ii] > 0);
    }

    assert(tab->num_param > 0);
    test_tab_num_param(param, status, tab);

  }

  print_status_outcome(status);

}

static void test_spec_norm(xillSpec *spec, int *status) {

  CHECK_STATUS_VOID(*status);

  int ii;
  int jj;

  double sum = 0.0;

  for (jj = 0; jj < spec->n_incl; jj++) {
    for (ii = 0; ii < spec->n_ener; ii++) {
      assert(spec->flu[jj][ii] >= 0);
      sum += spec->flu[jj][ii];
    }
  }

  printf("Sum = %.4e ", sum);

  if (sum < 1e-8) {
    *status = EXIT_FAILURE;
    printf("\n *** normalization of spectrum seems wrong \n");
    assert(sum >= 1e-8);
  }

}

static void testAutomaticLoadingTable(xillParam *param, int *status) {

  xillTable *tab = NULL;
  char *fname = get_init_xillver_table(&tab, param, status);

  printf("  -> loaded table for model type = %i: ", param->model_type);
  assert(fname != NULL);
  printf("%s\n\n", fname);

  assert(tab != NULL);

}

static void test_get_spec(int *status, char *specName, xillParam *param) {

  CHECK_STATUS_VOID(*status);

  printf("\n - loading spec for model %s (type = %i): ", specName,
         param->model_type);

  xillSpec *spec = get_xillver_spectra(param, status);
  if (*status != EXIT_SUCCESS) {
    printf("\n *** ERROR trying to load the spectrum \n");
  }

  assert(spec->n_ener > 0);
  assert(spec->n_incl > 0);

  test_spec_norm(spec, status);
  free_xill_spec(spec);

}

static void testEvaluateAllXillverSpecNonZero(int *status) {

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  xillParam *param;

  //putenv("DEBUG_RELXILL=1");

  param = get_std_param_xillver(status);
  test_get_spec(status, "xillver", param);

  param = get_std_param_xillver_dens(status);
  test_get_spec(status, "xillverD", param);

  param = get_std_param_xillver_co(status);
  test_get_spec(status, "xillverCO", param);

  param = get_std_param_xillver_ns(status);
  test_get_spec(status, "xillverNS", param);

  param = get_std_param_xillver_nthcomp(status);
  test_get_spec(status, "xillverCp", param);

  param = get_std_param_xillver_dens_nthcomp(status);
  test_get_spec(status, "xillverDCp", param);

  putenv("DEBUG_RELXILL=0");

  if (param != NULL) {
    param->dens = 17.0;
  }
  test_get_spec(status, "xillverDCp (dens=1e17)", param);

  free(param);

  print_status_outcome(status);

}

static void testAllXilltables(int *status) {

  CHECK_STATUS_VOID(*status);
  xillParam *param;

  PRINT_RELXILL_TEST_MSG(" output tables and parameters \n");

  putenv("DEBUG_RELXILL=1");

  param = get_std_param_xillver(status);
  testAutomaticLoadingTable(param, status);

  param = get_std_param_xillver_dens(status);
  testAutomaticLoadingTable(param, status);

  param = get_std_param_xillver_co(status);
  testAutomaticLoadingTable(param, status);

  param = get_std_param_xillver_ns(status);
  testAutomaticLoadingTable(param, status);

  param = get_std_param_xillver_dens_nthcomp(status);
  testAutomaticLoadingTable(param, status);

  free(param);

  putenv("DEBUG_RELXILL=0");

}

static void testLoadingNonExistingTable(int *status) {

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  char *nonExistingFilename = "no_table_has_this_name_1234.fits";

  int statusFailing = EXIT_SUCCESS;
  xillTable *tab = NULL;
  init_xillver_table(nonExistingFilename, &tab, &statusFailing);

  if (statusFailing == EXIT_SUCCESS) {
    *status = EXIT_FAILURE;
    printf("\n *** error : loading a non-existing table did not result in an error\n");
  }

  print_status_outcome(status);

}

static void testCheckForExistingTable(int *status) {

  CHECK_STATUS_VOID(*status);

  printf(" TEST: checking for table existence ");

  char *nonExistingTable = "no_table_has_this_name_1234.fits";
  if (checkIfTableExists(nonExistingTable, status) != 0) {
    *status = EXIT_FAILURE;
    printf(" TEST: *** error : test for a not existing table has failed\n");
  }

  char *existingTable = XILLTABLE_FILENAME;
  if (checkIfTableExists(existingTable, status) != 1) {
    *status = EXIT_FAILURE;
    printf(" TEST: *** error : trying to open table %s failed\n", existingTable);
  }

  print_status_outcome(status);

}

static void testLoadingAlternativeTable(int *status) {

  CHECK_STATUS_VOID(*status);

  printf(" TEST: loading alternative table ");

  char *existingTable = XILLTABLE_FILENAME;
  char *nonExistingTable = "no_table_has_this_name_1234.fits";

  char *tableName = getXilltableNameUsingAlternativeIfNotExisting(
      nonExistingTable, existingTable, status);

  if (checkIfTableExists(tableName, status) != 1) {
    *status = EXIT_FAILURE;
    printf(" TEST: *** error : trying loading alternative table %s instead of %s failed\n",
           existingTable, nonExistingTable);
  }

  print_status_outcome(status);

}

static void testNormfacBand(xillSpec **spec, double elo, double ehi, double prec) {

  double sum0 = calcSumInEnergyBand(spec[0]->flu[0], spec[0]->n_ener, spec[0]->ener, elo, ehi);
  double sum1 = calcSumInEnergyBand(spec[1]->flu[0], spec[1]->n_ener, spec[0]->ener, elo, ehi);

  double fraction = sum0 / sum1 - 1;

  assert(fabs(fraction) < prec);

}

static void testRenormXilltableDensLogxi(int *status) {

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  xillParam *param;

  xillSpec *spec[2];
  param = get_std_param_xillver_dens_nthcomp(status);
  spec[0] = get_xillver_spectra(param, status);

  test_get_spec(status, "xillverDCp", param);
  param->dens = 15.5;
  spec[1] = get_xillver_spectra(param, status);

  testNormfacBand(spec, 30.0, 80.0, 0.01);

  print_status_outcome(status);

}

void test_xilltables(int *status) {
  if (*status != EXIT_SUCCESS) {
    printf(" *** SKIP testing Xilltables as an error occured previously! \n\n");
    return;
  }

  printf("\n === Testing XILLVER with RELXILL Version %s === \n\n", PROJECT_VER);

  testCreateNewXilltable(status);

  testInitializeXillverTableCorrectly(status);

  testLoadingNonExistingTable(status);

  testLoadingAlternativeTable(status);

  testCheckForExistingTable(status);

  testAllXilltables(status);

  testEvaluateAllXillverSpecNonZero(status);

  testRenormXilltableDensLogxi(status);

  if (*status != EXIT_SUCCESS) {
    printf(" *** TESTING XILLVER TABLES NOT SUCCESSFUL \n");
    printf("\n *** Cleaning up and freeing cached structures\n");
  } else {
    printf(" *** TESTING XILLVER TABLES SUCCESSFUL \n");
  }

}
