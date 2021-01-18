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


#include "test_rellp.h"

#define PREC 1e-6

static double sumArray(double *arr, int n) {
  double sum = 0.0;
  for (int ii = 0; ii < n; ii++) {
    sum += arr[ii];
  }
  return sum;
}

static relParam *get_std_relParam_relxilllpDCp(int *status) {
  relParam *rel_param = NULL;
  xillParam *xill_param = NULL;
  get_std_param_relxilllpDCp(&rel_param, &xill_param, status);
  return rel_param;
}

static void test_pointSourceDecision(int *status) {

  printf(" *** testing: pointSourcedDecision\n");

  relParam *rel_param = get_std_relParam_relxilllpDCp(status);

  rel_param->height = 3.0;
  rel_param->htop = 3.0;
  assert(modelLampPostPointsource(rel_param) == 1);

  rel_param->htop = 10.0;
  assert(modelLampPostPointsource(rel_param) == 0);

  rel_param->htop = 1.0;
  assert(modelLampPostPointsource(rel_param) == 1);

  rel_param->htop = 0.0;
  assert(modelLampPostPointsource(rel_param) == 1);

  CHECK_RELXILL_DEFAULT_ERROR(status);

}

static void test_extendedGeometryHeight(int *status) {

  printf(" *** testing: extendedGeometryHeight\n");

  relParam *rel_param = get_std_relParam_relxilllpDCp(status);

  rel_param->height = 3.0;
  rel_param->htop = 10.0;

  extPrimSource *source = getExtendedJetGeom(rel_param, status);
  assert(*status == EXIT_SUCCESS);

  // test that hbase and htop are set correctly
  assert(fabs(source->heightArr[0] - rel_param->height) < PREC);
  assert(fabs(source->heightArr[source->nh] - rel_param->htop) < PREC);  // we have nh+1 elements

  // test that the height values are ascending
  for (int ii = 0; ii < source->nh; ii++) {
    assert(source->heightArr[ii] < source->heightArr[ii + 1]);
  }

  // test that values are sensible
  assert(source->heightArr[0] >= 1.0);
  assert(source->heightArr[source->nh] <= 1000.0);

}

static void test_calcEmisProfileLp(int *status) {

  printf(" *** testing: calcEmisProfileLp \n");

  relParam *rel_param = get_std_relParam_relxilllpDCp(status);

  rel_param->height = 3.0;
  rel_param->htop = 3.0;
  RelSysPar *sysPar = get_system_parameters(rel_param, status);
  assert(sysPar->emis->normFactorPrimSpec > 0);

  double sumEmisPoint = sumArray(sysPar->emis->emis, sysPar->emis->nr);
  assert(sumEmisPoint > 0);

  rel_param->height = 3.0;
  rel_param->htop = 10.0;
  assert(modelLampPostPointsource(rel_param) == 0);
  RelSysPar *sysPar2 = get_system_parameters(rel_param, status);

  assert(sysPar2->emis->normFactorPrimSpec > 0);
  double sumEmisExt = sumArray(sysPar2->emis->emis, sysPar->emis->nr);
  assert(sumEmisExt > 0);

  assert(sumEmisExt != sumEmisPoint);

  assert(sysPar2->emis->normFactorPrimSpec != sysPar->emis->normFactorPrimSpec);
  assert(sysPar2->emis->returnFracs->refl_frac != sysPar->emis->returnFracs->refl_frac);

  CHECK_RELXILL_DEFAULT_ERROR(status);

}

void test_rellp(int *status) {

  CHECK_STATUS_VOID(*status);

  char *testName = "EMISSIVTIY Profile Calculation";
  printf("\n === Testing %s with RELXILL Version %s === \n\n", testName, PROJECT_VER);

  test_extendedGeometryHeight(status);

  test_calcEmisProfileLp(status);

  test_pointSourceDecision(status);

  if (*status != EXIT_SUCCESS) {
    printf(" *** TESTING %s NOT SUCCESSFUL \n", testName);
  } else {
    printf(" *** TESTING %s SUCCESSFUL \n", testName);
  }

}
