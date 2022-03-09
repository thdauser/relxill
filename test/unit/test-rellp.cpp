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

#include "LocalModel.h"
#include "xspec_wrapper_lmodels.h"
#include "XspecSpectrum.h"
#include "common-functions.h"
#include "Rellp.h"

#define PREC 1e-6

TEST_CASE(" test printing of reflection strength", "[rellp]" ){

  DefaultSpec default_spec{};
  LocalModel lmod(ModelName::relxilllp);

  const char *env = "RELXILL_PRINT_DETAILS";

  lmod.set_par(XPar::refl_frac, 1);
  lmod.set_par(XPar::switch_switch_reflfrac_boost, 1);
  auto spec = default_spec.get_xspec_spectrum();

  lmod.eval_model(spec);

  setenv(env, "1", 1);

  lmod.eval_model(spec);

  unsetenv(env);

}

/*

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
  assert(sysPar2->emis->photon_fate_fractions->refl_frac != sysPar->emis->photon_fate_fractions->refl_frac);

  CHECK_RELXILL_DEFAULT_ERROR(status);

}
*/

