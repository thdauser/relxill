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

// Unit tests of the full relxillBB and relxillBB_NS models (evaluated through LocalModel): evaluation, caching,
// the boost parameter, and the environment variables RELXILL_RRAD_BBODY_TPROFILE and RELXILL_RRAD_BBODY_ORDER.
// Run with bin/tests "[rrad-bbody][model]".

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "XspecSpectrum.h"

extern "C" {
#include "relutility.h"
}

#define PREC 1e-6

TEST_CASE (
"relxillBB, relxillBB_NS: model evaluates without error"
,
"[rrad-bbody][model]"
)
 {
  const ModelName mod_name = GENERATE(ModelName::relxillBB, ModelName::relxillBB_NS);
  LocalModel lmod(mod_name);
  INFO("model " << lmod.get_model_string());

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  REQUIRE_NOTHROW(lmod.eval_model(spec));
}

TEST_CASE (
"relxillBB: second evaluation gives the identical result"
,
"[rrad-bbody][model]"
)
 {
  LocalModel lmod(ModelName::relxillBB);

  DefaultSpec default_spec1{};
  XspecSpectrum spec1 = default_spec1.get_xspec_spectrum();
  lmod.eval_model(spec1);

  DefaultSpec default_spec0{};
  XspecSpectrum spec0 = default_spec0.get_xspec_spectrum();
  lmod.eval_model(spec0);

  auto sum0 = calcSum(spec0.flux, spec0.num_flux_bins());
  REQUIRE(std::isnormal(sum0));
  auto sum1 = calcSum(spec1.flux, spec0.num_flux_bins());
  REQUIRE(fabs(sum0 - sum1) < PREC);
}

TEST_CASE (
"relxillBB: boost parameter, direct + reflected = total"
,
"[rrad-bbody][model]"
)
 {
  LocalModel lmod(ModelName::relxillBB);
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  lmod.set_par(XPar::boost, 1);
  lmod.eval_model(spec);
  double sum_total = calcSum(spec.flux, spec.num_flux_bins());
  REQUIRE(sum_total > 1e-8);

  lmod.set_par(XPar::boost, -1);
  lmod.eval_model(spec);
  double sum_refl = calcSum(spec.flux, spec.num_flux_bins());
  REQUIRE(sum_refl > 1e-8);

  lmod.set_par(XPar::boost, 0);
  lmod.eval_model(spec);
  double sum_direct = calcSum(spec.flux, spec.num_flux_bins());
  REQUIRE(sum_direct > 1e-8);

  // check that boost parameter works: refl+direct has to give total flux
  INFO(" total: " << sum_total << " direct: " << sum_direct << " refl: " << sum_refl
    << " (diff: " << fabs(sum_total - sum_direct - sum_refl) / sum_total << ")");
  REQUIRE(fabs(sum_total - sum_direct - sum_refl) / sum_total < 1e-4);
}

TEST_CASE (
"relxillBB: temperature profile (RELXILL_RRAD_BBODY_TPROFILE) changes the spectrum, 1 is the default"
,
"[rrad-bbody][model]"
)
 {
  LocalModel lmod(ModelName::relxillBB);
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  lmod.eval_model(spec);
  double sum_default = calcSum(spec.flux, spec.num_flux_bins());

  setenv("RELXILL_RRAD_BBODY_TPROFILE", "1", 1);
  lmod.eval_model(spec);
  double sum_1 = calcSum(spec.flux, spec.num_flux_bins());

  // check of setting it to 1 is the default
  REQUIRE(fabs(sum_default - sum_1) < 1e-6);

  setenv("RELXILL_RRAD_BBODY_TPROFILE", "2", 1);
  lmod.eval_model(spec);
  double sum_2 = calcSum(spec.flux, spec.num_flux_bins());

  unsetenv("RELXILL_RRAD_BBODY_TPROFILE");

  // check that 1 and 2 are different
  REQUIRE(fabs(sum_2 - sum_1) > 1e-6);
}

TEST_CASE (
"relxillBB: 2nd order returning radiation (RELXILL_RRAD_BBODY_ORDER=2) changes the spectrum"
,
"[rrad-bbody][model]"
)
 {
  LocalModel lmod(ModelName::relxillBB);
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  setenv("RELXILL_RRAD_BBODY_ORDER", "1", 1);

  lmod.eval_model(spec);
  double sum_1st = calcSum(spec.flux, spec.num_flux_bins());

  setenv("RELXILL_RRAD_BBODY_ORDER", "2", 1);
  lmod.eval_model(spec);
  double sum_2nd = calcSum(spec.flux, spec.num_flux_bins());

  REQUIRE((sum_2nd - sum_1st) > 1e-3);

  unsetenv("RELXILL_RRAD_BBODY_ORDER");
}
