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


// Unit tests of the full lamp post models (relxilllp, relxilllpCp) with the returning radiation of the corona, evaluated
// through LocalModel. Run with bin/tests "[rrad-corona][model]".

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "XspecSpectrum.h"
#include "common-functions.h"
#include "writeOutfiles.h"

#include <string>

extern "C" {
#include "relutility.h"
}

TEST_CASE(
  "relxilllp: returning radiation increases the flux",
  "[rrad-corona][model]"
) {
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllp};

  local_model.set_par(XPar::switch_switch_returnrad, 0);
  local_model.eval_model(spec);
  double model_no_rrad_flux = sum_flux(spec.flux, spec.num_flux_bins());

  local_model.set_par(XPar::switch_switch_returnrad, 1);
  local_model.eval_model(spec);
  double model_rrad_flux = sum_flux(spec.flux, spec.num_flux_bins());

  string fname = "test-spec-relxilllpret.dat";
  save_xillver_spectrum(spec.energy, spec.flux, spec.num_flux_bins(), fname.data());

  REQUIRE(model_rrad_flux > model_no_rrad_flux);
}

TEST_CASE(
  "relxilllp: RELXILL_RETURNRAD_SWITCH does not override the returnrad parameter",
  "[rrad-corona][model]"
) {
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllp};

  local_model.set_par(XPar::switch_switch_returnrad, 0);
  local_model.eval_model(spec);
  double model_no_rrad_flux = sum_flux(spec.flux, spec.num_flux_bins());

  const char *env = "RELXILL_RETURNRAD_SWITCH";
  setenv(env, "1", 1);
  local_model.eval_model(spec);
  double model_rrad_flux = sum_flux(spec.flux, spec.num_flux_bins());

  unsetenv(env);

  // env variable can not change if parameter value is given
  REQUIRE(fabs(model_rrad_flux - model_no_rrad_flux) <1e-8);
}

TEST_CASE(
  "relxilllp: changing Rin slightly near the ISCO changes the spectrum with returning radiation",
  "[rrad-corona][model]"
) {
  double spin = 0.995;
  double rin_scaled = kerr_rms(0.99025);
  double rin_scaled2 = kerr_rms(0.990);

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  const char *env_outfiles = "RELXILL_WRITE_OUTFILES";
  setenv(env_outfiles, "1", 1);

  LocalModel local_model{ModelName::relxilllp};
  local_model.set_par(XPar::a, spin);
  local_model.set_par(XPar::switch_switch_returnrad, 1);
  local_model.set_par(XPar::rout, 1000.0);

  local_model.set_par(XPar::rin, rin_scaled2);
  local_model.eval_model(spec);

  double sum1 = sum_flux(spec.flux, spec.num_flux_bins());

  local_model.set_par(XPar::rin, rin_scaled);
  local_model.eval_model(spec);
  double sum2 = sum_flux(spec.flux, spec.num_flux_bins());

  REQUIRE(sum1 > 1e-6);
  REQUIRE(abs(sum1 - sum2) > 1e-4);

  unsetenv(env_outfiles);
}

TEST_CASE(
  "relxilllp: spectrum with returning radiation is identical after changing a parameter and back (caching)",
  "[rrad-corona][model]"
) {
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllp};

  double a_default = 0.998;
  local_model.set_par(XPar::a, a_default);
  local_model.set_par(XPar::h, 3.);
  local_model.set_par(XPar::switch_switch_reflfrac_boost, 1);
  local_model.set_par(XPar::refl_frac, 1);

  local_model.eval_model(spec);
  double model1 = sum_flux(spec.flux, spec.num_flux_bins());

  // reset relxill total cache by this
  local_model.set_par(XPar::a, a_default * 0.99);
  local_model.eval_model(spec);

  // evaluate wih the same values again, should give the same result
  local_model.set_par(XPar::a, a_default);
  local_model.eval_model(spec);
  double model2 = sum_flux(spec.flux, spec.num_flux_bins());
  REQUIRE(fabs(model1 - model2)/model1 < 1e-6);
}

TEST_CASE(
  "relxilllpCp: returning radiation changes the spectrum also for spin 0",
  "[rrad-corona][model]"
) {
  LocalModel local_model{ModelName::relxilllpCp};

  local_model.set_par(XPar::a, 0.0);
  local_model.set_par(XPar::switch_switch_reflfrac_boost, 1);
  local_model.set_par(XPar::refl_frac, -1);
  local_model.set_par(XPar::afe, 1.0);
  local_model.set_par(XPar::h, 3.);
  local_model.set_par(XPar::rin, 10);
  local_model.set_par(XPar::switch_iongrad_type, 0);

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  // from v2.7 use renormalized relxill by default, here we need old normalization
  // because renormalized version does not keep the true flux
  setenv("RELXILL_RENORMALIZE", "0", 0);

  local_model.set_par(XPar::switch_switch_returnrad, 1);
  local_model.eval_model(spec);
  double model1 = sum_flux(spec.flux, spec.num_flux_bins());

  local_model.set_par(XPar::switch_switch_returnrad, 0);
  local_model.eval_model(spec);
  double model2 = sum_flux(spec.flux, spec.num_flux_bins());

  REQUIRE(fabs(model1 - model2)/model1 > 1e-3);

  unsetenv("RELXILL_RENORMALIZE");
}

TEST_CASE(
  "relxilllpCp: evaluates with ionization gradient (type 2) and returning radiation",
  "[rrad-corona][model]"
) {
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllpCp};

  local_model.set_par(XPar::a, 0.8);
  local_model.set_par(XPar::switch_switch_reflfrac_boost, 1);
  local_model.set_par(XPar::refl_frac, 1);
  local_model.set_par(XPar::afe, 10.0);
  local_model.set_par(XPar::h, 6.);
  local_model.set_par(XPar::switch_iongrad_type, 2);
  local_model.set_par(XPar::switch_switch_returnrad, 1);

  //  relParam* rel_param = local_model.get_rel_params();
  //  xillParam * xill_param = local_model.get_xill_params();

  local_model.eval_model(spec);
  double model1 = sum_flux(spec.flux, spec.num_flux_bins());
}
