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
#include "IonGradient.h"
#include "common-functions.h"

#define PREC 1e-6


TEST_CASE(" Ion Grad PL Index ", "[iongrad]") {

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  int status = EXIT_SUCCESS;

  LocalModel local_model{ModelName::relxilllpCp};

  local_model.eval_model(spec);

  xillParam* xill_param = local_model.get_xill_params();
  relParam* rel_param = local_model.get_rel_params();

  relline_spec_multizone *rel_profile = relbase(spec.energy, spec.num_flux_bins(), rel_param, nullptr, &status);

  IonGradient ion_gradient{rel_profile->rgrid, rel_profile->n_zones,xill_param->ion_grad_type};
  ion_gradient.calculate(rel_param, xill_param);

  for (int ii=0; ii<rel_profile->n_zones; ii++){
    REQUIRE(ion_gradient.dens[ii] >=15.0);
    REQUIRE(ion_gradient.dens[ii] <=22.0);
    REQUIRE(ion_gradient.lxi[ii] >=0.0);
    REQUIRE(ion_gradient.lxi[ii] <=4.7);
  }

}


TEST_CASE(" Test Alpha Model (writing output) ", "[iongrad]") {

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllpCp};

  const char* env_outfiles = "RELXILL_WRITE_OUTFILES";
  setenv(env_outfiles, "1", 1);
  setenv("RELXILL_NUM_RZONES","50",1);

  local_model.set_par(XPar::logn,19.0);
  local_model.set_par(XPar::h, 6);
  local_model.set_par(XPar::a, 0.9);
  local_model.set_par(XPar::switch_iongrad_type,2);
  local_model.eval_model(spec);

  unsetenv(env_outfiles);
}


TEST_CASE(" Exec single iongrad model, should change with xindex","[iongrad]") {
  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxilllpCp);

  auto spec = default_spec.get_xspec_spectrum();
  lmod.set_par(XPar::switch_iongrad_type, 1);
  lmod.set_par(XPar::iongrad_index,1.0);

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  double sum = sum_flux(spec.flux, spec.num_flux_bins());

  auto spec2 = default_spec.get_xspec_spectrum();
  lmod.set_par(XPar::iongrad_index,2.0);

  REQUIRE_NOTHROW(lmod.eval_model(spec2));
  double sum2 = sum_flux(spec2.flux, spec2.num_flux_bins());

  REQUIRE( sum > 1e-8);
  REQUIRE( sum2 > 1e-8);

  REQUIRE( fabs(sum - sum2) > 1e-8);

}
