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

#include <vector>
#include <filesystem>

TEST_CASE(" Execute relxillpAlpha", "[alpha]") {

  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxilllpAlpha);

  auto spec = default_spec.get_xspec_spectrum();

  REQUIRE_NOTHROW(lmod.eval_model(spec));
}

void set_default_par(LocalModel &lmod) {
  lmod.set_par(XPar::a, 0.998);
  lmod.set_par(XPar::h, 60.0);
  lmod.set_par(XPar::afe, 1.0);
  lmod.set_par(XPar::incl, 30.0);
  lmod.set_par(XPar::logxi, 3.1);
  lmod.set_par(XPar::logn, 15);
  lmod.set_par(XPar::gamma, 2.0);
  lmod.set_par(XPar::kte, 60);
  lmod.set_par(XPar::beta, 0.0);
  lmod.set_par(XPar::z, 0.0);
  lmod.set_par(XPar::switch_switch_reflfrac_boost, 1.0);
}

TEST_CASE(" Calculate L_source", "[alpha-lsource]") {

  LocalModel lmod_alpha(ModelName::relxilllpAlpha);
  set_default_par(lmod_alpha);
  lmod_alpha.set_par(XPar::distance, 1e3);  // distance of 1e3 kpc
  lmod_alpha.set_par(XPar::mass, 1e6);      // mass of 1e6 Msolar
  lmod_alpha.set_par(XPar::h, 300.0);  // to be in the Newtonian regime

  const double norm_factor_ergs = 1e-12; // Flux in the 0.01-1000keV band in erg/cm^2/s
  lmod_alpha.set_par(XPar::norm_flux_cgs, norm_factor_ergs);

  auto primary_source_params = PrimarySourceParameters(lmod_alpha.get_model_params());

  int status = EXIT_SUCCESS;
  RelSysPar *sys_par = get_system_parameters(lmod_alpha.get_rel_params(), &status);
  const double l_source = primary_source_params.luminosity_source_cgs((*sys_par->emis->photon_fate_fractions));

  const double REF_value_l_source = 1.196495197217232e+38;
  REQUIRE(fabs(REF_value_l_source / l_source - 1) < 0.01);

}

TEST_CASE(" Flux Normalization of the Continuum", "[alpha]") {
  auto default_spec = DefaultSpec(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, 3000);

  LocalModel lmod_alpha(ModelName::relxilllpAlpha);
  lmod_alpha.set_par(XPar::distance, 1e3);  // distance of 1e3 kpc
  lmod_alpha.set_par(XPar::mass, 1e6);      // mass of 1e6 Msolar


  lmod_alpha.set_par(XPar::refl_frac, 0.0);
  set_default_par(lmod_alpha);

  const double norm_factor_ergs = 1e-12; // Flux in the 0.01-1000keV band in erg/cm^2/s
  lmod_alpha.set_par(XPar::norm_flux_cgs, norm_factor_ergs);

  // make sure the primary continuum normalization is given by XPar::norm_factor
  auto spec = default_spec.get_xspec_spectrum();
  lmod_alpha.eval_model(spec);

  const double primary_energy_flux = spec.get_energy_flux();
  REQUIRE(fabs(primary_energy_flux - norm_factor_ergs) < 0.005);
  // printf(" ENERGY FLUX: %e \n", spec.get_energy_flux());

  lmod_alpha.set_par(XPar::refl_frac, -1.0);
  lmod_alpha.eval_model(spec);
  const double refl_energy_flux = spec.get_energy_flux();

  // TODO: once logxi is not a parameter any more, we need to output the calculated
  // logxi and then use it as input to the relxilllpCp model
  LocalModel lmod_std(ModelName::relxilllpCp);
  lmod_std.set_par(XPar::switch_iongrad_type, 2.0);
  lmod_std.set_par(XPar::refl_frac, 0.0);
  set_default_par(lmod_std);

  lmod_std.eval_model(spec);
  const double std_primary_energy_flux = spec.get_energy_flux();
  lmod_std.set_par(XPar::refl_frac, -1.0);
  lmod_std.eval_model(spec);
  const double std_refl_energy_flux = spec.get_energy_flux();

  const double flux_ratio_alpha = refl_energy_flux / primary_energy_flux;
  const double flux_ratio_std = std_refl_energy_flux / std_primary_energy_flux;

  // now test that the flux ratio to the reflection spectrum is the same as in
  REQUIRE(fabs(flux_ratio_alpha - flux_ratio_std) < 0.01);

}
