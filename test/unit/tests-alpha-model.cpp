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
  lmod.set_par(XPar::rout, 400);
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

TEST_CASE(" Calculate L_source", "[alpha]") {

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

double calculate_lxi(const LocalModel &lmod) {

  auto rel_param = lmod.get_rel_params();

  int status = EXIT_SUCCESS;
  RelSysPar *sys_par = get_system_parameters(rel_param, &status);

  auto radial_grid = RadialGrid(rel_param->rin, rel_param->rout, rel_param->num_zones, rel_param->height);
  IonGradient const ion_gradient{radial_grid, rel_param->ion_grad_type, 0};

  const double logn = lmod.get_model_params().get_par(XPar::logn);

  auto primary_source_params = PrimarySourceParameters(lmod.get_model_params());
  const double lxi_max = IonGradient::calculate_lxi_max_from_distance(*(sys_par->emis),
                                                                      primary_source_params,
                                                                      logn);

  return lxi_max;

}

TEST_CASE(" Calculate lxi ", "[alpha]") {

  int status = EXIT_SUCCESS;

  LocalModel lmod_alpha(ModelName::relxilllpAlpha);
  set_default_par(lmod_alpha);
  lmod_alpha.set_par(XPar::distance, 10);  // distance of 2 kpc
  lmod_alpha.set_par(XPar::mass, 20);      // mass of 20 Msolar
  lmod_alpha.set_par(XPar::h, 6.0);  // to be in the Newtonian regime
  lmod_alpha.set_par(XPar::norm_flux_cgs, 1e-10); // Flux in the 0.01-1000keV band in erg/cm^2/s
  lmod_alpha.set_par(XPar::logn, 18);

  const double ref_value_l_source = 2.5972e36;
  auto primary_source_params = PrimarySourceParameters(lmod_alpha.get_model_params());
  auto sys_par = get_system_parameters(lmod_alpha.get_rel_params(), &status);
  REQUIRE(1e-3 > fabs(
      ref_value_l_source /
          primary_source_params.luminosity_source_cgs((*sys_par->emis->photon_fate_fractions))
          - 1)
  );

  const double lxi_max = calculate_lxi(lmod_alpha);

  // ref value (not taking into account the delt_inc correction, so this is only estimated)
  const double ref_value_emis_lxi_max = 6.83e+21;

  const double logn_rin = lmod_alpha.get_model_params().get_par(XPar::logn);
  const double rin = lmod_alpha.get_model_params().get_par(XPar::rin);
  const double r_lximax = pow(11. / 9, 2) * rin;

  const double density_r_lximax = density_ss73_zone_a(r_lximax, rin) * pow(10, logn_rin);
  // const double ref_value_density = 5.5230109739368972e+19;
  // REQUIRE(fabs(ref_value_density/density_r_lximax - 1) < 1e-4);

  const double ref_approx_lxi_max = log10(4.0 * M_PI * ref_value_emis_lxi_max / density_r_lximax);

  // although we do not use the delta_inc correction, the value should be within 10% precise
  REQUIRE(fabs(lxi_max / ref_approx_lxi_max - 1) < 0.1);

}

TEST_CASE(" Flux Normalization of the Continuum", "[alpha-test]") {
  auto default_spec = DefaultSpec(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, 3000);

  LocalModel lmod_alpha(ModelName::relxilllpAlpha);
  set_default_par(lmod_alpha);
  lmod_alpha.set_par(XPar::distance, 1e5);  // distance of 100 Mpc
  lmod_alpha.set_par(XPar::mass, 1e6);      // mass of 1e6 Msolar


  lmod_alpha.set_par(XPar::refl_frac, 0.0);

  const double norm_factor_ergs = 1e-11; // Flux in the 0.01-1000keV band in erg/cm^2/s
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

  // calculate logxi from the alpha distance model and then use it as input to the relxilllpCp model
  // this should give the same result
  const double lxi_max_alpha = calculate_lxi(lmod_alpha);
  LocalModel lmod_std(ModelName::relxilllpCp);
  set_default_par(lmod_std);
  lmod_std.set_par(XPar::logxi, lxi_max_alpha);
  lmod_std.set_par(XPar::switch_iongrad_type, 2.0);
  lmod_std.set_par(XPar::refl_frac, 0.0);

  lmod_std.eval_model(spec);
  const double std_primary_energy_flux = spec.get_energy_flux();
  lmod_std.set_par(XPar::refl_frac, -1.0);
  lmod_std.eval_model(spec);
  const double std_refl_energy_flux = spec.get_energy_flux();

  const double flux_ratio_alpha = refl_energy_flux / primary_energy_flux;
  const double flux_ratio_std = std_refl_energy_flux / std_primary_energy_flux;

  // now test that the flux ratio to the reflection spectrum is the same as in
  REQUIRE(fabs(flux_ratio_alpha / flux_ratio_std - 1) < 0.01);

}
