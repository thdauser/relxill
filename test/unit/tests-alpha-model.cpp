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
  auto rel_params = lmod_alpha.get_rel_params();
  RelSysPar *sys_par = get_system_parameters(rel_params, &status);
  const double l_source = primary_source_params.luminosity_source_cgs((*sys_par->emis->photon_fate_fractions));

  const double REF_value_l_source = 1.196495197217232e+38;
  REQUIRE(fabs(REF_value_l_source / l_source - 1) < 0.01);

  // delete rel_params;

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

  delete rel_param;

  return lxi_max;

}

TEST_CASE(" Calculate lxi ", "[alpha]") {

  int status = EXIT_SUCCESS;

  LocalModel lmod_alpha(ModelName::relxilllpAlpha);
  set_default_par(lmod_alpha);
  lmod_alpha.set_par(XPar::distance, 10);  // distance of 2 kpc
  lmod_alpha.set_par(XPar::mass, 20);      // mass of 20 Msolar
  lmod_alpha.set_par(XPar::h, 6.0);  // to be in the Newtonian regime
  lmod_alpha.set_par(XPar::norm_flux_cgs, 1e-10); // Flux in the 0.1-1000keV band in erg/cm^2/s
  lmod_alpha.set_par(XPar::logn, 18);
  lmod_alpha.set_par(XPar::switch_switch_reflfrac_boost, 1);
  lmod_alpha.set_par(XPar::refl_frac, -1);

  const double ref_value_l_source = 2.5972e36;
  auto primary_source_params = PrimarySourceParameters(lmod_alpha.get_model_params());
  auto rel_params = lmod_alpha.get_rel_params();
  auto sys_par = get_system_parameters(rel_params, &status);
  REQUIRE(1e-3 >
      fabs(
          ref_value_l_source / primary_source_params.luminosity_source_cgs((*sys_par->emis->photon_fate_fractions)) - 1
      ));

  // delete rel_params;

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

LocalModel get_std_lmod_alpha(double norm_factor_ergs) {
  LocalModel lmod_alpha(ModelName::relxilllpAlpha);
  set_default_par(lmod_alpha);
  lmod_alpha.set_par(XPar::distance, 1e5);  // distance of 100 Mpc
  lmod_alpha.set_par(XPar::mass, 1e6);      // mass of 1e6 Msolar
  lmod_alpha.set_par(XPar::norm_flux_cgs, norm_factor_ergs);

  return lmod_alpha;
}

TEST_CASE(" Flux Normalization of the Continuum", "[alpha]") {
  auto default_spec = DefaultSpec(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, 3000);
  auto spec = default_spec.get_xspec_spectrum();

  const double norm_factor_ergs = 1e-11; // Flux in the 0.01-1000keV band in erg/cm^2/s
  auto lmod_alpha = get_std_lmod_alpha(norm_factor_ergs);

  SECTION("Test that flux normalization works  for the full model, regardless of the reflection fraction") {
    std::array<double, 3> const refl_frac{{0.0, 1.0, -1.34}};

    for (const double rf: refl_frac) {
      INFO(" for reflection fraction " << std::fixed << rf);
      lmod_alpha.set_par(XPar::refl_frac, rf);
      lmod_alpha.eval_model(spec);

      const double energy_flux = spec.get_energy_flux();
      REQUIRE(fabs(energy_flux - norm_factor_ergs) / norm_factor_ergs < 0.005);
    }
  }
}

TEST_CASE("refl/direct flux ratio identical to relxilllpCp Model", "[alpha]") {
  auto default_spec = DefaultSpec(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, 3000);
  auto spec = default_spec.get_xspec_spectrum();

  // (1) calculate alpha flux ratio
  const double norm_factor_ergs = 1e-11; // Flux in the 0.01-1000keV band in erg/cm^2/s
  auto lmod_alpha = get_std_lmod_alpha(norm_factor_ergs);
  lmod_alpha.set_par(XPar::refl_frac, 0.0);
  lmod_alpha.eval_model(spec);
  const double primary_energy_flux = spec.get_energy_flux();

  lmod_alpha.set_par(XPar::refl_frac, -1.0);
  lmod_alpha.eval_model(spec);
  const double refl_energy_flux = spec.get_energy_flux();

  // (2) calculate relxillpCp flux ratio, with logxi from the
  //  alpha distance model and then use it as input to the relxilllpCp model
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


TEST_CASE(" Test the refl_frac parameter works with the alpha model as well", "[alpha]") {
  auto default_spec = DefaultSpec(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, 3000);
  auto spec = default_spec.get_xspec_spectrum();

  LocalModel lmod_alpha(ModelName::relxilllpAlpha);
  set_default_par(lmod_alpha);
  lmod_alpha.set_par(XPar::logn, 18);
  lmod_alpha.set_par(XPar::distance, 1e5);  // distance of 100 Mpc
  lmod_alpha.set_par(XPar::mass, 1e6);      // mass of 1e6 Msolar

  const double norm_factor_ergs = 1e-11; // Flux in the 0.01-1000keV band in erg/cm^2/s
  lmod_alpha.set_par(XPar::norm_flux_cgs, norm_factor_ergs);

  lmod_alpha.set_par(XPar::refl_frac, 1.0);
  lmod_alpha.eval_model(spec);
  const double refl_energy_flux_1 = spec.get_energy_flux();
  const double lxi_max_1 = calculate_lxi(lmod_alpha);

  lmod_alpha.set_par(XPar::refl_frac, -1.0);
  lmod_alpha.eval_model(spec);
  const double refl_energy_flux_1_refl = spec.get_energy_flux();

  lmod_alpha.set_par(XPar::refl_frac, 2.0);
  lmod_alpha.eval_model(spec);
  const double refl_energy_flux_2 = spec.get_energy_flux();
  const double lxi_max_2 = calculate_lxi(lmod_alpha);
  lmod_alpha.set_par(XPar::refl_frac, -2.0);
  lmod_alpha.eval_model(spec);
  const double refl_energy_flux_2_refl = spec.get_energy_flux();

  // the energy flux should be different
  REQUIRE(fabs(refl_energy_flux_1 / refl_energy_flux_2 - 1) > 1e-4);

  // for larger boost, we have a higher ionization
  REQUIRE(lxi_max_2 > lxi_max_1);

  // as the  boost is twice as large, we expect also the reflection strength ROUGHLY to be twice as high
  // (not exactly, though, as the ionization changes!)
  REQUIRE(fabs(refl_energy_flux_2_refl / refl_energy_flux_1_refl - 2) < 1e-2);

  // the primary normalization should not change (it is fixed by norm_flux_cgs)
  const double prim_energy_flux_1 = refl_energy_flux_1 - refl_energy_flux_1_refl;
  const double prim_energy_flux_2 = refl_energy_flux_2 - refl_energy_flux_2_refl;
  REQUIRE(fabs(prim_energy_flux_1 / prim_energy_flux_2 - 1) < 1e-4);

}

TEST_CASE(" Test the distance parameter", "[alpha]") {
  auto default_spec = DefaultSpec(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, 3000);
  auto spec = default_spec.get_xspec_spectrum();

  LocalModel lmod_alpha(ModelName::relxilllpAlpha);
  set_default_par(lmod_alpha);
  lmod_alpha.set_par(XPar::distance, 1e5);  // distance of 100 Mpc
  lmod_alpha.set_par(XPar::mass, 1e6);      // mass of 1e6 Msolar
  lmod_alpha.set_par(XPar::refl_frac, -1.0);      // mass of 1e6 Msolar

  lmod_alpha.eval_model(spec);
  const double flux_1 = spec.get_energy_flux();
  const double lxi_1 = calculate_lxi(lmod_alpha);

  lmod_alpha.set_par(XPar::distance, 2e5);
  lmod_alpha.eval_model(spec);
  const double flux_2 = spec.get_energy_flux();
  const double lxi_2 = calculate_lxi(lmod_alpha);

  REQUIRE(fabs(lxi_1 / lxi_2 - 1) > 0.01);

}

TEST_CASE(" Test the mass parameter", "[alpha]") {
  auto default_spec = DefaultSpec(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, 3000);
  auto spec = default_spec.get_xspec_spectrum();

  LocalModel lmod_alpha(ModelName::relxilllpAlpha);
  set_default_par(lmod_alpha);
  lmod_alpha.set_par(XPar::distance, 1e5);  // distance of 100 Mpc
  lmod_alpha.set_par(XPar::mass, 1e6);      // mass of 1e6 Msolar
  lmod_alpha.set_par(XPar::refl_frac, -1.0);      // mass of 1e6 Msolar

  lmod_alpha.eval_model(spec);
  const double flux_1 = spec.get_energy_flux();
  const double lxi_1 = calculate_lxi(lmod_alpha);

  lmod_alpha.set_par(XPar::mass, 2e6);
  lmod_alpha.eval_model(spec);
  const double flux_2 = spec.get_energy_flux();
  const double lxi_2 = calculate_lxi(lmod_alpha);

  REQUIRE(fabs(lxi_1 / lxi_2 - 1) > 0.01);

  REQUIRE(fabs(flux_1 / flux_2 - 1) > 1e-6);

}

TEST_CASE(" Test the density gradient", "[alpha-test]") {
  auto default_spec = DefaultSpec(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, 3000);
  auto spec = default_spec.get_xspec_spectrum();

  LocalModel lmod_alpha(ModelName::relxilllpAlpha);
  set_default_par(lmod_alpha);
  lmod_alpha.set_par(XPar::distance, 1e5);  // distance of 100 Mpc
  lmod_alpha.set_par(XPar::mass, 1e6);      // mass of 1e6 Msolar
  lmod_alpha.set_par(XPar::logn, 15);

  const double norm_factor_ergs = 1e-11; // Flux in the 0.1-1000keV band in erg/cm^2/s
  lmod_alpha.set_par(XPar::norm_flux_cgs, norm_factor_ergs);

  // setenv("RELXILL_PRINT_DETAILS","1",1);
  lmod_alpha.set_par(XPar::refl_frac, 1.0);
  lmod_alpha.eval_model(spec);

}

TEST_CASE("CGS Flux does not change with different energy range of data", "[alpha]") {
  LocalModel lmod_alpha(ModelName::relxilllpAlpha);
  set_default_par(lmod_alpha);
  lmod_alpha.set_par(XPar::distance, 1e5);  // distance of 100 Mpc
  lmod_alpha.set_par(XPar::mass, 1e6);      // mass of 1e6 Msolar
  lmod_alpha.set_par(XPar::logn, 15);

  const double norm_factor_ergs = 1e-11; // Flux in the 0.01-1000keV band in erg/cm^2/s
  lmod_alpha.set_par(XPar::norm_flux_cgs, norm_factor_ergs);
  lmod_alpha.set_par(XPar::refl_frac, 0.0); // only the continuum

  auto default_spec_narrow =  DefaultSpec(3,10.0, 1000);
  auto spec_narrow_band = default_spec_narrow.get_xspec_spectrum();
  lmod_alpha.eval_model(spec_narrow_band);
  const double narrow_primary_energy_flux = spec_narrow_band.get_energy_flux();
  REQUIRE(fabs(narrow_primary_energy_flux - norm_factor_ergs) / norm_factor_ergs > 0.005);

  auto default_spec_narrow2 = DefaultSpec(3.0, 100.0, 400);
  auto spec_narrow_band_2 = default_spec_narrow2.get_xspec_spectrum();
  lmod_alpha.eval_model(spec_narrow_band_2);
  const double narrow_primary_energy_flux_2 = spec_narrow_band_2.get_energy_flux();
  REQUIRE(fabs(narrow_primary_energy_flux_2 - narrow_primary_energy_flux) / norm_factor_ergs > 0.005);

  // flux in the XspecSpectrum is given as integrated per bin (!)
  double const first_flux_bin = spec_narrow_band.flux[0] / (spec_narrow_band.energy[1] - spec_narrow_band.energy[0]);
  double const
      first_flux_bin_2 = spec_narrow_band_2.flux[0] / (spec_narrow_band_2.energy[1] - spec_narrow_band_2.energy[0]);

  REQUIRE(fabs(first_flux_bin - first_flux_bin_2) / first_flux_bin < 0.01);

  auto default_spec = DefaultSpec(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, 3000);
  auto spec = default_spec.get_xspec_spectrum();
  lmod_alpha.eval_model(spec);
  const double primary_energy_flux = spec.get_energy_flux();
  REQUIRE(fabs(primary_energy_flux - norm_factor_ergs) / norm_factor_ergs < 0.005);

}