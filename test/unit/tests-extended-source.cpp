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

#include "cfitsio.h"
#include "stdio.h"
#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "common-functions.h"
#include "Rellp_Extended.h"
#include "Rellp.h"
#include "Relphysics.h"

extern "C" {
#include "writeOutfiles.h"
}

#define PREC 1e-6



TEST_CASE("Default line model call  Test", "[extended]") {

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relline_ext};

  double height = 3.0;
  double r_ring = 5.0;
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::x, r_ring);

  local_model.eval_model(spec);
/*
 *  ONLY FOR LOCAL USAGE (to save fluxes or emissivities to files)
 *
  int n_bins = spec.num_flux_bins();
  double *flux = spec.flux;
  double *energy = spec.energy;

  string strname = "/home/nekrasov/git/tables/test-ext-line-h" + std::to_string(height)
  + "-d" + std::to_string(dist) + ".dat";
  char *foutName = new char[999];
  strcpy(foutName, strname.c_str());
  FILE *fp = fopen(foutName, "w+");
  if (fp != nullptr) {
    for (int ii = 0; ii < n_bins; ii++) {
      fprintf(fp, " %e \t %e \n", energy[ii], flux[ii]);
    }
  }
  fclose(fp);
*/
}


TEST_CASE("Xspec line model call Test", "[extended]") {

  int status = EXIT_SUCCESS;

  LocalModel local_model{ModelName::relline_ext};
  local_model.set_par(XPar::x, 5.0);
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  local_model.eval_model(spec);

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Default relxill model call Test", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxill_ext_ecut};
  double height = 5.0;
  double r_ring = 5.0;
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::x, r_ring);

  local_model.eval_model(spec);

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Default relxill Cp ext model call ", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxill_ext};
  double spin = 0.998;
  double height = 2.0;
  double r_ring = 0.0;
  double incl = 30.0;
  local_model.set_par(XPar::a, spin);
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::x, r_ring);
  local_model.set_par(XPar::incl, incl);

  local_model.eval_model(spec);

  REQUIRE(status == EXIT_SUCCESS);
}

TEST_CASE("Default call of relxill ext model with disk geometry", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxill_ext};
  double spin = 0.998;
  double height = 2.0;
  double slab_extent = 3.0;
  double incl = 30.0;
  local_model.set_par(XPar::a, spin);
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::x, slab_extent);
  local_model.set_par(XPar::incl, incl);
  local_model.set_par(XPar::switch_switch_source_geometry, 2);

  local_model.eval_model(spec);

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Checking the reflection fraction", "[extended]") {
  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  LocalModel local_model{ModelName::relxilllp};

  double spin = 0.998;
  double h = 10.1;
  local_model.set_par(XPar::a, spin);
  local_model.set_par(XPar::rin, -1.0);
  local_model.set_par(XPar::rout, 1000.0);
  local_model.set_par(XPar::h, h);

  local_model.eval_model(spec);
  relParam* rel_param = local_model.get_rel_params();
  auto sys_par = get_system_parameters(rel_param, &status);
  double refl = sys_par->emis->photon_fate_fractions->refl_frac;
  double frac_inf = sys_par->emis->photon_fate_fractions->f_inf;
  double frac_disk = sys_par->emis->photon_fate_fractions->f_ad;
  double frac_bhole = sys_par->emis->photon_fate_fractions->f_bh;
  printf("test lp reflfrac=%e \t fr_disk=%e \t fr_inf=%e \t fr_bh=%e \t sum=%e \t fr_disk/fr_inf=%e \n", refl,
         frac_disk, frac_inf, frac_bhole, frac_inf + frac_disk + frac_bhole, frac_disk / frac_inf);

  DefaultSpec default_spec_ext{};
  XspecSpectrum spec_ext = default_spec_ext.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_ext_ecut};

  local_model_ext.set_par(XPar::a, spin);
  local_model_ext.set_par(XPar::rin, -1.0);
  local_model_ext.set_par(XPar::rout, 1000.0);
  local_model_ext.set_par(XPar::h, h);
  local_model_ext.set_par(XPar::x, 0.0);
  local_model_ext.set_par(XPar::switch_switch_source_geometry, 1);

  local_model_ext.eval_model(spec_ext);
  relParam* rel_param_ext = local_model_ext.get_rel_params();
  auto sys_par_ext = get_system_parameters(rel_param_ext, &status);

  double refl_ext = sys_par_ext->emis->photon_fate_fractions->refl_frac;
  double fraction_inf = sys_par_ext->emis->photon_fate_fractions->f_inf;
  double fraction_disk = sys_par_ext->emis->photon_fate_fractions->f_ad;
  double fraction_bhole = sys_par_ext->emis->photon_fate_fractions->f_bh;
  //printf("%e \n", rel_param_ext->rout);
  printf("test ext reflfrac=%e \t fr_disk=%e \t fr_inf=%e \t fr_bh=%e \t sum=%e \t fr_disk/fr_inf=%e \n",
         refl_ext, fraction_disk, fraction_inf, fraction_bhole, fraction_inf + fraction_disk + fraction_bhole,
         fraction_disk / fraction_inf);

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Test to compare lp and ext emissivities", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  LocalModel local_model{ModelName::relxilllp};
  DefaultSpec default_spec_ext{};
  XspecSpectrum spec_ext = default_spec_ext.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_ext_ecut};
  // relParam* rel_param_ext = local_model_ext.get_rel_params();
  local_model_ext.set_par(XPar::x, 0.0);  // for this test we compare the lp emissivity and ext for ring=0

  double spins[3] = {0.0, 0.7, 0.998};
  double heights[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for(double height : heights){

      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::h, height);
      local_model.set_par(XPar::rin, -1.0);
      local_model.eval_model(spec);

      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::h, height);
      local_model_ext.set_par(XPar::rin, -1.0);
      local_model_ext.eval_model(spec_ext);
      relParam* rel_param = local_model.get_rel_params();
      relParam* rel_param_ext = local_model_ext.get_rel_params();
      auto sys_par = get_system_parameters(rel_param, &status);
      auto sys_par_ext = get_system_parameters(rel_param_ext, &status);
      int n_rad = sys_par->nr;
      double *intens = sys_par->emis->emis;
      int n_rad_ext = sys_par_ext->nr;
      double *intens_ext = sys_par_ext->emis->emis;

      assert(n_rad_ext == n_rad);
      for (int i = 0; i < n_rad; i++) {
        assert(sys_par_ext->re[i] == sys_par->re[i]); // check if the grids match
        // check if the emissivities are the same. Deviations of an order of a few % expected. Up to 10% in the most extreme cases
        // (h at EH)
        REQUIRE(abs(intens[i] / intens_ext[i] - 1.0) < 0.033);
        //printf("Intensities lp %.6e, ext %.6e, ratio %.6f at radius %.6f\n", intens[i],
        // intens_ext[i], intens[i] / intens_ext[i], sys_par_ext->re[i]);

      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}

TEST_CASE("Test to compare lp and ext line profiles for several spins and heights ", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  DefaultSpec default_spec_ext{};
  auto spec = default_spec.get_xspec_spectrum(); // when it is outside of "for" - energy shifts wrong
  auto spec_ext = default_spec_ext.get_xspec_spectrum(); // or is it?
  LocalModel local_model{ModelName::relline_lp};
  LocalModel local_model_ext{ModelName::relline_ext};
  local_model_ext.set_par(XPar::x, 0.0);

  double spins[3] = {0.0, 0.5, 0.998};
  double heights[3] = {2.2, 5.0, 10.0};

  for (double spin : spins){
    for (double height : heights){
      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::h, height);
      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::h, height);

      REQUIRE_NOTHROW(local_model.eval_model(spec));
      REQUIRE_NOTHROW(local_model_ext.eval_model(spec_ext));

      for (int i = 0; i < spec.num_flux_bins(); i++) {
        // to remove zero division
        if (spec.flux[i] > 0.0) {
          // currently due to slightly different energy shifts we can only get down to ~8% difference for these a,h
          REQUIRE(abs(spec.flux[i] - spec_ext.flux[i]) / spec.flux[i] < 0.075);
          // although these few % diff are only for the red tail of the line, on average it is below 1%
          if (spec.energy[i] >= 1.0 and spec.energy[i] <= 8.0) {
            REQUIRE(abs(spec.flux[i] - spec_ext.flux[i]) / spec.flux[i] < 0.003);
          }
        }
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Test to compare lp and ext fluxes for several spins and heights ", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  DefaultSpec default_spec_ext{};

  LocalModel local_model{ModelName::relxilllp};
  auto spec = default_spec.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_ext_ecut};
  auto spec_ext = default_spec_ext.get_xspec_spectrum();
  local_model_ext.set_par(XPar::x, 0.0);
  local_model_ext.set_par(XPar::refl_frac, -1);
  // with lensing it does not make sense to compare full fluxes, so we compare only reflected fluxes
  local_model.set_par(XPar::refl_frac, -1);
  double spins[3] = {0.0, 0.5, 0.998};
  double heights[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for (double height : heights){
      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::h, height);
      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::h, height);

      REQUIRE_NOTHROW(local_model.eval_model(spec));
      REQUIRE_NOTHROW(local_model_ext.eval_model(spec_ext));

      for (int i = 0; i < spec.num_flux_bins(); i++) {
        assert(spec.energy[i] == spec_ext.energy[i]);
        if (spec_ext.flux[i] != 0.0){
          // as with emissivities and line profiles, we expect the difference up to ~6% now
          REQUIRE(abs(spec.flux[i] / spec_ext.flux[i] - 1.0) < 0.06);
        }
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Test setting input parameters for extended source outside the allowed range","[extended]") {

  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxill_ext_ecut);

  auto spec = default_spec.get_xspec_spectrum();

  double height_below_horizon = 0.9;
  lmod.set_par(XPar::h, height_below_horizon);
  lmod.set_par(XPar::x, 0.01);

  REQUIRE_NOTHROW(lmod.eval_model(spec));

  // require a real
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

}

