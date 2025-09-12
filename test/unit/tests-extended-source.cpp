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



TEST_CASE("Default line model call ", "[extended]") {
  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relline_ext};

  double height = 3.0;

  double x = 1.0;
  double a = 0.9837;
  local_model.set_par(XPar::a, a);
  local_model.set_par(XPar::rin, -1);
  local_model.set_par(XPar::rout, 1000.0);
  local_model.set_par(XPar::r, height);
  local_model.set_par(XPar::theta, x);
  local_model.set_par(XPar::switch_switch_returnrad, 0);
  local_model.set_par(XPar::gamma, 2.0);

  // slab-specific
  //local_model.set_par(XPar::x_in, 0.0);
  //local_model.set_par(XPar::switch_switch_source_geometry, 2);

  REQUIRE_NOTHROW(local_model.eval_model(spec));
  // local_model.eval_model(spec);

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Default call of relxill ext models for different geometries ", "[extended]") {

  //int status = EXIT_SUCCESS;
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  // LocalModel local_model{ModelName::relxill_ext};
  // do some dummy param sets just to try that it works
  double spin = 0.99;
  double height = 2.0;
  double x = 10.0;
  double incl = 45.0;
  for (ModelName model : {ModelName::relline_ext, ModelName::relxill_ext_ecut, ModelName::relxill_ext}) {
    for (int type_geom : {1}) {
      LocalModel local_model{model};
      // local_model.set_par(XPar::switch_switch_source_geometry, type_geom);
      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::r, height);
      local_model.set_par(XPar::theta, x);
      local_model.set_par(XPar::incl, incl);

      REQUIRE_NOTHROW(local_model.eval_model(spec));
    }
  }
  //REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Checking the reflection fraction in the lp limit", "[extended]") {
  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  LocalModel local_model{ModelName::relxilllp};

  double spin = 0.998;
  double h = 3.0;
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
  local_model_ext.set_par(XPar::r, h);
  local_model_ext.set_par(XPar::theta, 0.0);
  //local_model_ext.set_par(XPar::switch_switch_source_geometry, 1);

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


// Does not work for spherical relxill, need modifications

TEST_CASE("Compare lp and ext emissivities for x = 0", "[extended]") {

  int status = EXIT_SUCCESS;
  double current_error = 0.045;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  LocalModel local_model{ModelName::relxilllp};

  DefaultSpec default_spec_ext{};
  XspecSpectrum spec_ext = default_spec_ext.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_ext_ecut};
  // relParam* rel_param_ext = local_model_ext.get_rel_params();
  local_model_ext.set_par(XPar::theta, 0.0);  // for this test we compare the lp emissivity and ext for ring=0

  double spins[3] = {0.8, 0.9, 0.998};
  double heights[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for(double height : heights){

      local_model.set_par(XPar::a, spin);

      local_model.set_par(XPar::h, height);
      local_model.set_par(XPar::rin, -1.0);
      local_model.eval_model(spec);

      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::r, height);
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
        // check if the emissivities are the same. Deviations of an order of a few % expected.
        // Up to 10% in the most extreme cases

        REQUIRE(abs(intens[i] / intens_ext[i] - 1.0) < current_error);
        //printf("Intensities lp %.6e, ext %.6e, ratio %.6f at radius %.6f\n", intens[i],
        // intens_ext[i], intens[i] / intens_ext[i], sys_par_ext->re[i]);
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}

TEST_CASE("Compare lp and ext line profiles for x = 0 for several spins and heights ", "[extended]") {

  int status = EXIT_SUCCESS;
  double current_error_peak = 0.002;
  double current_error_tail = 0.042; // almost same as for emissivity btw, which makes sense

  DefaultSpec default_spec{};
  DefaultSpec default_spec_ext{};
  auto spec = default_spec.get_xspec_spectrum(); // when it is outside of "for" - energy shifts wrong
  auto spec_ext = default_spec_ext.get_xspec_spectrum(); // or is it?
  LocalModel local_model{ModelName::relline_lp};
  LocalModel local_model_ext{ModelName::relline_ext};
  local_model_ext.set_par(XPar::theta, 0.0);

  double spins[3] = {0.8, 0.9, 0.998};
  double heights[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for (double height : heights){
      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::h, height);
      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::r, height);

      REQUIRE_NOTHROW(local_model.eval_model(spec));
      REQUIRE_NOTHROW(local_model_ext.eval_model(spec_ext));

      for (int i = 0; i < spec.num_flux_bins(); i++) {
        // to remove zero division
        if (spec.flux[i] > 0.0) {
          // currently due to slightly different energy shifts we can only get down to ~4% difference for these a,h
          REQUIRE(abs(spec.flux[i] - spec_ext.flux[i]) / spec.flux[i] < current_error_tail);
          // although these few % diff are only for the red tail of the line, on average it is below 1%
          if (spec.energy[i] >= 1.0 and spec.energy[i] <= 8.0) {
            REQUIRE(abs(spec.flux[i] - spec_ext.flux[i]) / spec.flux[i] < current_error_peak);
          }
        }
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Compare lp and ext reflected fluxes for x = 0 for several spins and heights ", "[extended]") {

  int status = EXIT_SUCCESS;
  double current_error = 0.047;

  DefaultSpec default_spec{};
  DefaultSpec default_spec_ext{};

  LocalModel local_model{ModelName::relxilllp};
  auto spec = default_spec.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_ext_ecut};
  auto spec_ext = default_spec_ext.get_xspec_spectrum();
  local_model_ext.set_par(XPar::theta, 0.0);
  local_model_ext.set_par(XPar::refl_frac, -1);
  // with lensing it does not make sense to compare full fluxes, so we compare only reflected fluxes
  local_model.set_par(XPar::refl_frac, -1);
  double spins[3] = {0.8, 0.9, 0.998};
  double heights[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for (double height : heights){
      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::h, height);
      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::r, height);

      REQUIRE_NOTHROW(local_model.eval_model(spec));
      REQUIRE_NOTHROW(local_model_ext.eval_model(spec_ext));

      for (int i = 0; i < spec.num_flux_bins(); i++) {
        assert(spec.energy[i] == spec_ext.energy[i]);
        if (spec_ext.flux[i] != 0.0){
          // as with emissivities and line profiles, we expect the difference up to current_error value now
          REQUIRE(abs(spec.flux[i] / spec_ext.flux[i] - 1.0) < current_error);
        }
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}

TEST_CASE("Setting input parameters for extended source outside of the allowed range ","[extended]") {

  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxill_ext);

  auto spec = default_spec.get_xspec_spectrum();

  // height below EH (will be replaced with spherical radius eventually)
  double height_below_horizon = 0.9;
  lmod.set_par(XPar::r, height_below_horizon);
  lmod.set_par(XPar::theta, 0.1); // so the source is inside of BH

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

  double negative_x = -10.0;
  // lmod.set_par(XPar::h, 3.0);
  lmod.set_par(XPar::theta, negative_x);

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

  /*
  // and for slab source
  lmod.set_par(XPar::switch_switch_source_geometry, 2);
  double too_small_x = 0.0;
  // lmod.set_par(XPar::h, 3.0);
  lmod.set_par(XPar::x, too_small_x);

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

  double who_missed_inequalities_lesson_in_school_math = 10.0;
  lmod.set_par(XPar::x_in, who_missed_inequalities_lesson_in_school_math);

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);
  */
}

