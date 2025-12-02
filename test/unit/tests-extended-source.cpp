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
// #include "stdio.h"
#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "common-functions.h"
#include "Rellp_Extended.h"
// #include "Rellp.h"
// #include "Relphysics.h"

extern "C" {
// #include "writeOutfiles.h"
}

#define PREC 1e-6


TEST_CASE (
"Default relline_ring model call "
,
"[extended]"
)
 {
  int status = EXIT_SUCCESS;

  const DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relline_ring};

  const double r_src = 3.0;

  const double theta = 1.0;
  const double a = 0.9837;
  local_model.set_par(XPar::a, a);
  local_model.set_par(XPar::rin, -1);
  local_model.set_par(XPar::rout, 1000.0);
  local_model.set_par(XPar::r_src, r_src);
  local_model.set_par(XPar::theta_src, theta);
  local_model.set_par(XPar::switch_switch_returnrad, 0);
  local_model.set_par(XPar::gamma, 2.0);

  REQUIRE_NOTHROW(local_model.eval_model(spec));

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
"Default relline_slab model call "
,
"[extended]"
)
{
  int status = EXIT_SUCCESS;

  const DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relline_slab};

  const double h = 2.0;

  const double x = 0.1;
  const double a = 0.998;
  local_model.set_par(XPar::a, a);
  local_model.set_par(XPar::rin, -1);
  local_model.set_par(XPar::rout, 1000.0);
  local_model.set_par(XPar::h, h);
  local_model.set_par(XPar::x, x);
  local_model.set_par(XPar::switch_switch_returnrad, 0);
  local_model.set_par(XPar::gamma, 2.0);

  REQUIRE_NOTHROW(local_model.eval_model(spec));


  //const auto spars = get_system_parameters(local_model.rel_params(), &status);
  // save_emis_profiles(spars);

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
"Default call of relxill ring models "
,
"[extended]"
)
{

  int status = EXIT_SUCCESS;
  const DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  for (const ModelName model : {ModelName::relxill_ring_ecut, ModelName::relxill_ring}) {
    const double incl = 45.0;
    const double theta = 10.0;
    const double r_src = 2.0;
    const double spin = 0.99;
    LocalModel local_model{model};
    local_model.set_par(XPar::a, spin);
    local_model.set_par(XPar::r_src, r_src);
    local_model.set_par(XPar::theta_src, theta);
    local_model.set_par(XPar::incl, incl);

    REQUIRE_NOTHROW(local_model.eval_model(spec));
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
"Default call of relxill slab model "
,
"[extended]"
)
{

  int status = EXIT_SUCCESS;
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  double spin = 0.99;
  double h = 2.0;
  double x = 5.0;

  LocalModel local_model{ModelName::relxill_slab};
  local_model.set_par(XPar::a, spin);
  local_model.set_par(XPar::h, h);
  local_model.set_par(XPar::x, x);

  REQUIRE_NOTHROW(local_model.eval_model(spec));
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
"Comparing the reflection fraction of a ring in the lp limit "
,
"[extended]"
)
 {
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
  auto sys_par = get_system_parameters(local_model.rel_params(), &status);
  double refl = sys_par->emis->photon_fate_fractions->refl_frac;
  double frac_inf = sys_par->emis->photon_fate_fractions->f_inf;
  double frac_disk = sys_par->emis->photon_fate_fractions->f_ad;
  double frac_bhole = sys_par->emis->photon_fate_fractions->f_bh;
  printf("test lp reflfrac=%e \t fr_disk=%e \t fr_inf=%e \t fr_bh=%e \t sum=%e \t fr_disk/fr_inf=%e \n", refl,
         frac_disk, frac_inf, frac_bhole, frac_inf + frac_disk + frac_bhole, frac_disk / frac_inf);

  DefaultSpec default_spec_ext{};
  XspecSpectrum spec_ext = default_spec_ext.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_ring_ecut};


  local_model_ext.set_par(XPar::a, spin);
  local_model_ext.set_par(XPar::rin, -1.0);
  local_model_ext.set_par(XPar::rout, 1000.0);

  local_model_ext.set_par(XPar::r_src, h); // h == r_src when theta == 0
  local_model_ext.set_par(XPar::theta_src, 0.0);

  local_model_ext.eval_model(spec_ext);
  auto sys_par_ext = get_system_parameters(local_model_ext.rel_params(), &status);

  double refl_ext = sys_par_ext->emis->photon_fate_fractions->refl_frac;
  double fraction_inf = sys_par_ext->emis->photon_fate_fractions->f_inf;
  double fraction_disk = sys_par_ext->emis->photon_fate_fractions->f_ad;
  double fraction_bhole = sys_par_ext->emis->photon_fate_fractions->f_bh;

  printf("test ring reflfrac=%e \t fr_disk=%e \t fr_inf=%e \t fr_bh=%e \t sum=%e \t fr_disk/fr_inf=%e \n",
         refl_ext, fraction_disk, fraction_inf, fraction_bhole, fraction_inf + fraction_disk + fraction_bhole,
         fraction_disk / fraction_inf);

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
"Comparing the reflection fraction of a slab in the lp limit "
,
"[extended]"
)
 {
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
  auto sys_par = get_system_parameters(local_model.rel_params(), &status);
  double refl = sys_par->emis->photon_fate_fractions->refl_frac;
  double frac_inf = sys_par->emis->photon_fate_fractions->f_inf;
  double frac_disk = sys_par->emis->photon_fate_fractions->f_ad;
  double frac_bhole = sys_par->emis->photon_fate_fractions->f_bh;
  printf("test lp reflfrac=%e \t fr_disk=%e \t fr_inf=%e \t fr_bh=%e \t sum=%e \t fr_disk/fr_inf=%e \n", refl,
         frac_disk, frac_inf, frac_bhole, frac_inf + frac_disk + frac_bhole, frac_disk / frac_inf);

  DefaultSpec default_spec_ext{};
  XspecSpectrum spec_ext = default_spec_ext.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_slab};


  local_model_ext.set_par(XPar::a, spin);
  local_model_ext.set_par(XPar::rin, -1.0);
  local_model_ext.set_par(XPar::rout, 1000.0);

  local_model_ext.set_par(XPar::h, h);
  local_model_ext.set_par(XPar::x, 0.1);

  local_model_ext.eval_model(spec_ext);
  auto sys_par_ext = get_system_parameters(local_model_ext.rel_params(), &status);

  double refl_ext = sys_par_ext->emis->photon_fate_fractions->refl_frac;
  double fraction_inf = sys_par_ext->emis->photon_fate_fractions->f_inf;
  double fraction_disk = sys_par_ext->emis->photon_fate_fractions->f_ad;
  double fraction_bhole = sys_par_ext->emis->photon_fate_fractions->f_bh;

  printf("test slab reflfrac=%e \t fr_disk=%e \t fr_inf=%e \t fr_bh=%e \t sum=%e \t fr_disk/fr_inf=%e \n",
         refl_ext, fraction_disk, fraction_inf, fraction_bhole, fraction_inf + fraction_disk + fraction_bhole,
         fraction_disk / fraction_inf);

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
    "Compare lp and ring primary fluxes for zero ring radius for several parameter combinations "
,
    "[extended]"
)
{

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  DefaultSpec default_spec_ext{};

  LocalModel local_model{ModelName::relxilllp};
  auto spec = default_spec.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_ring_ecut};
  auto spec_ext = default_spec_ext.get_xspec_spectrum();
  local_model_ext.set_par(XPar::theta_src, 0.0);
  local_model_ext.set_par(XPar::refl_frac, 0);
  // with lensing it does not make sense to compare full fluxes, so we compare only reflected fluxes
  local_model.set_par(XPar::refl_frac, 0);
  double spins[3] = {0.8, 0.9, 0.998};
  double sph_radii[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for (double r_src : sph_radii){
      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::h, r_src);
      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::r_src, r_src);
      REQUIRE_NOTHROW(local_model.eval_model(spec));
      REQUIRE_NOTHROW(local_model_ext.eval_model(spec_ext));

      for (int i = 0; i < spec.num_flux_bins(); i++) {
        assert(spec.energy[i] == spec_ext.energy[i]);
        if (spec_ext.flux[i] != 0.0) {
          double current_error = 0.01;
          // as with emissivities and line profiles, we expect the difference up to current_error value now

          REQUIRE(abs(spec.flux[i] / spec_ext.flux[i] - 1.0) < current_error);
        }
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
"Compare lp and ring emissivities for zero ring size "
,
"[extended]"
)
 {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  LocalModel local_model{ModelName::relxilllp};

  DefaultSpec default_spec_ext{};
  XspecSpectrum spec_ext = default_spec_ext.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_ring_ecut};
  // relParam* rel_param_ext = local_model_ext.get_rel_params();
  local_model_ext.set_par(XPar::theta_src, 0.0);

  double spins[3] = {0.8, 0.9, 0.998};
  double sph_radii[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for(double r_src : sph_radii){

      local_model.set_par(XPar::a, spin);

      local_model.set_par(XPar::h, r_src); // h == r_src when theta == 0
      local_model.set_par(XPar::rin, -1.0);
      local_model.eval_model(spec);

      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::r_src, r_src);
      local_model_ext.set_par(XPar::rin, -1.0);
      local_model_ext.eval_model(spec_ext);
      auto sys_par = get_system_parameters(local_model.rel_params(), &status);
      auto sys_par_ext = get_system_parameters(local_model_ext.rel_params(), &status);
      int n_rad = sys_par->nr;
      double *intens = sys_par->emis->emis;
      int n_rad_ext = sys_par_ext->nr;
      double *intens_ext = sys_par_ext->emis->emis;

      assert(n_rad_ext == n_rad);
      for (int i = 0; i < n_rad; i++) {
        double current_error = 0.045;
        assert(sys_par_ext->re[i] == sys_par->re[i]); // check if the grids match

        REQUIRE(abs(intens[i] / intens_ext[i] - 1.0) < current_error);
        //printf("Intensities lp %.6e, ext %.6e, ratio %.6f at radius %.6f\n", intens[i],
        // intens_ext[i], intens[i] / intens_ext[i], sys_par_ext->re[i]);
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
"Compare lp and slab emissivities in the lp limit "
,
"[extended]"
)
{

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  LocalModel local_model{ModelName::relxilllp};

  DefaultSpec default_spec_ext{};
  XspecSpectrum spec_ext = default_spec_ext.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_slab};
  local_model_ext.set_par(XPar::x, 0.1);

  double spins[3] = {0.8, 0.9, 0.998};
  double heights[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for(double h : heights){

      local_model.set_par(XPar::a, spin);

      local_model.set_par(XPar::h, h);
      local_model.set_par(XPar::rin, -1.0);
      local_model.eval_model(spec);

      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::h, h);
      local_model_ext.set_par(XPar::rin, -1.0);
      local_model_ext.eval_model(spec_ext);
      auto sys_par = get_system_parameters(local_model.rel_params(), &status);
      auto sys_par_ext = get_system_parameters(local_model_ext.rel_params(), &status);
      int n_rad = sys_par->nr;
      double *intens = sys_par->emis->emis;
      int n_rad_ext = sys_par_ext->nr;
      double *intens_ext = sys_par_ext->emis->emis;

      assert(n_rad_ext == n_rad);
      for (int i = 0; i < n_rad; i++) {
        double current_error = 0.06;
        assert(sys_par_ext->re[i] == sys_par->re[i]); // check if the grids match

        REQUIRE(abs(intens[i] / intens_ext[i] - 1.0) < current_error);
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
"Compare lp and ring line profiles for zero ring size for several parameter combinations "
,
"[extended]"
)
 {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  DefaultSpec default_spec_ext{};
  auto spec = default_spec.get_xspec_spectrum();
  auto spec_ext = default_spec_ext.get_xspec_spectrum();
  LocalModel local_model{ModelName::relline_lp};
  LocalModel local_model_ext{ModelName::relline_ring};
  local_model_ext.set_par(XPar::theta_src, 0.0);

  double spins[3] = {0.8, 0.9, 0.998};
  double sph_radii[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for (double r_src : sph_radii){
      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::h, r_src);
      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::r_src, r_src);

      REQUIRE_NOTHROW(local_model.eval_model(spec));
      REQUIRE_NOTHROW(local_model_ext.eval_model(spec_ext));

      for (int i = 0; i < spec.num_flux_bins(); i++) {
        if (spec.flux[i] > 0.0) {
          double current_error_tail = 0.042;
          // to remove zero division
          REQUIRE(abs(spec.flux[i] - spec_ext.flux[i]) / spec.flux[i] < current_error_tail);
          if (spec.energy[i] >= 1.0 and spec.energy[i] <= 8.0) {
            double current_error_peak = 0.002;
            // this is the most important range of energies
            REQUIRE(abs(spec.flux[i] - spec_ext.flux[i]) / spec.flux[i] < current_error_peak);
          }
        }
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}

TEST_CASE (
"Compare lp and slab line profiles in the lp limit for several parameter combinations "
,
"[extended]"
)
{

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  DefaultSpec default_spec_ext{};
  auto spec = default_spec.get_xspec_spectrum();
  auto spec_ext = default_spec_ext.get_xspec_spectrum();
  LocalModel local_model{ModelName::relline_lp};
  LocalModel local_model_ext{ModelName::relline_slab};
  local_model_ext.set_par(XPar::x, 0.1);

  double spins[3] = {0.8, 0.9, 0.998};
  double heights[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for (double h : heights){
      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::h, h);
      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::h, h);

      REQUIRE_NOTHROW(local_model.eval_model(spec));
      REQUIRE_NOTHROW(local_model_ext.eval_model(spec_ext));

      for (int i = 0; i < spec.num_flux_bins(); i++) {
        if (spec.flux[i] > 0.0) {
          double current_error_tail = 0.06;
          // to remove zero division
          REQUIRE(abs(spec.flux[i] - spec_ext.flux[i]) / spec.flux[i] < current_error_tail);
          if (spec.energy[i] >= 1.0 and spec.energy[i] <= 8.0) {
            double current_error_peak = 0.003;
            // this is the most important range of energies
            REQUIRE(abs(spec.flux[i] - spec_ext.flux[i]) / spec.flux[i] < current_error_peak);
          }
        }
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
"Compare lp and ring reflected fluxes for zero ring radius for several parameter combinations "
,
"[extended]"
)
 {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  DefaultSpec default_spec_ext{};

  LocalModel local_model{ModelName::relxilllp};
  auto spec = default_spec.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_ring_ecut};
  auto spec_ext = default_spec_ext.get_xspec_spectrum();
  local_model_ext.set_par(XPar::theta_src, 0.0);
  local_model_ext.set_par(XPar::refl_frac, -1);
  // with lensing it does not make sense to compare full fluxes, so we compare only reflected fluxes
  local_model.set_par(XPar::refl_frac, -1);
  double spins[3] = {0.8, 0.9, 0.998};
  double sph_radii[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for (double r_src : sph_radii){
      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::h, r_src);
      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::r_src, r_src);
      REQUIRE_NOTHROW(local_model.eval_model(spec));
      REQUIRE_NOTHROW(local_model_ext.eval_model(spec_ext));

      for (int i = 0; i < spec.num_flux_bins(); i++) {
        assert(spec.energy[i] == spec_ext.energy[i]);
        if (spec_ext.flux[i] != 0.0) {
          double current_error = 0.047;
          // as with emissivities and line profiles, we expect the difference up to current_error value now
          REQUIRE(abs(spec.flux[i] / spec_ext.flux[i] - 1.0) < current_error);
        }
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
"Compare lp and slab reflected fluxes in the lp limit for several parameter combinations "
,
"[extended]"
)
{

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  DefaultSpec default_spec_ext{};

  LocalModel local_model{ModelName::relxilllpCp};
  auto spec = default_spec.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxill_slab};
  auto spec_ext = default_spec_ext.get_xspec_spectrum();
  local_model_ext.set_par(XPar::x, 0.1);
  local_model_ext.set_par(XPar::refl_frac, -1);
  // with lensing it does not make sense to compare full fluxes, so we compare only reflected fluxes
  local_model.set_par(XPar::refl_frac, -1);
  double spins[3] = {0.8, 0.9, 0.998};
  double heights[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for (double h : heights){
      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::h, h);
      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::h, h);
      REQUIRE_NOTHROW(local_model.eval_model(spec));
      REQUIRE_NOTHROW(local_model_ext.eval_model(spec_ext));

      for (int i = 0; i < spec.num_flux_bins(); i++) {
        assert(spec.energy[i] == spec_ext.energy[i]);
        if (spec_ext.flux[i] != 0.0) {
          double current_error_high_en = 0.07; // high energies (>100 keV)
          double current_error_low_en = 0.01; // else
          if (spec.energy[i] > 100.0) {
            // a high energy tail is different for slab
            REQUIRE(abs(spec.flux[i] / spec_ext.flux[i] - 1.0) < current_error_high_en);
          } else {
            REQUIRE(abs(spec.flux[i] / spec_ext.flux[i] - 1.0) < current_error_low_en);
          }

        }
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE (
"Setting input parameters for ring source outside of the allowed range "
,
"[extended]"
)
 {

  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxill_ring);

  auto spec = default_spec.get_xspec_spectrum();

  double radius_below_horizon = 0.5;
  lmod.set_par(XPar::r_src, radius_below_horizon);
  lmod.set_par(XPar::theta_src, 5.0); // the source is inside BH and not on-axis

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

  double negative_angle = -10.0;
  lmod.set_par(XPar::theta_src, negative_angle);

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

}


TEST_CASE (
"Setting negative slab radius "
,
"[extended]"
)
{

  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxill_slab);

  auto spec = default_spec.get_xspec_spectrum();

  double negative_x = -5.0;
  lmod.set_par(XPar::x, negative_x);

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

}


TEST_CASE (
"Setting slab source inside of EH "
,
"[extended]"
)
{

  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxill_slab);

  auto spec = default_spec.get_xspec_spectrum();

  double h_below_horizon = 0.9;
  lmod.set_par(XPar::x, 3.0);
  lmod.set_par(XPar::h, h_below_horizon);

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

}


TEST_CASE (
"Setting outer edge of slab source outside table limits (extrapolation) "
,
"[extended]"
)
{

  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxill_slab);

  auto spec = default_spec.get_xspec_spectrum();

  double low_height = 0.5;
  lmod.set_par(XPar::x, 30.0);
  lmod.set_par(XPar::h, low_height);

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

}
