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

extern "C" {
#include "writeOutfiles.h"
}

#define PREC 1e-6



TEST_CASE("Default line model call  Test", "[extended]") {

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relline_lp_ext};

  double height = 3.0;
  double r_ring = 5.0;
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::r_ring, r_ring);

  local_model.eval_model(spec);
/*
 *  ONLY FOR LOCAL USAGE
 *
  int n_bins = spec.num_flux_bins();
  double *flux = spec.flux;
  double *energy = spec.energy;

  string strname = "/home/nekrasov/git/tables/test-ext-line-h" + std::to_string(height) + "-d" + std::to_string(dist) + ".dat";
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

  LocalModel local_model{ModelName::relline_lp_ext};
  local_model.set_par(XPar::r_ring, 5.0);
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  local_model.eval_model(spec);

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Default relxill model call Test", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllp_ext};
  double height = 3.0;
  double r_ring = 5.0;
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::r_ring, r_ring);

  local_model.eval_model(spec);

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Default relxill Cp ext model call ", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllpCp_ext};
  double height = 3.0;
  double r_ring = 5.0;
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::r_ring, r_ring);

  local_model.eval_model(spec);

  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Checking the reflection fraction", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllp};
  relParam* rel_param = local_model.get_rel_params();
  auto sys_par = get_system_parameters(rel_param, &status);

  double spin = 0.9;
  double h = 5.0;
  rel_param->a = spin;
  rel_param->height = h;
  //rel_param->rout= 1000.0;
  lpTable *cached_lp_table = get_lp_table(&status);
  calc_emis_jet_point_source(sys_par->emis, rel_param, rel_param->height, rel_param->beta,
                             cached_lp_table, &status);

  local_model.eval_model(spec);

  double refl = sys_par->emis->photon_fate_fractions->refl_frac;
  double frac_inf = sys_par->emis->photon_fate_fractions->f_inf;
  double frac_disk = sys_par->emis->photon_fate_fractions->f_ad;
  double frac_bhole = sys_par->emis->photon_fate_fractions->f_bh;
  //printf("%e \n", rel_param->rout);
  printf("test lp reflfrac=%e \t fr_disk=%e \t fr_inf=%e \t fr_bh=%e \t sum=%e \t fr_disk/fr_inf=%e \n", refl,
         frac_disk, frac_inf, frac_bhole, frac_inf + frac_disk + frac_bhole, frac_disk / frac_inf);

  DefaultSpec default_spec_ext{};
  XspecSpectrum spec_ext = default_spec_ext.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxilllp_ext};
  relParam* rel_param_ext = local_model_ext.get_rel_params();
  auto sys_par_ext = get_system_parameters(rel_param_ext, &status);

  rel_param_ext->a = spin;
  //rel_param_ext->height = 2.8962; // weird that rfrac coincides with lp only when I set slightly lower height than lp
  rel_param_ext->height = h;
  rel_param_ext->r_ring = 0.0;
  //rel_param_ext->rout= 1000.0;
  get_emis_extendedSource(sys_par_ext->emis, rel_param_ext, &status);

  local_model_ext.eval_model(spec_ext);

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
  relParam* rel_param = local_model.get_rel_params();

  LocalModel local_model_ext{ModelName::relxilllp_ext};
  relParam* rel_param_ext = local_model_ext.get_rel_params();
  rel_param_ext->r_ring = 0.0; // for this test we compare the lp emissivity and ext for ring=0

  double spins[3] = {-0.998, 0.0, 0.998};
  double heights[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for(double height : heights){
      rel_param->a = spin;
      rel_param->height = height;
      rel_param->rin = kerr_rms(spin);
      local_model.eval_model(spec);

      rel_param_ext->a = spin;
      rel_param_ext->height = height;
      rel_param_ext->rin = kerr_rms(spin);
      local_model_ext.eval_model(spec);
      auto sys_par = get_system_parameters(rel_param, &status);
      auto sys_par_ext = get_system_parameters(rel_param_ext, &status);
      int n_rad = sys_par->nr;
      double *rad = sys_par->re;
      double *intens = sys_par->emis->emis;
      int n_rad_ext = sys_par_ext->nr;
      double *rad_ext = sys_par_ext->re;
      double *intens_ext = sys_par_ext->emis->emis;

      assert(n_rad_ext == n_rad);
      assert(rad_ext == rad);
      for (int i = 0; i <= n_rad; i++) {
        assert(rad_ext[i] == rad[i]); // check if the grids are the same
        assert(intens[i] == intens_ext[i]); // check if the emissivities are the same
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Test to compare lp and ext fluxes for several spins and heights ", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};

  LocalModel local_model{ModelName::relxilllp};
  auto spec = default_spec.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxilllp_ext};
  auto spec_ext = default_spec.get_xspec_spectrum();
  local_model_ext.set_par(XPar::r_ring, 0.0);

  double spins[3] = {0.0, 0.5, 0.998};
  double heights[3] = {3.0, 5.0, 10.0};

  for (double spin : spins){
    for(double height : heights){
      local_model.set_par(XPar::a, spin);
      local_model.set_par(XPar::h, height);
      REQUIRE_NOTHROW(local_model.eval_model(spec));
      local_model_ext.set_par(XPar::a, spin);
      local_model_ext.set_par(XPar::h, height);
      REQUIRE_NOTHROW(local_model_ext.eval_model(spec_ext));
      for (int i = 0; i < spec.num_flux_bins(); i++) {
        REQUIRE(abs(spec.flux[i] - spec_ext.flux[i]) < PREC);
      }
    }
  }
  REQUIRE(status == EXIT_SUCCESS);
}


TEST_CASE("Test setting input parameters for extended source outside the allowed range","[extended]") {

  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxilllp_ext);

  auto spec = default_spec.get_xspec_spectrum();

  double height_below_horizon = 0.9;
  lmod.set_par(XPar::h, height_below_horizon);
  lmod.set_par(XPar::r_ring, 0.01);

  REQUIRE_NOTHROW(lmod.eval_model(spec));

  // require a real
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

}

