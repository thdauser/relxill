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

 // int status = EXIT_SUCCESS; // in fact this is a bool variable to check the program status

DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relline_lp_ext};

  double height = 2.0;
  double dist = 20.0;
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::d_ring, dist);

  local_model.eval_model(spec);
/*
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

  // int status = EXIT_SUCCESS;

  LocalModel local_model{ModelName::relline_lp_ext};
  local_model.set_par(XPar::d_ring, 5.0);
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  local_model.eval_model(spec);

}


TEST_CASE("Default relxill model call  Test", "[extended]") {

  // int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllp_ext};
  double height = 2.0;
  double dist = 20.0;
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::d_ring, dist);

  local_model.eval_model(spec);
  /*
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



TEST_CASE("Test to compare lp and ext", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  LocalModel local_model{ModelName::relxilllp};
  relParam* rel_param = local_model.get_rel_params();
  auto sys_par = get_system_parameters(rel_param, &status);

  // (rel_param).distance = 2.0;
  rel_param->a = 0.998;
  rel_param->height = 3.0;
  rel_param->emis1 = 3.0;
  rel_param->emis2 = 3.0;
  rel_param->rbr = 30.0;
  rel_param->lineE = 1.0;
  rel_param->return_rad = 0.0;
  lpTable *cached_lp_table = get_lp_table(&status);
  calc_emis_jet_point_source(sys_par->emis, rel_param, rel_param->height, rel_param->beta, cached_lp_table, &status);
  local_model.eval_model(spec);

  LocalModel local_model_ext{ModelName::relxilllp_ext};
  relParam* rel_param_ext = local_model_ext.get_rel_params();
  auto sys_par_ext = get_system_parameters(rel_param_ext, &status);

  // (rel_param).distance = 2.0;
  rel_param_ext->a = 0.998;
  rel_param_ext->height = 3.0;
  rel_param_ext->d_ring = 0.1;
  rel_param_ext->emis1 = 3.0;
  rel_param_ext->emis2 = 3.0;
  rel_param_ext->rbr = 30.0;
  rel_param_ext->lineE = 1.0;
  rel_param_ext->return_rad = 0.0;
  get_emis_extendedSource(sys_par_ext->emis, rel_param_ext, &status);
  local_model_ext.eval_model(spec);

  /*
  int n_rad = sys_par->nr;
  double *rad = sys_par->re;
  double *intens = sys_par->emis->emis;
  int n_rad_ext = sys_par_ext->nr;
  double *rad_ext = sys_par_ext->re;
  double *intens_ext = sys_par_ext->emis->emis;

  double refl_ext = sys_par_ext->emis->photon_fate_fractions->refl_frac;

  string strname = "/home/nekrasov/analysis/test-ext-emis-d" + std::to_string(rel_param_ext->d_ring) + ".dat";
  char *foutName = new char[999];
  strcpy(foutName, strname.c_str());
  FILE *fp = fopen(foutName, "w+");
  if (fp != nullptr) {
    for (int ii = 0; ii < n_rad; ii++) {
      fprintf(fp, " %e \t %e \t %e \t %e \t %e \n", rad[ii], intens[ii], rad_ext[ii], intens_ext[ii], intens[ii] / intens_ext[ii] / intens[200] * intens_ext[200]);
    }
  }
  fclose(fp);
  */
}



TEST_CASE("Checking the reflection fraction", "[extended]") {

  int status = EXIT_SUCCESS; // in fact this is a bool variable to check the program status

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum(); // afaiu it gives here a default energy, flux values and number of bins used

  LocalModel local_model{ModelName::relxilllp};
  relParam* rel_param = local_model.get_rel_params(); // getting parameters for the given model...
  auto sys_par = get_system_parameters(rel_param, &status); // also some parameters

  rel_param->a = 0.998;
  rel_param->height = 2.0;
 // rel_param->rout= 1000.0;
  lpTable *cached_lp_table = get_lp_table(&status);
  calc_emis_jet_point_source(sys_par->emis, rel_param, rel_param->height, rel_param->beta, cached_lp_table, &status);

  local_model.eval_model(spec);

  double refl = sys_par->emis->photon_fate_fractions->refl_frac;
  double frac_inf = sys_par->emis->photon_fate_fractions->f_inf;
  double frac_disk = sys_par->emis->photon_fate_fractions->f_ad;
  double frac_bhole = sys_par->emis->photon_fate_fractions->f_bh;
  printf("%e \n", rel_param->rout);
  printf("test reflfrac=%e \t fr_disk=%e \t fr_inf=%e \t fr_bh=%e \t sum=%e \t fr_disk/fr_inf=%e \n", refl, frac_disk, frac_inf, frac_bhole, frac_inf + frac_disk + frac_bhole, frac_disk / frac_inf);

  DefaultSpec default_spec_ext{};
  XspecSpectrum spec_ext = default_spec_ext.get_xspec_spectrum();
  LocalModel local_model_ext{ModelName::relxilllp_ext};
  relParam* rel_param_ext = local_model_ext.get_rel_params();
  auto sys_par_ext = get_system_parameters(rel_param_ext, &status);

  rel_param_ext->a = 0.998;
  rel_param_ext->height = 2.0;
  rel_param_ext->d_ring = 1.01; // setting some manual distance
 // rel_param_ext->rout= 1000.0;
  get_emis_extendedSource(sys_par_ext->emis, rel_param_ext, &status);

  local_model_ext.eval_model(spec_ext);

  double refl_ext = sys_par_ext->emis->photon_fate_fractions->refl_frac;
  double fraction_inf = sys_par_ext->emis->photon_fate_fractions->f_inf;
  double fraction_disk = sys_par_ext->emis->photon_fate_fractions->f_ad;
  double fraction_bhole = sys_par_ext->emis->photon_fate_fractions->f_bh;
  printf("%e \n", rel_param_ext->rout);
  printf("ext test reflfrac=%e \t fr_disk=%e \t fr_inf=%e \t fr_bh=%e \t sum=%e \t fr_disk/fr_inf=%e \n", refl_ext, fraction_disk, fraction_inf, fraction_bhole, fraction_inf + fraction_disk + fraction_bhole, fraction_disk / fraction_inf);

}



TEST_CASE("Test to compare lp and ext emissivities", "[extended]") {

  int status = EXIT_SUCCESS;

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  LocalModel local_model{ModelName::relxilllp};
  relParam* rel_param = local_model.get_rel_params();
  auto sys_par = get_system_parameters(rel_param, &status);

  rel_param->a = 0.998;
  rel_param->height = 30.0;
  //rel_param->rout = 1000.0;
  //rel_param->emis1 = 3.0;
  //rel_param->emis2 = 3.0;
  //rel_param->rbr = 30.0;
  //rel_param->lineE = 1.0;
  //rel_param->return_rad = 0.0;
  //get_emis_jet(sys_par->emis, rel_param, &status);
  local_model.eval_model(spec);


  LocalModel local_model_ext{ModelName::relxilllp_ext};
  relParam* rel_param_ext = local_model_ext.get_rel_params();
  auto sys_par_ext = get_system_parameters(rel_param_ext, &status);

  rel_param_ext->a = 0.998;
  rel_param_ext->height = 30.0;
  rel_param_ext->d_ring = 0.0;
  //rel_param_ext->rin = kerr_rms(rel_param_ext->a);
  //rel_param_ext->rout = 1000.0;
  //rel_param_ext->emis1 = 3.0;
  //rel_param_ext->emis2 = 3.0;
  //rel_param_ext->rbr = 30.0;
  //rel_param_ext->lineE = 1.0;
  //rel_param_ext->return_rad = 0.0;

  //get_emis_jet(sys_par_ext->emis, rel_param_ext, &status);
  local_model_ext.eval_model(spec);
  sys_par = get_system_parameters(rel_param, &status);
  sys_par_ext = get_system_parameters(rel_param_ext, &status);
  int n_rad = sys_par->nr;
  double *rad = sys_par->re;
  double *intens = sys_par->emis->emis;
  int n_rad_ext = sys_par_ext->nr;
  //double *rad_ext = sys_par_ext->re;
  double *intens_ext = sys_par_ext->emis->emis;

  double sum = 0.0;
  for (int i = 0; i <= n_rad; i++) {
    sum += intens[i];
  }
  printf("sum_emis = %e \n", sum);

  double sum_ext = 0.0;
  for (int i = 0; i <= n_rad_ext; i++) {
    sum_ext += intens_ext[i];
  }
  printf("sum_emis_ext = %e \n", sum_ext);
}

TEST_CASE("Test to compare lp and ext emissivities for several heights ", "[extended]") {
  int status = EXIT_SUCCESS;
  DefaultSpec default_spec{};

  LocalModel local_model{ModelName::relxilllp};
  auto spec = default_spec.get_xspec_spectrum();
  local_model.set_par(XPar::h, 3.0);
  REQUIRE_NOTHROW(local_model.eval_model(spec));

  LocalModel local_model_ext{ModelName::relxilllp_ext};
  auto spec_ext = default_spec.get_xspec_spectrum();
  local_model_ext.set_par(XPar::h, 2.0);
  local_model_ext.set_par(XPar::d_ring, 10.0);
  REQUIRE_NOTHROW(local_model_ext.eval_model(spec_ext));


  for (int i = 0; i < spec.num_flux_bins(); i++){
    if (spec.flux[i] - spec_ext.flux[i] > PREC){
      printf("fluxes are different: %e, %e", spec.flux[i], spec_ext.flux[i]);
    }
  }
}


TEST_CASE(" Test setting input parameters for extended source outside the allowed range","[text]") {

  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxilllp_ext);

  auto spec = default_spec.get_xspec_spectrum();

  double height_below_horizon = 0.9;
  lmod.set_par(XPar::h, height_below_horizon);
  lmod.set_par(XPar::d_ring, 0.01);

  REQUIRE_NOTHROW(lmod.eval_model(spec));

  // require a real
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

}

