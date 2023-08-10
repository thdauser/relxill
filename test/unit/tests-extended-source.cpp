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


TEST_CASE("Standard  Test", "[extended]") {

  int status = EXIT_SUCCESS; // in fact this is a bool variable to check the program status

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum(); // afaiu it gives here a default energy, flux values and number of bins used

  LocalModel local_model{ModelName::relxilllp_ext}; // defines a local model object (afaiu)

  // local_model.eval_model(spec); // evaluates a model with given specifications (spec)..?

//  REQUIRE(   sum_flux(spec.flux, spec.num_flux_bins() )  > PREC); // that's obvious - summary flux should be high enough

  auto table = get_es_table(&status); // getting a table with emissivity values - the most important part to manage here


  relParam* rel_param = local_model.get_rel_params(); // getting parameters for the given model...
  // auto rel_param = local_model.get_rel_params(); // the same as line above
  auto sys_par = get_system_parameters(rel_param, &status); // also some parameters

  // auto emis_profile = new_emisProfile(sys_par->re, sys_par->nr, &status);

  // (rel_param).distance = 2.0;
  rel_param->a = 0.998;
  rel_param->height = 3.0;
  rel_param->d_ring = 20.0; // setting some manual distance
  rel_param->emis1 = 3.0;
  rel_param->emis2 = 3.0;
  rel_param->rbr = 30.0;
  rel_param->lineE = 1.0;
  rel_param->return_rad = 0.0;
  // calc_extended_emis(sys_par->emis, rel_param, table, &status);
  get_emis_extendedSource(sys_par->emis, rel_param, &status);
  // with a given table and parameters it should give an emissivity for the extended source (with dist)

  // write_data_to_file("/userdata/data/nekrasov/git/tables/test-ext-emis.dat", sys_par->re, sys_par->emis->emis, sys_par->nr);
  int n_rad = sys_par->nr;
  double *rad = sys_par->re;
  double *intens = sys_par->emis->emis;
  // the part below loads data to file. The file is accessible only locally, so this part failes push to the server if not commented.
  // Need either to modify this part, or not to use it. The same is in the test below
/*
  string strname = "/home/nekrasov/git/tables/test-ext-emis-d" + std::to_string(rel_param->d_ring) + "_new.dat";
  char *foutName = new char[999];
  strcpy(foutName, strname.c_str());
  FILE *fp = fopen(foutName, "w+");
  if (fp != nullptr) {
    for (int ii = 0; ii < n_rad; ii++) {
      fprintf(fp, " %e \t %e \n", rad[ii], intens[ii]);
    }
  }
  fclose(fp);
*/
}


TEST_CASE("Default line model call  Test", "[extended]") {

 int status = EXIT_SUCCESS; // in fact this is a bool variable to check the program status

DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum(); // afaiu it gives here a default energy, flux values and number of bins used

  LocalModel local_model{ModelName::relline_lp_ext}; // defines a local model object (afaiu)

  double height = 2.0;
  double dist = 20.0;
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::d_ring, dist);

  local_model.eval_model(spec); // evaluates a model with given specifications (spec)..?

  int n_bins = spec.num_flux_bins();
  double *flux = spec.flux;
  double *energy = spec.energy;
/*
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



TEST_CASE("Xspec line model call  Test", "[extended]") {

  int status = EXIT_SUCCESS; // in fact this is a bool variable to check the program status

  // const double * params = []
  // xspec_C_wrapper_eval_model(ModelName::relline_lp_ext, parameter, flux, Nflux, energy);
  LocalModel local_model{ModelName::relline_lp_ext}; //parameter_values,
  local_model.set_par(XPar::d_ring, 5.0);
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  // XspecSpectrum spectrum{xspec_energy, xspec_flux, static_cast<size_t>(num_flux_bins)};

  local_model.eval_model(spec);
  // printf("%f, %f, %f", spec.flux[1000], spec.flux[1200], spec.flux[1600]);

}


TEST_CASE("Default relxill model call  Test", "[extended]") {

  int status = EXIT_SUCCESS; // in fact this is a bool variable to check the program status

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum(); // afaiu it gives here a default energy, flux values and number of bins used

  LocalModel local_model{ModelName::relxilllp_ext}; // defines a local model object (afaiu)

  double height = 2.0;
  double dist = 20.0;
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::d_ring, dist);

  local_model.eval_model(spec); // evaluates a model with given specifications (spec)..?

  int n_bins = spec.num_flux_bins();
  double *flux = spec.flux;
  double *energy = spec.energy;
  // printf(" %e \t %e \n", flux[100], energy[100]);
  // printf(" %e \t %e \n", flux[300], energy[300]);
  // printf(" %e \t %e \n", flux[600], energy[600]);
  /*
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

  int status = EXIT_SUCCESS; // in fact this is a bool variable to check the program status

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum(); // afaiu it gives here a default energy, flux values and number of bins used
  LocalModel local_model{ModelName::relxilllp}; // defines a local model object (afaiu)
  relParam* rel_param = local_model.get_rel_params(); // getting parameters for the given model...
  auto sys_par = get_system_parameters(rel_param, &status); // also some parameters

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
  // with a given table and parameters it should give an emissivity for the extended source (with dist)
  local_model.eval_model(spec);

  LocalModel local_model_ext{ModelName::relxilllp_ext}; // defines a local model object (afaiu)
  relParam* rel_param_ext = local_model_ext.get_rel_params(); // getting parameters for the given model...
  auto sys_par_ext = get_system_parameters(rel_param_ext, &status); // also some parameters

  // (rel_param).distance = 2.0;
  rel_param_ext->a = 0.998;
  rel_param_ext->height = 3.0;
  rel_param_ext->d_ring = 0.1; // setting some manual distance
  rel_param_ext->emis1 = 3.0;
  rel_param_ext->emis2 = 3.0;
  rel_param_ext->rbr = 30.0;
  rel_param_ext->lineE = 1.0;
  rel_param_ext->return_rad = 0.0;
  get_emis_extendedSource(sys_par_ext->emis, rel_param_ext, &status);
  local_model_ext.eval_model(spec);
  int n_rad = sys_par->nr;
  double *rad = sys_par->re;
  double *intens = sys_par->emis->emis;
  int n_rad_ext = sys_par_ext->nr;
  double *rad_ext = sys_par_ext->re;
  double *intens_ext = sys_par_ext->emis->emis;
/*
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

