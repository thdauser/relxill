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

extern "C" {
#include "writeOutfiles.h"
}

#define PREC 1e-6


TEST_CASE("Standard  Test", "[extended]") {

  int status = EXIT_SUCCESS; // in fact this is a bool variable to check the program status

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum(); // afaiu it gives here a default energy, flux values and number of bins used

  LocalModel local_model{ModelName::relxilllpCp}; // defines a local model object (afaiu)

  // local_model.eval_model(spec); // evaluates a model with given specifications (spec)..?

//  REQUIRE(   sum_flux(spec.flux, spec.num_flux_bins() )  > PREC); // that's obvious - summary flux should be high enough

  float *val;
  fitsfile *fptr = nullptr;
  auto table = get_es_table(&status); // getting a table with emissivity values - the most important part to manage here



  relParam* rel_param = local_model.get_rel_params(); // getting parameters for the given model...
  // auto rel_param = local_model.get_rel_params(); // the same as line above
  auto sys_par = get_system_parameters(rel_param, &status); // also some parameters

  // auto emis_profile = new_emisProfile(sys_par->re, sys_par->nr, &status);



  // (rel_param).distance = 2.0;
  rel_param->a = 0.99;
  rel_param->height = 5.0;
  rel_param->d_ring = 10.0; // setting some manual distance
  // calc_extended_emis(sys_par->emis, rel_param, table, &status);
  get_emis_extendedSource(sys_par->emis, rel_param, &status);
  // with a given table and parameters it should give an emissivity for the extended source (with dist)

  // write_data_to_file("/userdata/data/nekrasov/git/tables/test-ext-emis.dat", sys_par->re, sys_par->emis->emis, sys_par->nr);
  int n_rad = sys_par->nr;
  double *rad = sys_par->re;
  double *intens = sys_par->emis->emis;
  /*
  string strname = "test-ext-emis-d" + std::to_string(rel_param->d_ring) + "_v5.dat";
  char *foutName = new char[999];
  strcpy(foutName, strname.c_str());
  FILE *fp = fopen(foutName, "w+");
  for (int ii = 0; ii < n_rad; ii++) {
    fprintf(fp, " %e \t %e \n", rad[ii], intens[ii]);
  }
  fclose(fp);
*/
}


TEST_CASE("Default model call  Test", "[extended]") {

 int status = EXIT_SUCCESS; // in fact this is a bool variable to check the program status

DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum(); // afaiu it gives here a default energy, flux values and number of bins used

  LocalModel local_model{ModelName::relline_lp_ext}; // defines a local model object (afaiu)

  double height = 20.0;
  double dist = 3.0;
  local_model.set_par(XPar::h, height);
  local_model.set_par(XPar::d_ring, dist);

  local_model.eval_model(spec); // evaluates a model with given specifications (spec)..?

  int n_bins = spec.num_flux_bins();
  double *flux = spec.flux;
  double *energy = spec.energy;
  /*
  string strname = "test-ext-line-h" + std::to_string(height) + "-d" + std::to_string(dist) + ".dat";
  char *foutName = new char[999];
  strcpy(foutName, strname.c_str());
  FILE *fp = fopen(foutName, "w+");
  for (int ii = 0; ii < n_bins; ii++) {
    fprintf(fp, " %e \t %e \n", energy[ii], flux[ii]);
  }
  fclose(fp);
*/
}
