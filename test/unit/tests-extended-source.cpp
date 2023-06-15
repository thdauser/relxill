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

  local_model.eval_model(spec); // evaluates a model with given specifications (spec)..?

  REQUIRE(   sum_flux(spec.flux, spec.num_flux_bins() )  > PREC); // that's obvious - summary flux should be high enough

  float *val;
  fitsfile *fptr = nullptr;
  // fits_open_table(&fptr, ESTABLE_FILENAME, READONLY, &status);
  // fits_movabs_hdu(fptr, 2, nullptr, &status);
  // fits_read_key(fptr, TFLOAT, "SPIN", val, nullptr, &status);
  // printf("%.2f", val);

  auto table = get_es_table(&status); // getting a table with emissivity values - the most important part to manage here



  relParam* rel_param = local_model.get_rel_params(); // getting parameters for the given model...
  // auto rel_param = local_model.get_rel_params(); // the same as line above
  auto sys_par = get_system_parameters(rel_param, &status); // also some parameters

  // auto emis_profile = new_emisProfile(sys_par->re, sys_par->nr, &status);



  // (rel_param).distance = 2.0;
  rel_param->a = 0.99;
  rel_param->height = 5.0;
  rel_param->distance = 25.0; // setting some manual distance
  // calc_extended_emis(sys_par->emis, rel_param, table, &status);
  get_emis_extendedSource(sys_par->emis, rel_param, &status);
  // with a given table and parameters it should give an emissivity for the extended source (with dist)

  // write_data_to_file("/userdata/data/nekrasov/git/tables/test-ext-emis.dat", sys_par->re, sys_par->emis->emis, sys_par->nr);
  int n_rad = sys_par->nr;
  double *rad = sys_par->re;
  double *intens = sys_par->emis->emis;
  string strname = "/userdata/data/nekrasov/git/tables/test-ext-emis-d" + std::to_string(rel_param->distance) + "_v4.dat"; // std::round(rel_param->distance / 0.1) * 0.1
  char *foutName = new char[999];
  strcpy(foutName, strname.c_str());
  // std::to_string(rel_param->a) + std::to_string(rel_param->height) +
  FILE *fp = fopen(foutName, "w+");
  // fprintf(fp, " %s \t %s \n", "rad", "emis");
  printf(" %i", n_rad);
  for (int ii = 0; ii < n_rad; ii++) {
    fprintf(fp, " %e \t %e \n", rad[ii], intens[ii]);
  }
  fclose(fp);

}


///
// TEST_CASE("Standard  Test", "[extended]") {
//  param->distance = 5.0;
//  calc_extended_emis(emisProf, param, tab, status);
// }