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
#include "JedSad.h"


TEST_CASE(" Execute relxill jedsad", "[jedsad]") {
  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxill_jedsad);

  auto spec = default_spec.get_xspec_spectrum();
  lmod.eval_model(spec);
  REQUIRE_NOTHROW(lmod.eval_model(spec));

  double sum = sum_flux(spec.flux, spec.num_flux_bins());
  REQUIRE(sum > 1E-8);
// not yet working, model not implemented  REQUIRE(sum > 1e-8);

}

TEST_CASE(" Read table", "[jedsad]"){

  auto tab = &JedsadTable::instance();

 // if (tab.data.empty())
 //   tab.read_table();

  // check the table dimensions
  int ncols = tab->raw_data.size();
  REQUIRE(ncols == 2);

  int nvals = tab->num_data_vals();

  auto ref_val = 0.904671;
  REQUIRE(fabs(ref_val - tab->raw_data[0][9599]) < 1e-6);

  auto avg_gam = calcSum(tab->raw_data[0], nvals) / nvals;
  auto avg_ect = calcSum(tab->raw_data[1], nvals) / nvals;


  // printf("Avg Gam %e \n", avg_gam );
  // printf("Avg Ect %e \n", avg_ect );

  REQUIRE(avg_gam > 0);
  REQUIRE(avg_gam < 3.4);
  REQUIRE(avg_ect > 0);

}

TEST_CASE(" Parameter Values in the table", "[jedsad]") {

  auto tab = &JedsadTable::instance();

// values of the first and second entry for each parameter
  std::vector<double> vec = {10, 2, 0.01, 2, 0.01, 0.5, 0.5, 0.3, 0.01};
  std::vector<double> vec2 = {10, 2.464, 0.01, 2, 0.01, 0.5, 0.5, 0.3, 0.01};


  auto _p_vals = tab->param_vals;
  REQUIRE(_p_vals.size() == vec.size());

  for (int i = 0; i < _p_vals.size(); i++) {
    REQUIRE(fabs(_p_vals[i][0] - vec[i]) < 1e-6);
  }
  REQUIRE(fabs(_p_vals[1][1] - vec2[1]) < 1e-6);

  auto res1 = tab->interpolate(vec.data());
  vec2[1] = 5.68;
  vec2[6] = 0.5;
  auto res2 = tab->interpolate(vec2.data());
  vec2[6] = 0.55;
  auto res3 = tab->interpolate(vec2.data());
  vec2[6] = 0.6;
  auto res4 = tab->interpolate(vec2.data());
  //printf("res2 %.3e %.3e \n", res2[0], res2[1] );
  //printf("res3 %.3e %.3e \n", res3[0], res3[1] );
  //printf("res4 %.3e %.3e \n", res4[0], res4[1] );

  REQUIRE(res3[1] - fabs(0.5 * res2[1] + 0.5 * res4[1]) < 1e-6);

  std::vector<int> pind = {0, 5, 0, 0, 0, 0, 0, 0, 0};
  auto flat_pind = tab->get_data_index(pind);

  std::vector<int> pind2 = {0, 5, 0, 0, 0, 0, 1, 0, 0};
  auto flat_pind2 = tab->get_data_index(pind2);

  REQUIRE(flat_pind == flat_pind2 - 1);  // as the parameter 8 and 9 only have 1 entry
  // printf( "flat_pind %d \n", flat_pind );

  // based on the table  "Refl_XrB_DauserTest_lb3.fits"
  std::vector<double> fref1 = {0.0, 0.0};  // flat ind 0
  std::vector<double> fref2 = {2.063146, 2893.463};  // flat ind "flat_pind"

  for (int i = 0; i < 2; i++) {
    // check reference values
    REQUIRE(fabs(fref1[i] - res1[i]) < 1e-6);
    // check linear interpolation (in one parameter)
    REQUIRE(res3[i] - fabs(0.5 * res2[i] + 0.5 * res4[i]) < 1e-6);
  }

  // printf("res2 %.3e %.3e \n", res2[0], res2[1] );

}

