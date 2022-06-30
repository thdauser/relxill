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



static void test_xspec_lmod_call(ModelName model_name, const DefaultSpec& default_spec) {

  try {
    eval_xspec_lmod_default(model_name, default_spec);
     REQUIRE(sum_flux(default_spec.flux, default_spec.num_flux_bins) > 1e-6);
  } catch (ModelNotFound &e) {
    WARN("Skipping test as model not implemented");
  }

}

/*
 * TEST CASE
 */
TEST_CASE(" default spectrum class", "[basic]") {

  DefaultSpec default_spec{};

  REQUIRE(default_spec.energy[0]);
  REQUIRE(default_spec.energy[1] > default_spec.energy[0]);

  double emin = 0.5;
  double emax = 10.0;
  size_t nbins = 100;
  auto own_spec = DefaultSpec(emin, emax, nbins);

  REQUIRE(own_spec.energy[0] == emin);
  REQUIRE(own_spec.energy[1] > emin);
  REQUIRE(own_spec.num_flux_bins == nbins);
  REQUIRE(own_spec.energy[own_spec.num_flux_bins] == emax);  // energy has num+1 bins

}

/*
 * TEST CASE
 */
static void test_loading_default_parameters(ModelName model_name){
  const double* param = get_xspec_default_parameter_array(model_name);
  delete[] param;
}

TEST_CASE(" testing if local model is implemented", "[basic]") {

  XspecModelDatabase database{};

  for (const auto &elem: database.all_models()) {
    DYNAMIC_SECTION("  - model: " << elem.second.name()) {
      REQUIRE_NOTHROW(test_loading_default_parameters(elem.first));
    }
  }

}


/*
 * TEST CASE
 */
TEST_CASE(" Execute ALL local models", "[model]") {

  DefaultSpec default_spec{};
  XspecModelDatabase database{};

  for (const auto &elem: database.all_models()) {

    auto model_name_type = elem.first;
    auto model_name_string = elem.second.name();

    DYNAMIC_SECTION(" - model: " << model_name_string) {
      test_xspec_lmod_call(model_name_type, default_spec);
    }

  }

  free_cache(); // TODO: find a better place to free the cache

}

TEST_CASE(" Execute single model", "[single]") {
  DefaultSpec default_spec{};
  REQUIRE_NOTHROW(test_xspec_lmod_call(ModelName::relxilllp, default_spec) );

  LocalModel lmod(ModelName::relxilllp);

  auto spec = default_spec.get_xspec_spectrum();
  REQUIRE_NOTHROW(lmod.eval_model(spec));

  lmod.set_par(XPar::h,6.00);
  lmod.set_par(XPar::rin,-1.088889);
  REQUIRE_NOTHROW(lmod.eval_model(spec));

}

TEST_CASE(" Exec single model with LocalModel Structure","[model-change]") {
  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxilllpCp);

  auto spec = default_spec.get_xspec_spectrum();
  lmod.set_par(XPar::iongrad_index, 1.0);
  lmod.set_par(XPar::logn, 19.5);

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  double sum = sum_flux(spec.flux, spec.num_flux_bins());

  REQUIRE( sum > 1e-8);

}



static void require_file_exists(const string& fname){
  std::filesystem::path f{ fname };
  INFO("trying to find file: " +  fname );
  REQUIRE(std::filesystem::exists(f));
}

TEST_CASE(" Execute single model with output writing ", "[output]") {
  DefaultSpec default_spec{};
  const char* env_outfiles = "RELXILL_WRITE_FILES";

  setenv(env_outfiles, "1", 1);

  eval_xspec_lmod_default(ModelName::relxilllp, default_spec);
  unsetenv(env_outfiles);

  require_file_exists("__relxillOutput_rellineProfile.dat");
  require_file_exists("__relxillOutput_emisProfile.dat");
}


TEST_CASE(" Execute iongrad model with output writing ", "[output]") {
  DefaultSpec default_spec{};
  const char* env_outfiles = "RELXILL_WRITE_FILES";

  setenv(env_outfiles, "1", 1);

  LocalModel lmod(ModelName::relxilllpCp);
  auto spec = default_spec.get_xspec_spectrum();

  lmod.set_par(XPar::switch_iongrad_type, 2);
  lmod.eval_model(spec);

  unsetenv(env_outfiles);

  require_file_exists("__relxillOutput_iongrad.dat");
}



TEST_CASE(" Test setting input parameters outside the allowed range","[test]") {
  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxilllp);

  auto spec = default_spec.get_xspec_spectrum();

  double height_below_horizon = 0.9;
  lmod.set_par(XPar::h, height_below_horizon);

  REQUIRE_NOTHROW(lmod.eval_model(spec));

  // require a real
  REQUIRE(sum_flux(spec.flux, spec.num_flux_bins()) > 1e-8);

}




TEST_CASE(" Relconv applied outside defined energy range should reset to zero", "[new]") {

  const double test_val = 1000.0;

  auto input_spec_outside = DefaultSpec(0.001, 0.01, 1);
  input_spec_outside.flux[input_spec_outside.num_flux_bins/2] = test_val;

  auto spec_outside = input_spec_outside.get_xspec_spectrum();

  LocalModel lmod(ModelName::relconv);
  REQUIRE_NOTHROW(lmod.eval_model(spec_outside));

  // should be set to zero
  REQUIRE(sum_flux(spec_outside.flux, spec_outside.num_flux_bins()) < 1e-8);

  auto input_spec = DefaultSpec(0.01, 1, 100);
  input_spec.flux[input_spec.num_flux_bins/2] = test_val;

  auto spec = input_spec.get_xspec_spectrum();
  REQUIRE_NOTHROW(lmod.eval_model(spec));

  // should be very close to the test value  (as it is normalized)
  REQUIRE( fabs(sum_flux(spec.flux, spec.num_flux_bins()) - test_val) < 0.1);


}



