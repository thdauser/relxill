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

    Copyright 2020 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "catch2/catch.hpp"
#include "../../cppmodels.h"
#include "../../xspec_wrapper_lmodels.h"

#include <vector>

class DefaultSpec {
 public:
  DefaultSpec(double emin, double emax, size_t n_bins) : num_flux_bins{n_bins} {
    size_t n_energy = n_bins + 1;
    energy = new double[n_energy];
    flux = new double[n_bins];

    set_log_grid(energy, n_energy, emin, emax);
    set_input_flux();
  };

  DefaultSpec() : DefaultSpec(0.1, 1000.0, 3000) {
  };

  /* get a logarithmic grid from emin to emax with n_ener bins  */
  static void set_log_grid(double *ener, size_t n_ener, double emin, double emax) {
    for (int ii = 0; ii < n_ener - 1; ii++) {
      ener[ii] = 1.0 * ii / (static_cast<double>(n_ener) - 1.0) * (log(emax) - log(emin)) + log(emin);
      ener[ii] = exp(ener[ii]);
    }
    ener[n_ener - 1] = emax; // separate case (otherwise it is only approx. emax, due to the log/exp functions)
  }

 public:
  double *energy{nullptr};
  double *flux{nullptr};
  const size_t num_flux_bins;

 private:
  void set_input_flux() const {
    flux[num_flux_bins / 2] = 1.0;
  }
};


const double *get_xspec_default_parameter_array(ModelName model_name) {

  auto const model_parameters = ModelDatabase::instance().get(model_name).input_parameters();

  auto default_param_values = ModelParams();
  auto output_param_array = new double[model_parameters.size()];

  for (int ii = 0; ii < model_parameters.size(); ii++) {
    output_param_array[ii] = default_param_values[model_parameters[ii]];
  }

  return output_param_array;
}

double sum_flux(const double *flux, int nbins) {

  double sum = 0.0;
  for (int ii = 0; ii < nbins; ii++) {
    sum += flux[ii];
  }
  return sum;
}

//double sum_flux(const XspecSpectrum &spec) {
//  return sum_flux(spec.flux(), spec.num_flux_bins());
//}

static void test_xspec_lmod_call(ModelName model_name, DefaultSpec default_spec) {

  try {
    const double *xspec_parameters = get_xspec_default_parameter_array(model_name);

    xspec_C_wrapper_eval_model(model_name,
                               xspec_parameters,
                               default_spec.flux,
                               default_spec.num_flux_bins,
                               default_spec.energy);

    delete (xspec_parameters);

    REQUIRE(sum_flux(default_spec.flux, default_spec.num_flux_bins) > 1e-6);
  } catch (ModelNotFound &e) {
    WARN("Skipping test as model not implemented");
  }

}

//static void test_internal_lmod_call(ModelName model_name, const DefaultSpec &default_spec) {
//  LocalModel testModel{model_name};
//
//  XspecSpectrum spec{default_spec.energy, default_spec.flux, default_spec.num_flux_bins};
//
//  testModel.eval_model(spec);
//  REQUIRE(sum_flux(spec) >= 0.0);
//  REQUIRE(sum_flux(spec) > 1e-6);
//
//}

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
TEST_CASE(" testing if local model is implemented", "[basic]") {

  XspecModelDatabase database{};

  for (const auto &elem: database.all_models()) {
    DYNAMIC_SECTION("  - model: " << elem.second.name()) {
      REQUIRE_NOTHROW(get_xspec_default_parameter_array(elem.first));
    }
  }

}


/*
 * TEST CASE
 */
TEST_CASE(" Execute local models", "[model]") {

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
  test_xspec_lmod_call(ModelName::relxilllp, default_spec);
}

TEST_CASE(" Execute Old relxilllp Model", "[old]") {

  const double *xspec_parameters = get_xspec_default_parameter_array(ModelName::relxilllp);

  const int n_parameter = 12;
  int status = EXIT_SUCCESS;
  DefaultSpec default_spec_old{};
  tdrelxilllp(default_spec_old.flux,
              default_spec_old.num_flux_bins,
              default_spec_old.flux,
              xspec_parameters,
              n_parameter,
              &status);

  DefaultSpec default_spec_new{};
  xspec_C_wrapper_eval_model(ModelName::relxilllp,
                             xspec_parameters,
                             default_spec_new.flux,
                             default_spec_new.num_flux_bins,
                             default_spec_new.energy);

  double old_flux = sum_flux(default_spec_old.flux, default_spec_old.num_flux_bins);
  double new_flux = sum_flux(default_spec_new.flux, default_spec_new.num_flux_bins);

  std::cout << "old flux" << old_flux
            << "  --- new flux "
            << new_flux << std::endl;

  REQUIRE(fabs(old_flux - new_flux) < 1e-8);

  delete (xspec_parameters);

}