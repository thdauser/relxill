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

#include <vector>
#include <unordered_map>
#include <iostream>

class DefaultSpec {
 public:
  DefaultSpec(double emin, double emax, size_t n_bins) : num_flux_bins{n_bins} {
    size_t n_energy = n_bins + 1;
    energy = new double[n_energy];
    flux = new double[n_bins];

    set_log_grid(energy, n_energy, emin, emax);
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
};

const double *get_xspec_default_parameter_array(ModelName model_name) {

  auto const model_parameters = ModelDatabase::instance().get(model_name).input_parameters();

  auto default_param_values = ModelParams();
  //  auto output_param_array = new Array[model_parameters.size()];
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
double sum_flux(const XspecSpectrum &spec) {
  return sum_flux(spec.flux(), spec.num_flux_bins());
}

static void test_xspec_lmod_call(ModelName model_name, DefaultSpec default_spec) {
  const double *xspec_parameters = get_xspec_default_parameter_array(model_name);
  xspec_C_wrapper_eval_model(model_name,
                             xspec_parameters,
                             default_spec.flux,
                             default_spec.num_flux_bins,
                             default_spec.energy);

  REQUIRE(sum_flux(default_spec.flux, default_spec.num_flux_bins) >= 0.0);

  delete (xspec_parameters);
}

static void test_internal_lmod_call(ModelName model_name, const DefaultSpec &default_spec) {
  LocalModel testModel{model_name};

  XspecSpectrum spec{default_spec.energy, default_spec.flux, default_spec.num_flux_bins};

  testModel.eval_model(spec);
  REQUIRE(sum_flux(spec) >= 0.0);
  REQUIRE(sum_flux(spec) > 1e-6);

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
TEST_CASE(" Execute local models", "[model]") {

  const std::unordered_map<ModelName, std::string> all_models{
      {ModelName::relline, "relline"},
      {ModelName::relxill, "relxill"},
      {ModelName::relconv, "relconv"},
      {ModelName::xillver, "xillver"}
  };

  DefaultSpec default_spec{};

  for (const auto &elem: all_models) {

    auto model_name_type = elem.first;
    auto model_name_string = elem.second;

    DYNAMIC_SECTION(" testing model: " << model_name_string) {
      // std::cout << "- test model: " << elem.second << std::endl;

      DYNAMIC_SECTION(" xspec local model call ") {
        test_xspec_lmod_call(model_name_type, default_spec);
      }

      DYNAMIC_SECTION(" internal model call ") {
        test_internal_lmod_call(model_name_type, default_spec);
      }

    }

  }

  free_cache(); // TODO: find a better place to free the cache

}
