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
  DefaultSpec(double emin, double emax, size_t n_bins) {
    energy = Array(n_bins);
    flux = Array(n_bins);

    set_log_grid(energy, emin, emax);

    //   spec = new CppSpectrum(energy, flux);
  };

  DefaultSpec() : DefaultSpec(0.1, 1000.0, 3000) {
  };


  //  ~DefaultSpec(){
  //    delete spec;
  //  }

  /* get a logarithmic grid from emin to emax with n_ener bins  */
  static void set_log_grid(Array &ener, double emin, double emax) {
    size_t n_ener = ener.size();
    for (int ii = 0; ii < n_ener - 1; ii++) {
      ener[ii] = 1.0 * ii / (static_cast<double>(n_ener) - 1.0) * (log(emax) - log(emin)) + log(emin);
      ener[ii] = exp(ener[ii]);
    }
    ener[n_ener - 1] = emax; // separate case (otherwise it is only approx. emax, due to the log/exp functions)
  }

 public:
  //  CppSpectrum* spec{nullptr};
  Array energy;
  Array flux;

};

Array get_xspec_default_parameter_array(ModelName model_name) {

  auto const model_parameters = ModelDatabase::instance().get(model_name).input_parameters();

  auto default_param_values = ModelParams();
  //  auto output_param_array = new Array[model_parameters.size()];
  auto output_param_array = Array(model_parameters.size());

  for (int ii = 0; ii < output_param_array.size(); ii++) {
    output_param_array[ii] = default_param_values[model_parameters[ii]];
  }

  return output_param_array;
}


//double sum_flux(const CppSpectrum &spec){
//
//
//  return std::accumulate(std::begin(spec.flux()), std::end(spec.flux()), 0.0);
//
//}


static void test_xspec_lmod_call(ModelName model_name, DefaultSpec default_spec) {
  Array parameters = get_xspec_default_parameter_array(model_name);
  xspec_wrapper_eval_model(model_name, default_spec.energy, default_spec.flux, parameters);

  REQUIRE(default_spec.flux.max() >= 0.0);
}

static void test_internal_lmod_call(ModelName model_name, const DefaultSpec &default_spec) {
  LocalModel testModel{model_name};

  const auto energy = default_spec.energy;
  auto input_flux = default_spec.flux;
  CppSpectrum spec{energy, input_flux};

  REQUIRE(spec.flux().max() >= 0.0);
  testModel.eval_model(spec);
  auto flux = spec.flux();
  REQUIRE(flux.max() >= 0.0);

}

TEST_CASE(" default spectrum class", "[basic]") {

  DefaultSpec default_spec{};

  REQUIRE(default_spec.energy[0]);
  REQUIRE(default_spec.energy[1] > default_spec.energy[0]);
  // REQUIRE(default_spec.spec);


  double emin = 0.5;
  double emax = 10.0;
  size_t nbins = 100;
  auto own_spec = DefaultSpec(emin, emax, nbins);

  REQUIRE(own_spec.energy[0] == emin);
  REQUIRE(own_spec.energy[1] > emin);
  REQUIRE(own_spec.energy.size() == nbins);
  REQUIRE(own_spec.flux.size() == nbins);
  REQUIRE(own_spec.energy[own_spec.energy.size() - 1] == emax);

}

TEST_CASE(" Execute local models", "[model]") {

  const std::unordered_map<ModelName, std::string> all_models{
      {ModelName::relline, "relline"}
      //      {ModelName::relxill, "relxill"},
      //      {ModelName::relconv, "relconv"},
      //      {ModelName::xillver, "xillver"}
  };

  DefaultSpec default_spec{};

  for (const auto &elem: all_models) {

    auto model_name_type = elem.first;
    auto model_name_string = elem.second;

    DYNAMIC_SECTION(" testing model: " << model_name_string) {
      std::cout << "- test model: " << elem.second << std::endl;

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



//LmodTest eval_local_model(ModelName name){
//  LmodTest inp{};
//
//  auto const model_definition = ModelDatabase::instance().get(ModelName::relxill);
//
//  LocalModel model{model_definition.model_info()};
//  model.eval_model(inp.energy, inp.flux);
//
//  return inp;
//}
//
//TEST_CASE(" Execute relxill ") {
//  LmodTest lmod = eval_local_model(ModelName::relxill);
//  REQUIRE(lmod.flux[0] >= 0.0);
//}

