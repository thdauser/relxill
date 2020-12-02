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

#ifndef RELXILL_TEST_UNIT_TESTS_EXECMODEL_H_
#define RELXILL_TEST_UNIT_TESTS_EXECMODEL_H_

#include "catch2/catch.hpp"
#include "../../cppmodels.h"

#include <vector>
#include <unordered_map>
#include <iostream>



//TEST_CASE(" Execute local models", "[model]") {
//
//
//  LmodTest inp{};
//
//  for (const auto &elem: all_models) {
//
//    DYNAMIC_SECTION(" testing model: " << elem.second) {
//      std::cout << "- model: " << elem.second << std::endl;
//      xspec_wrapper_eval_model(elem.first, inp.energy, inp.flux, inp.parameter);
//      REQUIRE(inp.flux[0] >= 0.0);
//    }
//
//  }
//
//
//}

class TestSpec {
 public:
  TestSpec() :
      test_spec{CppSpectrum(energy, flux)} {
  };

 public:
  const Array energy{0.1, 1.0, 10.0};
  Array flux{0.0, 0.0, 0.0};
  CppSpectrum test_spec;
};

TEST_CASE(" Spectrum Class") {

  TestSpec test_spec{};
  CppSpectrum spec(test_spec.energy, test_spec.flux);

  DYNAMIC_SECTION(" test initial array without operations ") {
    REQUIRE(typeid(spec.energy()) == typeid(Array));
    REQUIRE(typeid(spec.flux()) == typeid(Array));

    REQUIRE(typeid(spec.energy_double()[0]) == typeid(double));
    REQUIRE(typeid(spec.flux_double()[0]) == typeid(double));
  }

}

TEST_CASE(" Execute local models", "[model]") {

  const std::unordered_map<ModelName, std::string> all_models{
      {ModelName::relxill, "relxill"},
      {ModelName::relline, "relline"},
      {ModelName::relconv, "relconv"},
      {ModelName::xillver, "xillver"}
  };

  TestSpec test_spec{};

  CppSpectrum spec(test_spec.energy, test_spec.flux);

  for (const auto &elem: all_models) {

    auto model_name_type = elem.first;
    auto model_name_string = elem.second;

    DYNAMIC_SECTION(" testing model: " << model_name_string) {
      LocalModel testModel{model_name_type};
      std::cout << "- test model: " << elem.second << std::endl;
      testModel.eval_model(spec);
      REQUIRE(spec.flux()[0] >= 0.0);
    }

  }

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

#endif //RELXILL_TEST_UNIT_TESTS_EXECMODEL_H_
