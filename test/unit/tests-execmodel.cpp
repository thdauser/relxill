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


class TestSpec {
 public:
  TestSpec() = default;
 public:
  const Array energy{0.1, 1.0, 3.0, 10.0};
  Array flux{0.0, 1.0, 0.0};
};


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

