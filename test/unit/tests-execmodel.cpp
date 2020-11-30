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


TEST_CASE(" Execute local models", "[model]") {

  struct TestSpec {
    const Array energy{0.1, 1.0, 10.0};
    Array flux{0.0, 0.0, 0.0};
  } testSpec;

  const std::unordered_map<ModelName, std::string> all_models{
      //    {ModelName::relxill, "relxill"},
      //    {ModelName::relline, "relline"},
      {ModelName::relconv, "relconv"}
      //    {ModelName::xillver, "xillver"}
  };

  for (const auto &elem: all_models) {

    DYNAMIC_SECTION(" testing model: " << elem.second) {
      LocalModel testModel{ModelDatabase.instance().get(elem.first)};
      std::cout << "- model: " << elem.second << std::endl;
      testModel.eval_model(testSpec.energy, testSpec.flux);
      REQUIRE(testSpec.flux[0] >= 0.0);
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
