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

class LmodTest {

 public:
  const Array energy{0.1, 1.0, 10.0};
  Array flux{0.0, 0.0, 0.0};
  Array parameter{0.0, 0.0, 0.0};
  const std::unordered_map<ModelName, std::string> all_models{
      {ModelName::relline, "relline"},
      {ModelName::relxill, "relxill"},
      {ModelName::relconv, "relconv"},
      {ModelName::xillver, "xillver"}
  };

};

TEST_CASE(" Execute local models", "[model]") {

  LmodTest inp{};

  for (const auto &elem: inp.all_models) {

    DYNAMIC_SECTION(" testing model: " << elem.second) {
      std::cout << "- model: " << elem.second << std::endl;
      eval_model_xspec(elem.first, inp.energy, inp.flux, inp.parameter);
      REQUIRE(inp.flux[0] >= 0.0);
    }

  }
}

#endif //RELXILL_TEST_UNIT_TESTS_EXECMODEL_H_
