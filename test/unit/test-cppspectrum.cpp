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
  Array flux{1.0, 1.0, 10.0, 0.0};
};

/*
 * Testing of cppspectrum.h
 */
TEST_CASE(" Spectrum Class", "[basic]") {

  TestSpec test_spec{};
  CppSpectrum spec(test_spec.energy, test_spec.flux);

  DYNAMIC_SECTION(" test initial array without operations ") {
    REQUIRE(typeid(spec.energy()) == typeid(Array));
    REQUIRE(typeid(spec.flux()) == typeid(Array));

    REQUIRE(spec.energy_double());
    REQUIRE(spec.flux_double());
  }

  DYNAMIC_SECTION(" is the given number of energy bins correct (equal flux bins)") {
    REQUIRE(spec.nener_bins() == spec.flux().size() - 1);
    REQUIRE(spec.nener_bins() == spec.energy().size() - 1);
  }

  DYNAMIC_SECTION(" are Array and Double-Arrays the same") {
    for (int ii = 0; ii < spec.nener_bins(); ii++) {
      REQUIRE(spec.energy()[ii] == spec.energy_double()[ii]);
      REQUIRE(spec.flux()[ii] == spec.flux_double()[ii]);
    }

  }

  DYNAMIC_SECTION(" is the energy_double() returning double type") {
    REQUIRE(typeid(spec.energy_double()[0]) != typeid(int));
    REQUIRE(typeid(spec.flux_double()[0]) != typeid(int));

    REQUIRE(typeid(spec.energy_double()[0]) == typeid(double));
    REQUIRE(typeid(spec.flux_double()[0]) == typeid(double));
  }

}

TEST_CASE(" Spectrum Class:  shifting of the energy grid ") {

  TestSpec test_spec{};
  CppSpectrum spec(test_spec.energy, test_spec.flux);

  double ener0 = spec.energy()[0];

  DYNAMIC_SECTION(" shifting produces valid energy arrays ") {
    spec.shift_energy_grid_1keV(0.5, 0.0);
    REQUIRE(spec.energy_double());
    REQUIRE(spec.flux_double());
  }

  DYNAMIC_SECTION(" shifting by line energy produces correct energies? ") {
    spec.shift_energy_grid_1keV(0.5, 0.0);
    double shifted_ener = spec.energy()[0];
    REQUIRE(shifted_ener == ener0 * 2);
    double shifted_ener_double = spec.energy_double()[0];
    REQUIRE(shifted_ener == shifted_ener_double);
  }

  DYNAMIC_SECTION(" shifting by z correct energies? ") {
    spec.shift_energy_grid_1keV(1.0, 0.5);
    double shifted_ener = spec.energy()[0];
    REQUIRE(shifted_ener == ener0 * 1.5);
  }

}
