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

#include <vector>
#include <unordered_map>
#include <iostream>

class SimpleSpec {

 public:
  double energy[4] = {0.1, 1.0, 3.0, 10.0};
  double flux[3] = {1.0, 1.0, 10.0};
  int n_flux_bins = 3;
};

/*
 * Testing of cppspectrum.h
 */
TEST_CASE(" Spectrum Class", "[basic]") {

  SimpleSpec test_spec{};
  XspecSpectrum spec(test_spec.energy, test_spec.flux, test_spec.n_flux_bins);

  DYNAMIC_SECTION(" test initial array without operations ") {
    REQUIRE(spec.energy());
    REQUIRE(spec.flux());

    REQUIRE(spec.energy());
    REQUIRE(spec.flux());
  }

  DYNAMIC_SECTION(" is the given number of flux bins and energy bins correctly set") {
    REQUIRE(spec.num_flux_bins() == test_spec.n_flux_bins);
    REQUIRE(spec.n_energy() == test_spec.n_flux_bins + 1);
  }

}

TEST_CASE(" Spectrum Class:  shifting of the energy grid ") {

  SimpleSpec test_spec{};
  XspecSpectrum spec(test_spec.energy, test_spec.flux, test_spec.n_flux_bins);

  double ener0 = spec.energy()[0];

  DYNAMIC_SECTION(" shifting by line energy produces correct energies? ") {
    spec.shift_energy_grid_1keV(0.5);
    double shifted_ener = spec.energy()[0];
    REQUIRE(shifted_ener == ener0 * 2);
  }

  DYNAMIC_SECTION(" shifting by z correct energies? ") {
    spec.shift_energy_grid_1keV(1.0);
    spec.shift_energy_grid_redshift(0.5);
    double shifted_ener = spec.energy()[0];
    REQUIRE(shifted_ener == ener0 * 1.5);
  }

}
