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

    Copyright 2024 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "XspecSpectrum.h"

#include <chrono>

TEST_CASE(" caching of energy grid", "[cache]") {

  DefaultSpec default_spec{};

  const double elo = 0.1;
  const double eup1 = 100.0;
  const double eup2 = 150.0;
  const int nbins = 1000;

  DefaultSpec def_spec1{elo, eup1, nbins};
  auto spec1 = def_spec1.get_xspec_spectrum();

  DefaultSpec def_spec2{elo, eup2, nbins};
  auto spec2 = def_spec2.get_xspec_spectrum();

  LocalModel lmod(ModelName::relxilllpCp);

  // first evaluation
  lmod.eval_model(spec1);

  // not cached
  auto tstart = std::chrono::steady_clock::now();
  lmod.set_par(XPar::a, 0.942);
  lmod.eval_model(spec1);
  auto time_elapsed_msec = std::chrono::duration_cast<std::chrono::milliseconds>
      (std::chrono::steady_clock::now() - tstart).count();

  printf("Time elapsed for model evaluation:      %.3fmsec\n", static_cast<double>(time_elapsed_msec));

  // should be cached (100 evaluations)
  const int n_eval_cache = 100;
  tstart = std::chrono::steady_clock::now();
  for (int ii = 0; ii < n_eval_cache; ii++) {
    lmod.eval_model(spec2);
  }
  auto time_elapsed_msec_cache = std::chrono::duration_cast<std::chrono::milliseconds>
      (std::chrono::steady_clock::now() - tstart).count();
  printf("Time elapsed for energy grid change (per eval):    %.3fmsec\n",
         static_cast<double>(time_elapsed_msec_cache) / n_eval_cache);

  // require that model evaluation takes longer than n_eval_cache*cached
  REQUIRE(time_elapsed_msec > time_elapsed_msec_cache);

  double flu1 = spec1.get_energy_flux();
  double flu2 = spec2.get_energy_flux();

  REQUIRE(fabs(flu1 / flu2 - 1) > 1e-6);

}

