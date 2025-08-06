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

    Copyright 2025 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "subs_speed_test.h"

// ------------------------- //
int main(int argc, char *argv[]) {

  const int num_evaluations = 100;

  const int num_zones = 25;
  setenv("RELXILL_NUM_RZONES", std::to_string(num_zones).c_str(), 1);
  printf(" testing relxill_jedsad for number of zones %i \n", num_zones);

  auto tstart = std::chrono::steady_clock::now();

  static_assert(num_evaluations > 1);
  eval_local_model_param_range(ModelName::relxill_jedsad, XPar::rj, 5.0, 10.0, num_evaluations);


  auto time_elapsed_msec = std::chrono::duration_cast<std::chrono::milliseconds>
      (std::chrono::steady_clock::now() - tstart).count();

  printf("Time elapsed:              %.2fsec\n", static_cast<double>(time_elapsed_msec) * 0.001);
  printf("Time per model evaluation: %.0fmsec ", static_cast<double>(time_elapsed_msec) / num_evaluations);


  return EXIT_SUCCESS;
}
