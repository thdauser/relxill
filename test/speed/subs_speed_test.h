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

#ifndef RELXILL_TEST_SPEED_SUBS_SPEED_TEST_H_
#define RELXILL_TEST_SPEED_SUBS_SPEED_TEST_H_

#include "LocalModel.h"
#include "XspecSpectrum.h"

#include <chrono>

extern "C" {
#include "relutility.h"
}

void eval_local_model_param_range(ModelName model_name, XPar param, double pmin, double pmax, int npar);

void eval_model_energy_grid_change(ModelName model_name, double eup_min, double eup_max, int npar);

void eval_model_relat_param_changes(ModelName model_name, const int num_evaluations);

void eval_model_xillver_param_changes(ModelName model_name, const int num_evaluations);


#endif //RELXILL_TEST_SPEED_SUBS_SPEED_TEST_H_
