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

#ifndef RELXILL_TEST_UNIT_COMMON_FUNCTIONS_H_
#define RELXILL_TEST_UNIT_COMMON_FUNCTIONS_H_

#include "XspecSpectrum.h"
#include "LocalModel.h"


double sum_flux(const double *flux, int nbins);

void eval_xspec_lmod_default(ModelName model_name, const DefaultSpec& default_spec);

#endif //RELXILL_TEST_UNIT_COMMON_FUNCTIONS_H_
