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

#include "common-functions.h"

double sum_flux(const double *flux, int nbins) {

  double sum = 0.0;
  for (int ii = 0; ii < nbins; ii++) {
    sum += flux[ii];
  }
  return sum;
}

void eval_xspec_lmod_default(ModelName model_name, const DefaultSpec& default_spec){

  const double *xspec_parameters = get_xspec_default_parameter_array(model_name);

  xspec_C_wrapper_eval_model(model_name,
                               xspec_parameters,
                               default_spec.flux,
                               default_spec.num_flux_bins,
                               default_spec.energy);

  delete[] xspec_parameters;

}
