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

void eval_local_model_param_range(ModelName model_name, XPar param, double pmin, double pmax, int npar) {

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();


  LocalModel local_model(model_name);

  for (int ii = 0; ii < npar; ii++) {
    double param_value = (pmax - pmin) * (static_cast<double>(ii) / static_cast<double>(npar)) + pmin;
    local_model.set_par(param, param_value);
    local_model.eval_model(spec);
  }

}

void eval_model_energy_grid_change(ModelName model_name, double eup_min, double eup_max, int npar) {

  const double elo = 0.1;
  const size_t nbins = 2000;

  LocalModel local_model(model_name);

  for (int ii = 0; ii < npar; ii++) {
    double eup = (eup_max - eup_min) * (static_cast<double>(ii) / static_cast<double>(npar));
    DefaultSpec default_spec{elo, eup, nbins};
    XspecSpectrum spec = default_spec.get_xspec_spectrum();
    local_model.eval_model(spec);
  }
}


void eval_model_relat_param_changes(ModelName model_name, const int num_evaluations) {
  XPar rel_param = XPar::a;
  eval_local_model_param_range(model_name, rel_param, 0.5, 0.998, num_evaluations);
}

void eval_model_xillver_param_changes(ModelName model_name, const int num_evaluations) {
  XPar rel_param = XPar::afe;
  // only choose a small difference here such that the table does not need to be re-loaded
  eval_local_model_param_range(model_name, rel_param, 3.1, 3.2, num_evaluations);
}
