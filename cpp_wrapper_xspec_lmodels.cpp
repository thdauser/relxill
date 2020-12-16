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

#include "cppmodels.h"
#include "cppparameters.h"

extern "C"
void lmodcpprelline(const double *energy, int Nflux, const double *parameter,
                    int spectrum, double *flux, double *fluxError, const char *init) {
  xspec_C_wrapper_eval_model(ModelName::relline, parameter, flux, Nflux, energy);
}

extern "C"
void lmodcpprelxill(const double *energy, int Nflux, const double *parameter,
                    int spectrum, double *flux, double *fluxError, const char *init) {
  xspec_C_wrapper_eval_model(ModelName::relxill, parameter, flux, Nflux, energy);
}
