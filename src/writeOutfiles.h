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
#ifndef WRITEOUTFILES_H_
#define WRITEOUTFILES_H_

#include "common.h"

void write_binned_data_to_file(char *foutName, const double *rad, double *intens, int n_rad);
void write_data_to_file(const char *foutName, double *rad, double *intens, int n_rad);

void save_relline_radial_flux_profile(double *rad, double *intens, int n_rad);

void save_xillver_spectrum(const double *ener, double *flu, int n_ener, char *fname);
void save_relline_profile(rel_spec *spec);
void save_emis_profiles(RelSysPar *sysPar);

#endif
