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

#ifndef XILLSPEC_H_
#define XILLSPEC_H_

#include "relutility.h"
#include "common.h"

double norm_factor_semi_infinite_slab(double incl_deg);
void norm_xillver_spec(xillSpec *spec, double incl);

/** the main routine for the xillver table: returns a spectrum for the given parameters
 *  (decides if the table needs to be initialized and/or more data loaded          */
xillSpec *get_xillver_spectra(xillParam *param, int *status);


xillSpec *new_xill_spec(int n_incl, int n_ener, int *status);
void free_xill_spec(xillSpec *spec);

#endif //XILLSPEC_H_
