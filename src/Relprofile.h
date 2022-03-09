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
#ifndef RELPROFILE_H_
#define RELPROFILE_H_

extern "C" {
#include "common.h"
}

void calc_relline_profile(relline_spec_multizone *spec, RelSysPar *sysPar, int *status);

RelSysPar *new_relSysPar(int nr, int ng, int *status);

void free_relSysPar(RelSysPar *sysPar);

RelSysPar *get_system_parameters(relParam *param, int *status);

void renorm_relline_profile(relline_spec_multizone *spec, relParam *rel_param, const int *status);

void init_relline_spec_multizone(relline_spec_multizone **spec, relParam *param, xillTable *xill_tab, const double *radial_zones,
                                 double **pt_ener, int n_ener, int *status);

void free_cached_relTable(void);
void free_relprofile_cache(void);
void free_cache_syspar(void);
#endif
