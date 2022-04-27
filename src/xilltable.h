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

    Copyright 2022 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef XILLTABLE_H_
#define XILLTABLE_H_

#include "relutility.h"
#include "common.h"



void renorm_xill_spec(float *spec, int n, double lxi, double dens);

enum xillTableIds get_xilltable_id(int model_id, int prim_type);

xillTableParam *get_xilltab_param(xillParam *param, int *status);

/** get a new and empty rel table (structure will be allocated)  */
xillTable *new_xillTable(int num_param, int *status);

/* destroy the relline table structure */
void free_xillTable(xillTable *tab);

xillSpec *interp_xill_table(xillTable *tab, xillTableParam *param, const int *ind, int *status);

int get_xilltab_param_index(xillTable *tab, int ind);
float *get_xilltab_paramvals(const xillTableParam *param, int *status);
int *get_xilltab_indices_for_paramvals(const xillTableParam *param, xillTable *tab, int *status);

void check_xilltab_cache(const char *fname, const xillTableParam *param, xillTable *tab, const int *ind, int *status);

void free_cached_xillTable(void);

void init_xillver_table(const char *filename, xillTable **inp_tab, int *status);

const char *get_init_xillver_table(xillTable **tab, int model_type, int prim_type, int *status);

void print_xilltable_parameters(const xillTable *tab, char *const *xilltab_parname);


fitsfile *open_fits_table_stdpath(const char *filename, int *status);

int checkIfTableExists(const char *filename, int *status);

int is_6dim_table(int model_type);

int is_ns_model(int model_type);
int is_co_model(int model_type);

void free_xill_spec(xillSpec *spec);

#endif /* XILLTABLE_H_ */
