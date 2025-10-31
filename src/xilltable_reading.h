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

#include "common.h"



void renorm_xill_spec(float *spec, int n, double lxi, double dens);

enum xillTableIds get_xilltable_id(int model_id, int prim_type);

int get_xilltab_param_index(const xillTable* tab, int unique_param_identifier);
xillTableParam* get_xilltab_param(const xillParam* param, int* status);

/** get a new and empty rel table (structure will be allocated)  */
xillTable *new_xillTable(int num_param, int *status);

/* destroy the relline table structure */
void free_xillTable(xillTable** tab);

void free_all_cached_xillTables(void);

void init_xillver_table(const char *filename, xillTable **inp_tab, int *status);

const char *get_init_xillver_table(xillTable **tab, int model_type, int prim_type, int *status);

void print_xilltable_parameters(const xillTable *tab, char *const *xilltab_parname);


fitsfile *open_fits_table_stdpath(const char *filename, int *status);

int checkIfTableExists(const char *filename, int *status);

int is_6dim_table(int model_type);

void xilltable_fits_load_single_spec(const char* fname,
                                     fitsfile** fptr,
                                     const xillTable* tab,
                                     double defDensity,
                                     double defLogxi,
                                     int nn,
                                     int ii,
                                     int jj,
                                     int kk,
                                     int ll,
                                     int mm,
                                     int* status);


float* get_xillspec_table(const xillTable* tab, int i0, int i1, int i2, int i3, int i4, int i5);

#endif /* XILLTABLE_H_ */
