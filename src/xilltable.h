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
#ifndef XILLTABLE_H_
#define XILLTABLE_H_

#include "relutility.h"
#include "common.h"

// currently the number of different parameters that can be given in a table
#define N_PARAM_MAX 10  // has to be as long as the NAME_ array
#define PARAM_GAM 0
#define PARAM_AFE 1
#define PARAM_ACO 1  // caveat: internally we use it same as AFE
#define PARAM_LXI 2
#define PARAM_ECT 3
#define PARAM_KTE 3  // caveat: internally we treat kTe as Ecut
#define PARAM_DNS 4
#define PARAM_KTB 5
#define PARAM_FRA 6
#define PARAM_INC 7

#define NAME_GAM "Gamma"
#define NAME_AFE "A_Fe"
#define NAME_LXI "logXi"
#define NAME_ECT "Ecut"
#define NAME_KTE "kTe"
#define NAME_KTB "kTbb"
#define NAME_DNS "Dens"
#define NAME_ACO "A_CO"
#define NAME_FRA "Frac"
#define NAME_INC "Incl"


/** name of the XILLVER table */
#define XILLTABLE_FILENAME "xillver-a-Ec5.fits"
#define XILLTABLE_NTHCOMP_FILENAME "xillverCp_v3.3.fits"
#define XILLTABLE_NS_FILENAME "xillverNS-2.fits"
#define XILLTABLE_CO_FILENAME "xillverCO.fits"

/** get a new and empty rel table (structure will be allocated)  */
xillTable *new_xillTable(int num_param, int *status);

/* destroy the relline table structure */
void free_xillTable(xillTable *tab);

xillSpec *interp_xill_table(xillTable *tab, xillParam *param, const int *ind, int *status);

int get_xilltab_param_index(xillTable *tab, int ind);
float *get_xilltab_paramvals(xillParam *param, int *status);
int *get_xilltab_indices_for_paramvals(xillParam *param, xillTable *tab, int *status);

void check_xilltab_cache(const char *fname, xillParam *param, xillTable *tab, const int *ind, int *status);

void free_cached_xillTable(void);

void init_xillver_table(const char *filename, xillTable **inp_tab, int *status);

const char *get_init_xillver_table(xillTable **tab, xillParam *param, int *status);

void print_xilltable_parameters(const xillTable *tab, char *const *xilltab_parname);


fitsfile *open_fits_table_stdpath(const char *filename, int *status);

int checkIfTableExists(const char *filename, int *status);

int is_6dim_table(int model_type);

int is_ns_model(int model_type);
int is_co_model(int model_type);

#endif /* XILLTABLE_H_ */
