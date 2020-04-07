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

    Copyright 2019 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef XILLTABLE_H_
#define XILLTABLE_H_

#include "relbase.h"
#include "relutility.h"
#include "relmodels.h"

#define XILLTABLE_N_PARAM 5

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
#define XILLTABLE_DENS_FILENAME "xillverD-5.fits"
#define XILLTABLE_NTHCOMP_FILENAME "xillver-comp.fits"
#define XILLTABLE_NS_FILENAME "xillverNS-2.fits"
#define XILLTABLE_CO_FILENAME "xillverCO.fits"

/** get a new and empty rel table (structure will be allocated)  */
xillTable *new_xillTable(int num_param, int *status);

/* destroy the relline table structure */
void free_xillTable(xillTable* tab);

/** the main routine for the xillver table: returns a spectrum for the given parameters
 *  (decides if the table needs to be initialized and/or more data loaded          */
xill_spec* get_xillver_spectra(xillParam* param, int* status);

xill_spec* new_xill_spec(int n_incl, int n_ener, int* status);
void free_xill_spec(xill_spec* spec);

void free_cached_xillTable(void);

void init_xillver_table(char *filename, xillTable **inp_tab, xillParam *param, int *status);

char* get_init_xillver_table(xillTable** tab, xillParam* param, int* status);

void print_xilltable_parameters(const xillTable *tab, char *const *xilltab_parname);

void norm_xillver_spec(xill_spec *spec, double incl);

int is_6dim_table(int model_type);

#endif /* XILLTABLE_H_ */
