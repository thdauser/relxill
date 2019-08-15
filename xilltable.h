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
#define XILLTABLE_N_GAM 13
#define XILLTABLE_N_AFE 4
#define XILLTABLE_N_LXI 15
#define XILLTABLE_N_ECT 11
#define XILLTABLE_N_INCL 10

#define XILLTABLE_DENS_N_GAM 13
#define XILLTABLE_DENS_N_AFE 5
#define XILLTABLE_DENS_N_LXI 15
#define XILLTABLE_DENS_N_DENS 9
#define XILLTABLE_DENS_N_INCL 10

#define XILLTABLE_NTHCOMP_N_GAM 13
#define XILLTABLE_NTHCOMP_N_AFE 4
#define XILLTABLE_NTHCOMP_N_LXI 15
#define XILLTABLE_NTHCOMP_N_KTE 12
#define XILLTABLE_NTHCOMP_N_INCL 10

// #define XILLTABLE_num_param_vals { XILLTABLE_N_GAM, XILLTABLE_N_AFE,  XILLTABLE_N_LXI, XILLTABLE_N_ECT,  XILLTABLE_N_INCL}


/** name of the XILLVER table */
#define XILLTABLE_FILENAME "xillver-a-Ec5.fits"
#define XILLTABLE_DENS_FILENAME "xillverD-5.fits"
#define XILLTABLE_NTHCOMP_FILENAME "xillver-comp.fits"

/** get a new and empty rel table (structure will be allocated)  */
xillTable* new_xillTable(int n_gam, int n_afe, int n_lxi, int n_ect, int n_incl, int* status);

/* destroy the relline table structure */
void free_xillTable(xillTable* tab);

/** the main routine for the xillver table: returns a spectrum for the given parameters
 *  (decides if the table needs to be initialized and/or more data loaded          */
xill_spec* get_xillver_spectra(xillParam* param, int* status);

xill_spec* new_xill_spec(int n_incl, int n_ener, int* status);
void free_xill_spec(xill_spec* spec);

void free_cached_xillTable(void);

char* get_init_xillver_table(xillTable** tab, xillParam* param, int* status);

void norm_xillver_spec(xill_spec* spec, double incl);

#endif /* XILLTABLE_H_ */
