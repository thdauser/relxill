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

    Copyright 2016 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELBASE_H_
#define RELBASE_H_

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <float.h>

#include "fitsio.h"

#include "common.h"

#include "relutility.h"
#include "reltable.h"
#include "rellp.h"


/*********** DEFINE STATEMENTS *********/

#define PARAM_DEFAULT -DBL_MAX

#define version_major 0
#define version_minor 1
#define version_build 0



/** path to all RELXILL tables */
#define RELXILL_TABLE_PATH "/home/thomas/data/relline_tables/"

/** dimensions of the RELLINE table */
#define RELTABLE_NA 30
#define RELTABLE_NMU0 22
#define RELTABLE_NR 100
#define RELTABLE_NG 20
#define RELTABLE_MAX_R 1000.0
#define RELTABLE_FILENAME "rel_table_v0.4e.fits"

/** parameters for interpolation an interagration **/
#define N_FRAD 1000      // values of radial bins (from rmin to rmax)
#define N_ZONES 10       // number of radial zones (as each zone is convolved with the input spectrum N_ZONES < N_FRAD)


/****** TYPE DEFINITIONS ******/

typedef struct{
	int save_g_ind;


	double cache_rad_relb_fun;
	double cache_val_relb_func[2];

	double cache_bin_ener;
	int cached_relbf;


	double re;
	double gmax;
	double gmin;
	double del_g;
	double emis;

	int ng;
	double** trff;
	double** cosne;
	double* gstar;
} str_relb_func;

/****** FUNCTION DEFINITIONS ******/

/* the relbase function calculating the basic relativistic line shape for a given parameter setup*/
void relbase(double* ener, const int n_ener, double* photar, relParam* param, int* status);

/** calculate the relline profile(s) for all given zones **/
void relline_profile(rel_spec* spec, relSysPar* sysPar, int* status);

void free_cached_tables(void );

relSysPar* new_relSysPar(int nr, int ng, int* status);

void free_relSysPar(relSysPar* sysPar);

#endif /* RELBASE_H_ */
