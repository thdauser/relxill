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

#include "relutility.h"
#include "reltable.h"


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
#define RELTABLE_FILENAME "rel_table_v0.4e.fits"



/****** TYPE DEFINITIONS ******/

typedef struct{
	int model_type;
	double a;
	double incl;
	double emis1;
	double emis2;
	double rbr;
	double rin;
	double rout;
	double lineE;
	double z;
	double height;
	double gamma;
} relParam;




/****** FUNCTION DEFINITIONS ******/

/* the relbase function calculating the basic relativistic line shape for a given parameter setup*/
void relbase(const double* ener, const int n_ener, double* photar, const relParam* param, int* status);

void free_cached_tables(void );


#endif /* RELBASE_H_ */
