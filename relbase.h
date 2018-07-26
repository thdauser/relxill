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

#include "common.h"

#include "fitsio.h"

#include "relutility.h"
#include "reltable.h"
#include "rellp.h"
#include "xilltable.h"


/*********** DEFINE STATEMENTS *********/

#define PARAM_DEFAULT 0.0

#define version_major 1
#define version_minor 1
#define version_build 0
#define version_dev "pre-9"

// TODO: allow to set it by an environment variable
/** path to all RELXILL tables */
#define RELXILL_TABLE_PATH "./"

/** dimensions of the RELLINE table */
#define RELTABLE_NA 25
#define RELTABLE_NR 100
#define RELTABLE_NG 40
#define RELTABLE_MAX_R 1000.0
#define RELTABLE_FILENAME "rel_table_v0.5a.fits"
#define RELTABLE_NMU0 30


/** dimensions of the LP table */
#define LPTABLE_NA 20
#define LPTABLE_NH 250
#define LPTABLE_FILENAME "rel_lp_table_v0.5b.fits"
#define LPTABLE_NR 100


/** parameters for interpolation an interagration **/
#define N_FRAD 1000      // values of radial bins (from rmin to rmax)
#define N_ZONES 10       // number of radial zones (as each zone is convolved with the input spectrum N_ZONES < N_FRAD)
#define N_ZONES_MAX 100  // maximal number of radial zones

/** parameters for the convolution **/
#define N_ENER_CONV  4096  // number of bins for the convolution, not that it needs to follow 2^N because of the FFT
#define EMIN_RELXILL 0.00035  // minimal energy of the convolution (in keV)
#define EMAX_RELXILL 2000.0 // minimal energy of the convolution (in keV)
#define EMIN_XILLVER 0.01
#define EMAX_XILLVER 1000.0

/** minimal and maximal energy for reflection strength calculation **/
#define RSTRENGTH_EMIN 20.0
#define RSTRENGTH_EMAX 40.0

/** file to (additionally) store the version number **/

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

	int limb_law;

	int ng;
	double** trff;
	double** cosne;
	double* gstar;
} str_relb_func;

typedef struct{
	int n_ener;
	double* ener;
	double* flux;
} out_spec;

typedef struct{
	int nzones;   // number of zones actually stored there
	int n_cache;  // number of array (nzones <= n_cache !!)
	int n_ener;
	double*** fft_xill;  // dimensions [n_cache,2,n_ener]
	double*** fft_rel;   // dimensions [n_cache,2,n_ener]
	xill_spec** xill_spec;
	out_spec* out_spec;
} specCache;

/****** FUNCTION DEFINITIONS ******/

/* the relbase function calculating the basic relativistic line shape for a given parameter setup*/
rel_spec* relbase(double* ener, const int n_ener,relParam* param, xillTable* xill_tab, int* status);

/** calculate the relline profile(s) for all given zones **/
void relline_profile(rel_spec* spec, relSysPar* sysPar, int* status);

void save_relline_profile(rel_spec* spec);

void free_cached_tables(void );

relSysPar* new_relSysPar(int nr, int ng, int* status);

void free_relSysPar(relSysPar* sysPar);

void relxill_kernel(double* ener_inp, double* spec_inp, int n_ener_inp, xillParam* xill_param, relParam* rel_param, int* status );

void relconv_kernel(double* ener_inp, double* spec_inp, int n_ener_inp, relParam* rel_param, int* status );

/** function adding a primary component with the proper norm to the flux **/
void add_primary_component(double* ener, int n_ener, double* flu, relParam* rel_param, xillParam* xill_param, int* status);

void free_rel_spec(rel_spec* spec);
rel_spec* new_rel_spec(int nzones, const int n_ener, int*status);

/** caching routines **/
void init_specCache(specCache** spec, int* status);
void free_specCache(void);
void free_fft_cache(double*** sp,int n1, int n2);
void free_out_spec(out_spec* spec);
out_spec* init_out_spec(int n_ener, double* ener, int* status);
int redo_relbase_calc(relParam* param, relParam* ca_para);
int redo_xillver_calc(relParam* rel_param, xillParam* xill_param, relParam* ca_rel, xillParam* ca_xill);
int comp_xill_param(xillParam* cpar, xillParam* par);
void set_cached_xill_param(xillParam* par,xillParam** ca_par, int* status);
void set_cached_rel_param(relParam* par, relParam** ca_rel_param, int* status);


#endif /* RELBASE_H_ */
