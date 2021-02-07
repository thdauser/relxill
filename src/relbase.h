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
#ifndef RELBASE_H_
#define RELBASE_H_

#ifndef _GNU_SOURCE
#define _GNU_SOURCE  // seems necessary for standard Makefile
#endif

#include "config.h"

#include "common.h"

#include "relutility.h"
#include "reltable.h"
#include "rellp.h"
#include "xilltable.h"
#include "relcache.h"
#include "relprofile.h"


/*********** DEFINE STATEMENTS *********/

#define PARAM_DEFAULT 0.0

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

/** dimensions of the RR table */
#define RETURNRAD_TABLE_NR 50
#define RETURNRAD_TABLE_NG 20
#define RETURNRAD_TABLE_FILENAME "table_returnRad_v20200506.fits"


/** parameters for interpolation an interagration **/
#define N_FRAD 1000      // values of radial bins (from rmin to rmax)
#define N_ZONES 10       // number of radial zones (as each zone is convolved with the input spectrum N_ZONES < N_FRAD)
#define N_ZONES_MAX 50  // maximal number of radial zones

/** parameters for the convolution **/
#define N_ENER_CONV  4096  // number of bins for the convolution, not that it needs to follow 2^N because of the FFT
#define EMIN_RELXILL 0.00035  // minimal energy of the convolution (in keV)
#define EMAX_RELXILL 2000.0 // minimal energy of the convolution (in keV)
#define EMIN_XILLVER_NORMALIZATION 0.1
#define EMAX_XILLVER_NORMALIZATION 1000.0
#define EMIN_XILLVER 0.01
#define EMAX_XILLVER EMAX_XILLVER_NORMALIZATION
#define N_ENER_XILLVER 3000

/** minimal and maximal energy for reflection strength calculation **/
#define RSTRENGTH_EMIN 20.0
#define RSTRENGTH_EMAX 40.0

/****** TYPE DEFINITIONS ******/

typedef struct {
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
  double **trff;
  double **cosne;
  double *gstar;
} str_relb_func;

/****** FUNCTION DEFINITIONS ******/

/* get the current version number */
void get_version_number(char **vstr, int *status);

/* the relbase function calculating the basic relativistic line shape for a given parameter setup*/
rel_spec *relbase(double *ener, const int n_ener, relParam *param, xillTable *xill_tab, int *status);



rel_spec *relbase_multizone(double *ener,
                            const int n_ener,
                            relParam *param,
                            xillTable *xill_tab,
                            double *radialZones,
                            int nzones,
                            int *status);

void free_cached_tables(void);


void relconv_kernel(double *ener_inp, double *spec_inp, int n_ener_inp, relParam *rel_param, int *status);


/** function adding a primary component with the proper norm to the flux **/
void add_primary_component(double *ener,
                           int n_ener,
                           double *flu,
                           relParam *rel_param,
                           xillParam *xill_param,
                           int *status);

void free_rel_spec(rel_spec *spec);
rel_spec *new_rel_spec(int nzones, const int n_ener, int *status);

void fft_conv_spectrum(double *ener, const double *fxill, const double *frel, double *fout, int n,
                       int re_rel, int re_xill, int izone, specCache *cache, int *status);
double calcFFTNormFactor(const double *ener, const double *fxill, const double *frel, const double *fout, int n);

/** caching routines **/
specCache *init_global_specCache(int *status);
void free_specCache(specCache* spec_cache);
void free_fft_cache(double ***sp, int n1, int n2);
void free_out_spec(OutSpec *spec);

OutSpec *init_out_spec(int n_ener, const double *ener, int *status);

int redo_xillver_calc(relParam *rel_param, xillParam *xill_param, relParam *ca_rel, xillParam *ca_xill);
int redo_relbase_calc(relParam *rel_param, relParam *ca_rel_param);

void set_cached_rel_param(relParam *par, relParam **ca_rel_param, int *status);

int did_xill_param_change(xillParam *cpar, xillParam *par);

void free_cache(void);


void convolveSpectrumFFTNormalized(double *ener, const double *fxill, const double *frel, double *fout, int n,
                                   int re_rel, int re_xill, int izone, specCache *local_spec_cache, int *status);

void get_std_relxill_energy_grid(int *n_ener, double **ener, int *status);

void renorm_xill_spec(float *spec, int n, double lxi, double dens);

double calcNormWrtXillverTableSpec(const double *flux, const double *ener, const int n, int *status);

void set_stdNormXillverEnerygrid(int *status);
EnerGrid *get_stdXillverEnergygrid(int *status);

#endif /* RELBASE_H_ */
