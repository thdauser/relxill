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
#ifndef RELBASE_H_
#define RELBASE_H_

#ifndef _GNU_SOURCE
#define _GNU_SOURCE  // seems necessary for standard Makefile
#endif

#include "Rellp.h"
#include "Relreturn_Corona.h"
#include "Relcache.h"
#include "Relprofile.h"
#include "Xillspec.h"

extern "C" {
#include "xilltable.h"
#include "config.h"
#include "common.h"
#include "relutility.h"
#include "reltable.h"
}


/*********** DEFINE STATEMENTS *********/

/** minimal and maximal energy for reflection strength calculation **/
#define RSTRENGTH_EMIN 20.0
#define RSTRENGTH_EMAX 40.0

// minimal and maximal values allowed for the convolution
#define RELCONV_EMIN 0.01
#define RELCONV_EMAX 1000.0

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
relline_spec_multizone *relbase(double *ener, const int n_ener, relParam *param, int *status);

relline_spec_multizone* relbase_profile(double *ener, int n_ener, relParam *param,
                                        RelSysPar *sysPar,
                                        xillTable *xill_tab,
                                        const double *radialZones,
                                        int nzones,
                                        int *status);


relline_spec_multizone *relbase_multizone(double *ener,
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
                           xillParam *xill_input_param,
                           RelSysPar *sys_par,
                           int *status);

void free_rel_spec(relline_spec_multizone *spec);
relline_spec_multizone *new_rel_spec(int nzones, const int n_ener, int *status);

double calcFFTNormFactor(const double *ener, const double *fxill, const double *frel, const double *fout, int n);

/** caching routines **/
specCache *init_global_specCache(int *status);
void free_specCache(specCache *spec_cache);
void free_fft_cache(double ***sp, int n1, int n2);
void free_spectrum(spectrum *spec);

spectrum *new_spectrum(int n_ener, const double *ener, int *status);

int redo_xillver_calc(const relParam *rel_param, const xillParam *xill_param,
                      const relParam *ca_rel, const xillParam *ca_xill);
int redo_relbase_calc(const relParam *rel_param, const relParam *ca_rel_param);

void set_cached_rel_param(const relParam *par, relParam **ca_rel_param, int *status);

int did_xill_param_change(const xillParam *cpar, const xillParam *par);

void free_cache(void);

void convolveSpectrumFFTNormalized(double *ener, const double *fxill, const double *frel, double *fout, int n,
                                   int re_rel, int re_xill, int izone, specCache *local_spec_cache, int *status);


double calcXillverNormFromPrimarySpectrum(const double *flux, const double *ener, const int n, int *status);

void set_stdNormXillverEnerygrid(int *status);
// EnerGrid *get_coarse_xillver_energrid();

double *calc_normalized_xillver_primary_spectrum(const double *ener, int n_ener,
                                                 const relParam *rel_param, const xillTableParam *xill_param, int *status);
int convert_relxill_to_xillver_model_type(int relxill_model_type, int *status);

#endif //RELBASE_H_
