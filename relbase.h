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

    Copyright 2020 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELBASE_H_
#define RELBASE_H_

#define _GNU_SOURCE  // seems necessary for standard Makefile

#include "common.h"

#include "relutility.h"
#include "reltable.h"
#include "rellp.h"
#include "xilltable.h"
#include "relcache.h"


/*********** DEFINE STATEMENTS *********/

#define PARAM_DEFAULT 0.0

#define version_major 1
#define version_minor 3
#define version_build 11
#define version_dev "-1"

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

/** calculate the relline profile(s) for all given zones **/
void relline_profile(rel_spec *spec, relSysPar *sysPar, int *status);

void save_relline_profile(rel_spec *spec, int* status);

void save_radial_profile(char *foutName, double *rad, double *intens, int n_rad);

rel_spec *relbase_multizone(double *ener,
                            const int n_ener,
                            relParam *param,
                            xillTable *xill_tab,
                            double *radialZones,
                            int nzones,
                            int *status);

void free_cached_tables(void);

relSysPar *new_relSysPar(int nr, int ng, int *status);

void free_relSysPar(relSysPar *sysPar);

void relxill_kernel(double *ener_inp,
                    double *spec_inp,
                    int n_ener_inp,
                    xillParam *xill_param,
                    relParam *rel_param,
                    int *status);

void relconv_kernel(double *ener_inp, double *spec_inp, int n_ener_inp, relParam *rel_param, int *status);

relSysPar *get_system_parameters(relParam *param, int *status);

/** function adding a primary component with the proper norm to the flux **/
void add_primary_component(double *ener,
                           int n_ener,
                           double *flu,
                           relParam *rel_param,
                           xillParam *xill_param,
                           int *status);

void free_rel_spec(rel_spec *spec);
rel_spec *new_rel_spec(int nzones, const int n_ener, int *status);

/** caching routines **/
specCache *init_globalSpecCache(int *status);
void free_specCache(void);
void free_fft_cache(double ***sp, int n1, int n2);
void free_out_spec(out_spec *spec);

out_spec *init_out_spec(int n_ener, const double *ener, int *status);

int redo_xillver_calc(relParam *rel_param, xillParam *xill_param, relParam *ca_rel, xillParam *ca_xill);
int redo_relbase_calc(relParam *rel_param, relParam *ca_rel_param);

void set_cached_rel_param(relParam *par, relParam **ca_rel_param, int *status);

int comp_xill_param(xillParam *cpar, xillParam *par);

/** free the CLI cache **/
void free_cache(void);

void get_xillver_angdep_spec(double *o_xill_flux,
                             int n_ener,
                             double *ener,
                             double *rel_dist,
                             xill_spec *xill_spec,
                             int *status);

void convolveSpectrumFFTNormalized(double *ener, const double *fxill, const double *frel, double *fout, int n,
                                   int re_rel, int re_xill, int izone, specCache *local_spec_cache, int *status);

void get_std_relxill_energy_grid(int *n_ener, double **ener, int *status);

void renorm_xill_spec(double *spec, int n, double lxi, double dens);

double calcNormWrtXillverTableSpec(const double *flux, const double *ener, const int n, int *status);

void set_stdNormXillverEnerygrid(int *status);
EnerGrid *get_stdXillverEnergygrid(int *status);

#endif /* RELBASE_H_ */
