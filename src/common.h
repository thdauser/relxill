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
#ifndef COMMON_H_
#define COMMON_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <fitsio.h>

#include "complex.h" // needed by fftw3 (has to be included before fftw3.h)
#include "fftw/fftw3.h"   // assumes installation in heasoft


/*********** DEFINE STATEMENTS *********/

/** define Emissivity Model Type **/
#define EMIS_TYPE_BKN 1
#define EMIS_TYPE_LP 2
#define EMIS_TYPE_ALPHA 3
#define EMIS_TYPE_CONST 4
/***************************************/

/** define primary spectrum Type **/
#define PRIM_SPEC_NONE 0
#define PRIM_SPEC_ECUT 1
#define PRIM_SPEC_NTHCOMP 2
#define PRIM_SPEC_BB 3
/***************************************/

/** define Ionization Gradient Type**/
#define ION_GRAD_TYPE_CONST 0
#define ION_GRAD_TYPE_PL 1
#define ION_GRAD_TYPE_ALPHA 2
/***************************************/


enum xillTableIds {
  XILLTABLE_ID_STANDARD,
  XILLTABLE_ID_NTHCOMP,
  XILLTABLE_ID_NS,
  XILLTABLE_ID_CO,
};

#define AD_ROUT_MAX 1000

/****** TYPE DEFINITIONS ******/

typedef struct {
  int model_type;
  int emis_type;
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
  double d_offaxis;
  double htop;
  double gamma;
  double beta;
  int limb;
  int do_renorm_relline;
  int num_zones;
  int return_rad;
  double
      return_rad_flux_correction_factor; // calculated ratio from xillver reflection to primary spectrum (energy flux)
  double xillver_gshift_corr_fac; // correction factor of flux boost from xillver wrt to g^Gamma (for g=1.5)
} relParam;

typedef struct {
  double gam;
  double afe;
  double lxi;
  double ect;
  double incl;
  double dens;
  double frac_pl_bb;
  double kTbb;
  enum xillTableIds xilltable_id;
  int model_type;
} xillTableParam;

typedef struct {
  double gam;
  double afe;
  double lxi;
  double ect;
  double incl;
  double dens;
  double frac_pl_bb;
  int model_type;
  int prim_type;
  double z;
  double refl_frac;
  int interpret_reflfrac_as_boost;
  double iongrad_index;
  double kTbb;
  int ion_grad_type;
  double boost;
  double shiftTmaxRRet; // temperature shift of Tmax, this should not be a free parameter in the end
} xillParam;


/** the XILLVER table structure */
typedef struct {

  float *elo;
  float *ehi;
  int n_ener;

  int num_param;        // number of parameters (basically the dimension of the table)
  int *num_param_vals;   // number of values given for each parameter

  // information on the inclination is stored separately (the last parameter in the table)
  int n_incl;
  float *incl;

  /* need to identify the meaning of each parameter here [index in the array]
   * (see routine "get_xilltab_paramvals" */
  int *param_index;  // lenth is N_PARAM_MAX
  char **param_names; // name of each parameter

  float **param_vals;    // array to store the parameter values (as given in the table)

  float **data_storage;   // storage of a n-dim table (n_elements spectra with n_ener bins each)
  int num_elements;

} xillTable;

typedef struct {
  double refl_frac;
  double f_bh;
  double f_ad;
  double f_inf;
} lpReflFrac;

/** the emissivity profile (del and del_inc are only of interest in the LP geometry) **/
typedef struct {

  double *re;
  int nr;
  double *emis;       // intensity on the surface of the accretion disc
  double *del_emit;   // angle under which the photon is emitted from the primary source
  double *del_inc;    // angle the photon hits the accretion disk (in the rest frame of the disk)

  lpReflFrac *photon_fate_fractions;
  double normFactorPrimSpec;  // determined from the f_inf and g_inf, the factor to multiply the direct radiation with

} emisProfile;

typedef struct {
  int nr;
  int ng;

  double *re;
  double *gmin;
  double *gmax;
  double *gstar;
  double *d_gstar;  // bin width for each gstar value

  double ***trff;
  double ***cosne;

  emisProfile *emis;

  double del_ad_max;  // delta of the photon where it would hit 1000rg (the outer edge of the disk in the relline table)

  int limb_law;

} RelSysPar;

/** angles (cosne) and their distribution over the radial zones **/
typedef struct {
  int n_cosne;
  int n_zones;
  double *cosne;
  double **dist;      // [n_zones][n_cosne]
} RelCosne;

typedef struct {
  int n_ener;
  int n_zones;
  double *rgrid;      // length=n_zones + 1
  double *ener;       // length=n_ener +1
  double **flux;      // [n_zones][n_ener]
  RelCosne *rel_cosne;
} relline_spec_multizone;

typedef struct {
  double *ener;  // has n_ener+1 elements
  double *incl;  // has n_incl elements
  double **flu;  // [n_incl,n_ener+1]
  int n_ener;
  int n_incl;
} xillSpec;

typedef struct {
  int n_ener;
  double *ener;
  double *flux;
} OutSpec;

typedef struct {
  int nzones;   // number of zones actually stored there
  int n_cache;  // number of array (nzones <= n_cache !!)
  int n_ener;
  double* conversion_factor_energyflux; // conversion from photons/bin to keV/keV
  double ***fft_xill;  // dimensions [n_cache,2,n_ener]
  double ***fft_rel;   // dimensions [n_cache,2,n_ener]

  fftw_complex** fftw_xill;  // dimensions [n_cache,n_ener]
  fftw_complex** fftw_rel;   // dimensions [n_cache,n_ener]

  fftw_complex* fftw_backwards_input;  // [nener]
  double* fftw_output;  // [nener]
  fftw_plan plan_c2r;


  xillSpec **xill_spec;
  OutSpec *out_spec;
} specCache;

typedef struct {

  double *ener;
  int nbins;

} EnerGrid;

typedef struct {

  int nbins;
  double *ener;  // has nbins+1
  double *flux; // has length

} Spectrum;


/******************************/
/* define the c_donthcomp function here */
void c_donthcomp(double *ear, int ne, double *param, double *photar);

#endif /* COMMON_H_ */
