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

// assumes FFTW3 installation in heasoft
#if __has_include("fftw3.h")  // quotes needed for GCC < 10
#  include "fftw3.h"
#else
#  include "fftw/fftw3.h"
#endif



/*********** DEFINE STATEMENTS *********/

// make sure M_PI is defined, also with old/strange compiler options
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


/** define Emissivity Model Type **/
#define EMIS_TYPE_NONE 0
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


/**** DEFINES **/
#define MOD_TYPE_RELLINE 1
#define MOD_TYPE_RELLINELP 2

#define MOD_TYPE_RELCONV 11
#define MOD_TYPE_RELCONVLP 12

#define MOD_TYPE_XILLVER 0
#define MOD_TYPE_XILLVER_NTHCOMP 100

#define MOD_TYPE_RELXILL -1
#define MOD_TYPE_RELXILLLP -2

/** density models **/
#define MOD_TYPE_RELXILLDENS -10
#define MOD_TYPE_RELXILLLPDENS -11
#define MOD_TYPE_XILLVERDENS -100

/** ion grad models **/
#define MOD_TYPE_RELXILLLPION -21

/** CO models **/
#define MOD_TYPE_RELXILLCO -200
#define MOD_TYPE_XILLVERCO -210

/** Neutron Star / BB models **/
#define MOD_TYPE_RELXILLNS -30
#define MOD_TYPE_XILLVERNS -101

// unpublished models
#define MOD_TYPE_RELXILLBBRET -300
#define MOD_TYPE_RELXILLALPHA -2000
#define MOD_TYPE_RELXILLLPALPHA -2001



#define PARAM_DEFAULT 0.0

/** path to all RELXILL tables */
#define RELXILL_TABLE_PATH "./"


/** dimensions of the RR table */
#define RETURNRAD_TABLE_NR 50
#define RETURNRAD_TABLE_NG 20
#define RETURNRAD_TABLE_FILENAME "table_returnRad_v20220301.fits"


/** parameters for interpolation an interagration **/
#define N_FRAD 1000      // values of radial bins (from rmin to rmax)
#define N_ZONES 10       // number of radial zones (as each zone is convolved with the input spectrum N_ZONES < N_FRAD)
#define N_ZONES_IONGRAD 25  // default number of radial zones for iongrad models
#define N_ZONES_MAX 50  // maximal number of radial zones


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
#define XILLTABLE_NTHCOMP_FILENAME "xillverCp_v3.4.fits"
#define XILLTABLE_NS_FILENAME "xillverNS-2.fits"
#define XILLTABLE_CO_FILENAME "xillverCO.fits"

// useful constants
#define CONVERT_KEV2ERG 1.6021773000008302e-09  // 1 keV in erg

enum xillTableIds {
  XILLTABLE_ID_STANDARD,
  XILLTABLE_ID_NTHCOMP,
  XILLTABLE_ID_NS,
  XILLTABLE_ID_CO,
};

#define AD_ROUT_MAX 1000

/****** TYPE DEFINITIONS ******/

typedef struct{
  double* rgrid;
  int n_zones;
  double* corrfac_flux;
  double* corrfac_gshift;
} rradCorrFactors;

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
  int ion_grad_type;
  rradCorrFactors* rrad_corr_factors;
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
  int prim_type;
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
  double boost;
  double shiftTmaxRRet; // temperature shift of Tmax, this should not be a free parameter in the end
  double mass_msolar;
  double distance;
  double luminosity_primary_source;
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
  int *param_index;  // length is N_PARAM_MAX
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
  double f_inf_rest;  // f_inf in the rest frame of the system (not taking beta of the LP into account)
} lpReflFrac;

/** the emissivity profile (del and del_inc are only of interest in the LP geometry) **/
typedef struct {

  double *re;
  int nr;
  double *emis;       // intensity on the surface of the accretion disc
  double *del_emit;   // angle under which the photon is emitted from the primary source
  double *del_inc;    // angle the photon hits the accretion disk (in the rest frame of the disk)

  lpReflFrac *photon_fate_fractions;
  double
      normFactorPrimSpec;  // determined from the f_inf_rest and g_inf, the factor to multiply the direct radiation with

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
} spectrum;

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
  spectrum *out_spec;
} specCache;

typedef struct {

  double *ener;
  int nbins;

} EnerGrid;

typedef struct {

  int nbins;
  double *ener;  // has num_flux_bins+1
  double *flux; // has length

} TestSpectrum;


/******************************/
/* define the c_donthcomp function here */
void c_donthcomp(const double *ear, int ne, double *param, double *photar);

#endif /* COMMON_H_ */
