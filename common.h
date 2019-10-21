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
#ifndef COMMON_H_
#define COMMON_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <float.h>

/*********** DEFINE STATEMENTS *********/

/** define Emissivity Model Type **/
#define EMIS_TYPE_BKN 1
#define EMIS_TYPE_LP 2
/***************************************/

/** define primary spectrum Type **/
#define PRIM_SPEC_NONE 0
#define PRIM_SPEC_ECUT 1
#define PRIM_SPEC_NTHCOMP 2
/***************************************/

/** define Ionization Gradient Type**/
#define ION_GRAD_TYPE_CONST 0
#define ION_GRAD_TYPE_PL 1
#define ION_GRAD_TYPE_ALPHA 2
/***************************************/


#define AD_ROUT_MAX 1000

/****** TYPE DEFINITIONS ******/

typedef struct{
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
	double gamma;
	double beta;
	int limb;
	int do_renorm_relline;
	int num_zones;
} relParam;

typedef struct{
	double gam;
	double afe;
	double lxi;
	double ect;
	double incl;
	double z;
	double refl_frac;
	double dens;
	double ion_grad_index;
	int ion_grad_type;
	int fixReflFrac;
	int model_type;
	int prim_type;
} xillParam;


typedef struct{
	double* lxi;
	double* fx;
	double* r;
	double* del_emit;
	int nbins;
} ion_grad;

/** the XILLVER table structure */
typedef struct{

	float* elo;
	float* ehi;
	int n_ener;

	int  num_param;        // number of parameters (basically the dimension of the table)
	int* num_param_vals;   // number of values given for each parameter

	/* need to identify the meaning of each parameter here [index in the array]
     * [param->gam, param->afe, param->lxi, param->ect, param->dens, param->incl] */
    int* param_index;

    float** param_vals;    // array to store the parameter values (as given in the table)

    float******* data_storage;  // storage of a 6-dim table
    void* dat;                  // pointer to point at the start of the table such that it has the correct dimensionality

}xillTable;


typedef struct{
   int nr;
   int ng;

   double* re;
   double* gmin;
   double* gmax;
   double* gstar;
   double* d_gstar;  // bin width for each gstar value

   double*** trff;
   double*** cosne;

   /** the emissivity profile (del and del_inc are only of interest in the LP geometry) **/
	double* emis;       // intensity on the surface of the accretion disc
	double* del_emit;   // angle under which the photon is emitted from the primary source
	double* del_inc;    // angle the photon hits the accretion disk (in the rest frame of the disk)

	double del_ad_risco; // delta of the photon where it would hit the ISCO (irrespective of Rin)
	double del_ad_rmax;  // delta of the photon where it would hit 1000rg (the outer edge of the disk in the relline table)

	int limb_law;

} relSysPar;


typedef struct{
	double refl_frac;
	double refl_frac_norm;
	double f_bh;
	double f_ad;
	double f_inf;
	double f_ad_norm;
} lpReflFrac;

/** angles (cosne) and their distribution over the radial zones **/
typedef struct{
	int n_cosne;
	int n_zones;
	double* cosne;
	double** dist;      // [n_zones][n_cosne]
} rel_cosne;

typedef struct{
	int n_ener;
	int n_zones;
	double* rgrid;      // length=n_zones + 1
	double* ener;       // length=n_ener +1
	double** flux;      // [n_zones][n_ener]
	rel_cosne* rel_cosne;
} rel_spec;


typedef struct{
	double* ener;  // has n_ener+1 elements
	double* incl;  // has n_incl elements
	double** flu;  // [n_incl,n_ener+1]
	int n_ener;
	int n_incl;
}xill_spec;


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


/******************************/
/* define the c_donthcomp function here */
void c_donthcomp(double *ear, int ne, double* param, double *photar);

#endif /* COMMON_H_ */
