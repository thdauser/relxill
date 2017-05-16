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
	int fixReflFrac;
	int model_type;
	int prim_type;
} xillParam;


/** the XILLVER table structure */
typedef struct{

	float* elo;
	float* ehi;
	int n_ener;

	int num_param;

	float* gam;
	float* afe;
	float* lxi;
	float* ect;
	float* incl;

	int n_gam;
	int n_afe;
	int n_lxi;
	int n_ect;
	int n_incl;

	int n_ca_gam;
	int n_ca_afe;
	int n_ca_lxi;
	int n_ca_ect;
	int n_ca_incl;

	float****** dat;

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

} relSysPar;


typedef struct{
	double refl_frac;
	double f_bh;
	double f_ad;
	double f_inf;
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


/******************************/
/* define the c_donthcomp function here */
void c_donthcomp(double *ear, int ne, double* param, double *photar);

#endif /* COMMON_H_ */
