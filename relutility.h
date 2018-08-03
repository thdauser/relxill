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

    Copyright 2018 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELUTILITY_H_
#define RELUTILITY_H_

#include "relbase.h"
#include "xilltable.h"


/****** DEFINE FUNCTION DEFINITIONS ******/

#define RELXILL_ERROR(msg,status) (relxill_error(__func__, msg,status))

#define CHECK_RELXILL_ERROR(msg,status) (check_relxill_error(__func__, msg,status))

#define CHECK_STATUS_RET(status, retval) \
  if (EXIT_SUCCESS!=status) return(retval);

#define CHECK_STATUS_VOID(status) \
  if (EXIT_SUCCESS!=status) return;

#define CHECK_STATUS_BREAK(status) \
  if (EXIT_SUCCESS!=status) break;

#define CHECK_MALLOC_VOID_STATUS(a,status) \
	if (NULL==a) { \
		RELXILL_ERROR("memory allocation failed",status); \
		return;\
	}

#define CHECK_MALLOC_RET_STATUS(a,status,retval) \
	if (NULL==a) { \
		RELXILL_ERROR("memory allocation failed",status); \
		return retval;\
	}

#define CHECK_STATUS_RET(status, retval) \
  if (EXIT_SUCCESS!=status) return(retval);


/**************************/
/** Function Definitions **/
/**************************/

/** linear interpolation in 1 dimension **/
double interp_lin_1d(double ifac_r, double rlo, double rhi);

/** log interpolation in 1 dimension **/
double interp_log_1d(double ifac_r, double rlo, double rhi);

/** linear interpolation in 2 dimensions **/
double interp_lin_2d(double ifac1, double r1lo, double r1hi, double ifac2, double r2lo, double r2hi);
double interp_lin_2d_float(double ifac1, double ifac2, float r11, float r12, float r21, float r22);

/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_log_grid(double* ener, int n_ener, double emin, double emax);

/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_lin_grid(double* ener, int n_ener, double emin, double emax);

/* get the current version number */
void get_version_number(char** vstr, int* status);

/* print relxill error message */
void relxill_error(const char* const func, const char* const msg, int* status);

/* check and print relxill error message */
void check_relxill_error(const char* const func, const char* const msg, int* status);

/* check and print relxill error message */
void check_relxill_error(const char* const func, const char* const msg, int* status);

/* inverse binary search */
int inv_binary_search(double* arr,int n,double val);

/* inverse binary search */
int inv_binary_search_float(float* arr,int n,float val);

/* binary search */
int binary_search_float(float* arr,int n,float val);

int binary_search(double* arr,int n,double val);

/** calculate the reflection fraction **/
lpReflFrac* calc_refl_frac(relSysPar* sysPar, relParam* param, int* status);

/** trapez integration around a single bin (returns only r*dr*PI!) **/
double trapez_integ_single(double* re, int ii, int nr);

/* calculate the radius of marginal stability */
double kerr_rms(double a);

/* get the rplus value (size if the black hole event horizon */
double kerr_rplus(double a);

/** test if it is a relxill flavour model **/
int is_relxill_model(int model_type);

/** check if we are currently debugging the model **/
int is_debug_run( void );

/** get a radial grid on the accretion disk in order to calculate a relline for each zone **/
void get_rzone_grid(double rmin, double rmax, double* rgrid, int nzones, double h, int* status);

void get_rgrid(double* ener, int n_ener, double emin, double emax);

/** convert gstar to energy */
double gstar2ener(double g, double gmin, double gmax, double ener);

/** calculate the doppler factor for a moving primary source **/
double doppler_factor(double del, double bet);

/** calculates g = E/E_i in the lamp post geometry (see, e.g., 27 in Dauser et al., 2013, MNRAS) **/
double gi_potential_lp(double r, double a, double h, double bet, double del);

/** print the xillver spectrum   **/
void save_xillver_spectrum(double* ener, double* flu, int n_ener, char* fname);

/* A simple implementation of the FFT taken
   from http://paulbourke.net/miscellaneous/dft/
   (uses the Radix-2 Cooley-Tukey algorithm) */
void FFT_R2CT(short int dir,long m,double *x,double *y);

/** rebin spectrum to a given energy grid length of ener is nbins+1       **/
void rebin_spectrum(double* ener, double* flu, int nbins, double* ener0, double* flu0, int nbins0);

/** calculate the gravitational redshift **/
double grav_redshift(relParam* param);

/** get the relxill table path (dynamically from env variable)  **/
char* get_relxill_table_path( void );

/** get the number of zones **/
int get_num_zones(int model_type, int emis_type);

/** check if it is the xillver model **/
int is_xill_model(int model_type);

void get_nthcomp_param( double* nthcomp_param, double gam, double kte, double z);

int do_renorm_model(relParam* rel_param);

/** check if we should return the relline/relconv physical norm from ENV **/
int do_not_normalize_relline( void );

#endif /* RELUTILITY_H_ */
