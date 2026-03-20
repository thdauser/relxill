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
#ifndef RELUTILITY_H_
#define RELUTILITY_H_

#include "common.h"
#include "envoptions.h"


/**************************/
/** Function Definitions **/
/**************************/


/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_log_grid(double *ener, int n_ener, double emin, double emax);

/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_lin_grid(double *ener, int n_ener, double emin, double emax);

/* inverse binary search */
int inv_binary_search(const double *arr, int n, double val);

/* inverse binary search */
int inv_binary_search_float(const float *arr, int n, float val);

/* binary search */
int binary_search_float(const float *arr, int n, float val);

int binary_search(const double *arr, int n, double val);

/** trapez integration around a single bin (returns only r*dr*PI!) **/
double trapez_integ_single(const double *re, int ii, int nr);

double trapez_integ_single_rad_ascending(const double *re, int ii, int nr);

/** test if it is a relxill flavor model **/
int is_relxill_model(int model_type);

void getLogGrid(double *ener, int n_ener, double emin, double emax);

/** convert gstar to energy */
double gstar2ener(double g, double gmin, double gmax, double ener);

/** get the number of zones **/
int get_num_zones(int model_type, int emis_type, int ion_grad_type);

void get_nthcomp_param(double *nthcomp_param, double gam, double kte, double z);

int do_renorm_model(relParam * rel_param);

// check for the model type
int is_iongrad_model(int ion_grad_type);

int is_ns_model(int model_type);

int is_co_model(int model_type);

int is_xill_model(int model_type);

int is_alpha_model(int model_type);

int is_returnrad_bb_model(int model_type);

int is_lp_primary_source(int emis_type);

int is_ring_primary_source(int emis_type);

int is_slab_primary_source(int emis_type);

int is_any_primary_source(int emis_type);


double calcSum(const double *array, int n_array);

double calcSumInEnergyBand(const double *array, int n_array, const double *ener, double valLo, double valHi);

void setArrayToZero(double *arr, int n);

EnerGrid *new_EnerGrid(int *status);

TestSpectrum *new_Spectrum(int *status);

void free_Spectrum(TestSpectrum * *spec);

TestSpectrum *getNewSpec(double emin, double emax, int nbins, int *status);

void invertArray(double *vals, int n);

double get_ipol_factor_radius(double rlo, double rhi, double del_inci, double radius);

void get_ipol_factor(const float value, const float *arr, const int n_arr, int *ind, double *ifac);

void get_ipol_factor_double(const double value, const double *arr, const int n_arr, int *ind, double *ifac);

void get_fine_radial_grid(double rin, double rout, double *re, int nr);

void print_version_number(void);

void fits_movenam_hdu_cpp(fitsfile *fptr, int type, const char *extname, int extver, int *status);

double convert_to_isco(double rin_relxill, double a);

#endif /* RELUTILITY_H_ */