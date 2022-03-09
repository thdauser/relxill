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
#ifndef RELUTILITY_H_
#define RELUTILITY_H_

#include "common.h"


/****** DEFINE FUNCTION DEFINITIONS ******/

#define RELXILL_ERROR(msg, status) (relxill_error(__func__, msg,status))

#define CHECK_RELXILL_ERROR(msg, status) (check_relxill_error(__func__, msg,status))

#define CHECK_RELXILL_DEFAULT_ERROR(status) (check_relxill_error(__func__, "function evaluation failed",status))

#define CHECK_STATUS_RET(status, retval) \
 if (EXIT_SUCCESS!=status) return(retval)

#define CHECK_STATUS_VOID(status)  \
 if (EXIT_SUCCESS!=status) return

#define CHECK_STATUS_BREAK(status) \
 if (EXIT_SUCCESS!=status) break

#define CHECK_MALLOC_VOID_STATUS(a, status) \
    if (NULL==a) { \
        RELXILL_ERROR("memory allocation failed",status); \
        return;\
    }

#define CHECK_MALLOC_RET_STATUS(a, status, retval) \
    if (NULL==a) { \
        RELXILL_ERROR("memory allocation failed",status); \
        return retval;\
    }

#define PRINT_RELXILL_TEST_MSG_DEFAULT() (print_relxill_test_msg(__func__,""))
#define PRINT_RELXILL_TEST_MSG(msg) (print_relxill_test_msg(__func__,msg))

/**************************/
/** Function Definitions **/
/**************************/

/** linear interpolation in 1 dimension **/
double interp_lin_1d(double ifac_r, double rlo, double rhi);
double interp_lin_1d_float(double ifac_r, float rlo, float rhi);

/** log interpolation in 1 dimension **/
double interp_log_1d(double ifac_r, double rlo, double rhi);

/** linear interpolation in 2 dimensions **/
double interp_lin_2d(double ifac1, double ifac2, double r11, double r12, double r21, double r22);
double interp_lin_2d_float(double ifac1, double ifac2, float r11, float r12, float r21, float r22);

/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_log_grid(double *ener, int n_ener, double emin, double emax);

/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_lin_grid(double *ener, int n_ener, double emin, double emax);

/* print relxill error message */
void relxill_error(const char *const func, const char *const msg, int *status);

/* print a standardized warning message */
void relxill_warning(const char *const msg);

/* check and print relxill error message */
void check_relxill_error(const char *const func, const char *const msg, int *status);

/** check and report FITS error   */
void relxill_check_fits_error(const int *status);

void print_relxill_test_msg(const char *const func, const char *const msg);

void print_relxill_test_result(int status);

/* inverse binary search */
int inv_binary_search(const double *arr, int n, double val);

/* inverse binary search */
int inv_binary_search_float(const float *arr, int n, float val);

/* binary search */
int binary_search_float(const float *arr, int n, float val);

int binary_search(const double *arr, int n, double val);

double relat_abberation(double del, double beta);

/** trapez integration around a single bin (returns only r*dr*PI!) **/
double trapez_integ_single(const double *re, int ii, int nr);
double trapez_integ_single_rad_ascending(const double *re, int ii, int nr);

/* calculate the radius of marginal stability */
double kerr_rms(double a);

/* get the rplus value (size if the black hole event horizon */
double kerr_rplus(double a);

/** test if it is a relxill flavor model **/
int is_relxill_model(int model_type);

/** check if we are currently debugging the model **/
int is_debug_run(void);

/** get a radial grid on the accretion disk in order to calculate a relline for each zone **/
double *get_rzone_grid(double rmin, double rmax, int nzones, double h, int *status);

void getLogGrid(double *ener, int n_ener, double emin, double emax);

/** convert gstar to energy */
double gstar2ener(double g, double gmin, double gmax, double ener);

/** calculate the doppler factor for a moving primary source **/
double doppler_factor(double del, double bet);

/** calculates g = E/E_i in the lamp post geometry (see, e.g., 27 in Dauser et al., 2013, MNRAS) **/
double gi_potential_lp(double r, double a, double h, double bet, double del);


/* A simple implementation of the FFT taken
   from http://paulbourke.net/miscellaneous/dft/
   (uses the Radix-2 Cooley-Tukey algorithm) */
void FFT_R2CT(short int dir, long m, double *x, double *y);

/** rebin spectrum to a given energy grid length of ener is nbins+1       **/
void rebin_spectrum(const double *ener, double *flu, int nbins, const double *ener0, const double *flu0, int nbins0);

/** calculate the gravitational redshift **/
double grav_redshift(const relParam *param);

/** get the relxill table path (dynamically from env variable)  **/
char *get_relxill_table_path(void);

/** get the number of zones **/
int get_num_zones(int model_type, int emis_type, int ion_grad_type);

/** check if it is the xillver model **/
int is_xill_model(int model_type);

void get_nthcomp_param(double *nthcomp_param, double gam, double kte, double z);

int do_renorm_model(relParam *rel_param);

/** check if we should return the relline/relconv physical norm from ENV **/
int do_not_normalize_relline(void);

/** is is a model for which we want to calculate the ionization gradient? **/
int is_iongrad_model(int ion_grad_type);

/** for x0 descending and xn ascending, calculate the mean at xn from y0 **/
void inv_rebin_mean(double *x0, double *y0, int n0, double *xn, double *yn, int nn, int *status);

void rebin_mean_flux(double *x0, double *y0, int n0, double *xn, double *yn, int nn, int *status);

// energy shift from the primary source to the observer
double calc_g_inf(double height, double a);

double calcSum(const double *array, int n_array);
double calcSumInEnergyBand(const double *array, int n_array, double *ener, double valLo, double valHi);

void normSpec(double *spec, int n_ener);

void setArrayToZero(double *arr, int n);
void multiplyArray(double *arr, int n, double factor);

EnerGrid *new_EnerGrid(int *status);

Spectrum *new_Spectrum(int *status);
void free_Spectrum(Spectrum **spec);
Spectrum *getNewSpec(double emin, double emax, int nbins, int *status);

int shouldOutfilesBeWritten(void);

void invertArray(double *vals, int n);

double get_ipol_factor_radius(double rlo, double rhi, double del_inci, double radius);

void get_ipol_factor(const float value, const float* arr, const int n_arr, int *ind, double *ifac);

void get_fine_radial_grid(double rin, double rout, double *re, int nr);

int shouldAuxInfoGetPrinted(void);

void print_version_number(void);

#endif /* RELUTILITY_H_ */
