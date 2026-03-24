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

#include "relutility.h"
#include "common.h"
#include "config.h"

int version_number_printed = 0;

void print_version_number(void) {
  if (version_number_printed == 0) {
    printf(" *** loading RELXILL model (version %s) *** \n", PROJECT_VER);
    version_number_printed = 1;
  }
}


int is_xill_model(int model_type) {
  if ((model_type == MOD_TYPE_XILLVERDENS) || (model_type == MOD_TYPE_XILLVER)
      || (model_type == MOD_TYPE_XILLVERNS) || (model_type == MOD_TYPE_XILLVERCO)
      || (model_type == MOD_TYPE_XILLVER_NTHCOMP)) {
    return 1;
  } else {
    return 0;
  }
}

int is_iongrad_model(int ion_grad_type) {
  if (ion_grad_type != ION_GRAD_TYPE_CONST) {
    return 1;
  } else {
    return 0;
  }
}

int is_alpha_model(int model_type) {
  if (model_type == MOD_TYPE_RELXILLLPALPHA) {
    return 1;
  } else {
    return 0;
  }
}

int is_ns_model(int model_type) {
  if ((model_type == MOD_TYPE_RELXILLNS) || (model_type == MOD_TYPE_XILLVERNS)) {
    return 1;
  } else {
    return 0;
  }
}

int is_co_model(int model_type) {
  if ((model_type == MOD_TYPE_RELXILLCO) || (model_type == MOD_TYPE_XILLVERCO)) {
    return 1;
  } else {
    return 0;
  }
}

int is_lp_primary_source(const int emis_type) {
  if (emis_type == EMIS_TYPE_LP) {
    return 1;
  }
  return 0;
}

int is_ring_primary_source(const int emis_type) {
  if (emis_type == EMIS_TYPE_RING) {
    return 1;
  }
  return 0;
}

int is_slab_primary_source(const int emis_type) {
  if (emis_type == EMIS_TYPE_SLAB) {
    return 1;
  }
  return 0;
}

int is_any_primary_source(const int emis_type) {
  if (is_lp_primary_source(emis_type) || is_ring_primary_source(emis_type) || is_slab_primary_source(emis_type)) {
    return 1;
  }
  return 0;
}

/** check and report FITS error   */


/** @brief FLOAT search for value "val" in array "arr" (sorted ASCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] *

 	 @description return -1 upon error (if array too small), otherwise always the index will always
 	 be within the given length of the array, meaning 0 also for val<arr[0] or n-1 for val>arr[n-1]
 	 */
int binary_search_float(const float *arr, int n, float val) {
  if (n <= 1) return -1;

  int klo = 0;
  int khi = n - 1;
  int k;
  while ((khi - klo) > 1) {
    k = (khi + klo) / 2;
    if (arr[k] > val) {
      khi = k;
    } else {
      klo = k;
    }
  }
  return klo;
}

/**  search for value "val" in array "arr" (sorted ASCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int binary_search(const double *arr, int n, double val) {
  if (n <= 1) return -1;

  int klo = 0;
  int khi = n - 1;
  int k;
  while ((khi - klo) > 1) {
    k = (khi + klo) / 2;
    if (arr[k] > val) {
      khi = k;
    } else {
      klo = k;
    }
  }
  return klo;
}

/**  FLOAT search for value "val" in array "arr" (sorted DESCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int inv_binary_search_float(const float *arr, int n, float val) {
  if (n <= 1) return -1;

  int klo = 0;
  int khi = n - 1;
  int k;
  while ((khi - klo) > 1) {
    k = (khi + klo) / 2;
    if (arr[k] < val) {
      khi = k;
    } else {
      klo = k;
    }
  }
  return klo;
}

/**  search for value "val" in array "arr" (sorted DESCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int inv_binary_search(const double *arr, int n, double val) {
  if (n <= 1) return -1;

  int klo = 0;
  int khi = n - 1;
  int k;
  while ((khi - klo) > 1) {
    k = (khi + klo) / 2;
    if (arr[k] < val) {
      khi = k;
    } else {
      klo = k;
    }
  }
  return klo;
}

/** test if it is a relxill flavour model **/
int is_relxill_model(int model_type) {
  if (model_type < 0) {
    return 1;
  } else {
    return 0;
  }
}

int is_returnrad_bb_model(int model_type) {
  if (model_type == MOD_TYPE_RELXILLBBRET || model_type == MOD_TYPE_RELXILLBBRETXILL) {
    return 1;
  }
  return 0;
}

/** trapez integration around a single bin
 *  caveat: only returns half of the full 2*PI*r*dr due to computational speed
 *  radius is descending **/
double trapez_integ_single(const double *re, int ii, int nr) {
  double dr;
  // dr is defined such that the full disk is covered once, with NO overlapping bins
  if (ii == 0) {
    dr = 0.5 * (re[ii] - re[ii + 1]);
  } else if (ii == nr - 1) {
    dr = 0.5 * (re[ii - 1] - re[ii]);
  } else {
    dr = 0.5 * (re[ii - 1] - re[ii + 1]);
  }
  return re[ii] * dr * M_PI;
}

double trapez_integ_single_rad_ascending(const double *re, int ii, int nr) {
  double dr;
  // dr is defined such that the full disk is covered once, with NO overlapping bins
  if (ii == 0) {
    dr = 0.5 * (re[ii + 1] - re[ii]);
  } else if (ii == nr - 1) {
    dr = 0.5 * (re[ii] - re[ii - 1]);
  } else {
    dr = 0.5 * (re[ii + 1] - re[ii - 1]);
  }
  return re[ii] * dr * M_PI;
}

/** convert gstar to energy */
double gstar2ener(double g, double gmin, double gmax, double ener) {
  return (g * (gmax - gmin) + gmin) * ener;
}


/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_log_grid(double *ener, int n_ener, double emin, double emax) {
  int ii;
  assert(ener != NULL);
  for (ii = 0; ii < n_ener; ii++) {
    ener[ii] = 1.0 * ii / (n_ener - 1) * (log(emax) - log(emin)) + log(emin);
    ener[ii] = exp(ener[ii]);
  }
}

double get_ipol_factor_radius(double rlo, double rhi, double del_inci, double radius) {
  double inter_r;
  // for larger angles logarithmic interpolation works slightly better
  if (del_inci / M_PI * 180.0 <= 75.0) {
    inter_r = (radius - rlo) / (rhi - rlo);
  } else {
    inter_r = (log(radius) - log(rlo)) /
              (log(rhi) - log(rlo));
  }
  return inter_r;
}

void get_ipol_factor(const float value, const float *arr, const int n_arr, int *ind, double *ifac) {
  (*ind) = binary_search_float(arr, n_arr, (float) value);
  (*ifac) = (value - arr[*ind]) /
            (arr[*ind + 1] - arr[*ind]);
}

void get_ipol_factor_double(const double value, const double *arr, const int n_arr, int *ind, double *ifac) {
  (*ind) = binary_search(arr, n_arr, value);
  (*ifac) = (value - arr[*ind]) /
            (arr[*ind + 1] - arr[*ind]);
}


/** get the number of zones on which we calculate the relline-spectrum **/
int get_num_zones(int model_type, int emis_type, int ion_grad_type) {
  int env_n_zones = get_from_env_relline_num_zones();

  // set the number of zones in radial direction (1 for relline/conv model, N_ZONES for m_cache_xill models)
  if (is_iongrad_model(ion_grad_type)) {
    if (env_n_zones > 0) {
      if ((env_n_zones > 9) && (env_n_zones <= N_ZONES_MAX)) {
        return env_n_zones;
      } else {
        printf(" *** warning: value of %i for RELXILL_NUM_RZONES not within required interval of [10,%i] \n",
               env_n_zones,
               N_ZONES_MAX);
      }
    }
    return N_ZONES_IONGRAD;
  } else if (is_relxill_model(model_type) && is_any_primary_source(emis_type)) {
    if (env_n_zones > 0) {
      if (env_n_zones <= N_ZONES_MAX) {
        return env_n_zones;
      }
      printf(" *** warning: value of %i for RELXILL_NUM_RZONES not within required interval of [1,%i] \n",
             env_n_zones,
             N_ZONES_MAX);
    }
    return N_ZONES;
  } else {
    return 1;
  }
}

/** rebin spectrum to a given energy grid
 *  length of ener is num_flux_bins+1       **/
int do_renorm_model(relParam *rel_param) {
  int renorm;

  if (is_relxill_model(rel_param->model_type)) {
    if (is_any_primary_source(rel_param->emis_type) || is_returnrad_bb_model(rel_param->model_type)
        || do_not_normalize_relline()) {
      renorm = 0;
    } else {
      renorm = 1;
    }
  } else {
    if (do_not_normalize_relline()) {
      renorm = 0;
    } else {
      renorm = 1;
    }
  }

  return renorm;
}

void get_nthcomp_param(double *nthcomp_param, double gam, double kte, double z) {
  nthcomp_param[0] = gam;
  nthcomp_param[1] = kte;
  nthcomp_param[2] = 0.05; // ktbb
  nthcomp_param[3] = 1.0; // inp_type
  nthcomp_param[4] = z;
}


/*  get the fine radial grid */
void get_fine_radial_grid(double rin, double rout, double *re, int nr) {
  double r1 = 1.0 / sqrt(rout);
  double r2 = 1.0 / sqrt(rin);
  int ii;
  for (ii = 0; ii < nr; ii++) {
    re[ii] = ((double) (ii)) * (r2 - r1) / (nr - 1) + r1;
    re[ii] = pow(1.0 / (re[ii]), 2);
    assert(re[ii] > 1.0);
  }
}


void setArrayToZero(double *arr, int n) {
  if (arr != NULL) {
    for (int jj = 0; jj < n; jj++) {
      arr[jj] = 0.0;
    }
  }
}

/* Function: calculate the sum of an array
 *
 */
double calcSum(const double *array, int n_array) {
  double testSum = 0.0;
  for (int jj = 0; jj < n_array; jj++) {
    testSum += array[jj];
  }
  return testSum;
}

/*
 * ener has n_array+1 bins
 */
double calcSumInEnergyBand(const double *array, int n_array, const double *ener, double valLo, double valHi) {
  double testSum = 0.0;
  for (int jj = 0; jj < n_array; jj++) {
    if (ener[jj] >= valLo && ener[jj + 1] <= valHi) {
      testSum += array[jj];
    }
  }
  return testSum;
}

/*
 * @brief calculate the energy flux (keV/cm^2/s) for an Xspec spectrum (photons/bin)
 * @note ener has n_ener+1 bins
 */
double calc_ener_flux_from_spectrum(const double *spec_binned, int n_ener, const double *ener) {
  double ener_flu = 0.0;
  for (int jj = 0; jj < n_ener; jj++) {
    ener_flu += spec_binned[jj] * 0.5 * (ener[jj] + ener[jj + 1]);
  }
  return ener_flu;
}


void normSpec(double *spec, int n_ener) {
  double sumSpec = 0.0;

  for (int ii = 0; ii < n_ener; ii++) {
    sumSpec += spec[ii];
  }
  for (int ii = 0; ii < n_ener; ii++) {
    spec[ii] /= sumSpec;
  }
}


TestSpectrum *new_Spectrum(int *status) {
  TestSpectrum *spec = malloc(sizeof(TestSpectrum));
  CHECK_MALLOC_RET_STATUS(spec, status, spec)
  return spec;
}

void free_Spectrum(TestSpectrum **spec) {
  if (*spec != NULL) {
    free((*spec)->ener);
    free((*spec)->flux);
    free(*spec);
  }
}

TestSpectrum *getNewSpec(double emin, double emax, int nbins, int *status) {
  double *ener = malloc(sizeof(double) * (nbins + 1));
  double *flux = malloc(sizeof(double) * nbins);

  TestSpectrum *spec = new_Spectrum(status);

  get_log_grid(ener, nbins + 1, emin, emax);

  setArrayToZero(flux, nbins);

  spec->nbins = nbins;
  spec->ener = ener;
  spec->flux = flux;

  return spec;
}

void invertArray(double *vals, int n) {
  double storage[n];
  for (int ii = 0; ii < n; ii++) {
    storage[n - ii - 1] = vals[ii];
  }

  for (int ii = 0; ii < n; ii++) {
    vals[ii] = storage[ii];
  }
}


void fits_movenam_hdu_cpp(fitsfile *fptr, int type, const char *extname, int extver, int *status) {
  fits_movnam_hdu(fptr, type, (char *) extname, extver, status);
}

double convert_to_isco(double rin_relxill, double a) {
  return 0;
}