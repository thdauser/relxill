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

#include "relutility.h"

#include "relphysics.h"
#include "relbase.h"
#include "writeOutfiles.h"

/** linear interpolation in 1 dimension **/
double interp_lin_1d(double ifac_r, double rlo, double rhi) {
  return ifac_r * rhi + (1.0 - ifac_r) * rlo;
}

/** linear interpolation in 1 dimension **/
double interp_lin_1d_float(double ifac_r, float rlo, float rhi) {
  return ifac_r * rhi + (1.0 - ifac_r) * rlo;
}


double interp_log_1d(double ifac_r, double rlo, double rhi) {
  return exp(ifac_r * log(rhi) + (1.0 - ifac_r) * log(rlo));
}

/** linear interpolation in 2 dimensions **/
double interp_lin_2d(double ifac1, double ifac2, double r11, double r12, double r21, double r22) {
  return (1.0 - ifac1) * (1.0 - ifac2) * r11 +
      (ifac1) * (1.0 - ifac2) * r12 +
      (1.0 - ifac1) * (ifac2) * r21 +
      (ifac1) * (ifac2) * r22;
}

double interp_lin_2d_float(double ifac1, double ifac2, float r11, float r12, float r21, float r22) {
  return (1.0 - ifac1) * (1.0 - ifac2) * r11 +
      (ifac1) * (1.0 - ifac2) * r12 +
      (1.0 - ifac1) * (ifac2) * r21 +
      (ifac1) * (ifac2) * r22;
}

void relxill_error(const char *const func, const char *const msg, int *status) {
  *status = EXIT_FAILURE;
  printf(" *** error in relxill (%s): %s!\n", func, msg);
}

void relxill_warning(const char *const msg) {
  printf(" *** warning from relxill: %s!\n", msg);
}

int is_xill_model(int model_type) {
  if ((model_type == MOD_TYPE_XILLVERDENS) || (model_type == MOD_TYPE_XILLVER)
      || (model_type == MOD_TYPE_XILLVERNS) || (model_type == MOD_TYPE_XILLVERCO)
      || (model_type == MOD_TYPE_XILLVERDENS_NTHCOMP)
      || (model_type == MOD_TYPE_XILLVER_NTHCOMP)) {
    return 1;
  } else {
    return 0;
  }
}

// ion-gradient model, which is not set to constant ionization
//  - in case ion_grad_type=constant, it is working as a normal model
int is_iongrad_model(int model_type, int ion_grad_type) {
  if ((model_type == MOD_TYPE_RELXILLLPION) && (ion_grad_type != ION_GRAD_TYPE_CONST)) {
    return 1;
  } else {
    return 0;
  }
}

/** check and report FITS error   */
void relxill_check_fits_error(const int *status) {
  if (*status != EXIT_SUCCESS) {
    char errtext[30];
    fits_get_errstatus(*status, errtext);
    printf("   cfitsio error: %s \n", errtext);
  }
}


void check_relxill_error(const char *const func, const char *const msg, int *status) {
  if (*status != EXIT_SUCCESS) {
    *status = EXIT_FAILURE;
    printf(" *** error in relxill (%s): %s!\n", func, msg);
  }
}

void print_relxill_test_msg(const char *const func, const char *const msg) {
  printf(" TEST: %s %s", func, msg);
}

void print_relxill_test_result(int status) {
  if (status == EXIT_SUCCESS) {
    printf("   --> OK\n");
  } else {
    printf("   --> ERROR\n");
    CHECK_RELXILL_DEFAULT_ERROR(&status);
  }
}

/**  FLOAT search for value "val" in array "arr" (sorted ASCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
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

int is_returnrad_model(int model_type) {
  if (model_type == MOD_TYPE_RELXILLBBRET) {
    return 1;
  } else {
    return 0;
  }
}

/** trapez integration around a single bin
 *  caveat: only returns half of the full 2*PI*r*dr due to computational speed**/
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

/** convert gstar to energy */
double gstar2ener(double g, double gmin, double gmax, double ener) {
  return (g * (gmax - gmin) + gmin) * ener;
}

/** get a radial grid on the accretion disk in order to calculate a relline for each zone **/
double *get_rzone_grid(double rmin, double rmax, int nzones, double h, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  double *rgrid = malloc(sizeof(double) * (nzones + 1));
  CHECK_MALLOC_RET_STATUS(rgrid, status, rgrid)

  if (nzones == 1) {
    rgrid[0] = rmin;
    rgrid[1] = rmax;
  } else {

    double r_transition = rmin;
    int indr = 0;

    // if h > rmin we choose a log grid for r<h
    if (h > rmin) {

      r_transition = h;

      get_log_grid(rgrid, nzones + 1, rmin, rmax);
      indr = binary_search(rgrid, nzones + 1, r_transition);

      r_transition = rgrid[indr];

    }

    if (indr < nzones) {

      double rlo = r_transition;
      double rhi = rmax; // rgrid[nzones];
      // add 1/r for larger radii
      int ii;
      for (ii = indr; ii < nzones + 1; ii++) {
        rgrid[ii] = 1.0 * (ii - indr) / (nzones - indr) * (1.0 / rhi - 1.0 / rlo) + 1.0 / rlo;
        rgrid[ii] = fabs(1.0 / rgrid[ii]);
      }

    }

  }

  return rgrid;
}

/** get the relxill table path (dynamically from env variable)  **/
char *get_relxill_table_path(void) {
  char *path;
  path = getenv("RELXILL_TABLE_PATH");
  if (path != NULL) {
    return path;
  } else {
    return RELXILL_TABLE_PATH;
  }
}

/** check if we are currently debugging the model **/
int is_debug_run(void) {
  char *env;
  env = getenv("DEBUG_RELXILL");
  if (env != NULL) {
    int debug = (int) atof(env);
    if (debug == 1) {
      return 1;
    }
  }
  return 0;
}


/** check if we are currently debugging the model **/
int shouldAuxInfoGetPrinted(void) {

  char *env = getenv("RELXILL_OUTPUT_DETAILS");
  if (env != NULL) {
    int envval = (int) strtod(env, NULL);
    if (envval == 1) {
      return 1;
    }
  }
  return 0;
}


/** check if we are currently debugging the model **/
int shouldOutfilesBeWritten(void) {
  char *env;
  env = getenv("RELXILL_OUTPUT_FILES");
  if (env != NULL) {
    int envval = (int) strtod(env, NULL);
    if (envval == 1) {
      return 1;
    }
  }
  return 0;
}


/** check if we should return the relline/relconv physical norm from ENV **/
int do_not_normalize_relline(void) {
  char *env;
  env = getenv("RELLINE_PHYSICAL_NORM");
  if (env != NULL) {
    int phys_norm = atof(env);
    if (phys_norm == 1) {
      return 1;
    }
  }
  return 0;
}

/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_log_grid(double *ener, int n_ener, double emin, double emax) {
  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    ener[ii] = 1.0 * ii / (n_ener - 1) * (log(emax) - log(emin)) + log(emin);
    ener[ii] = exp(ener[ii]);
  }
}

/* get a logarithmic grid from emin to emax with n_ener bins  */
void getLogGrid(double *ener, int n_ener, double emin, double emax) {
  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    ener[ii] = 1.0 * ii / (n_ener - 1) * (1.0 / emax - 1.0 / emin) + 1.0 / emin;
    ener[ii] = fabs(1.0 / ener[ii]);
  }
}

/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_lin_grid(double *ener, int n_ener, double emin, double emax) {
  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    ener[ii] = 1.0 * ii / (n_ener - 1) * (emax - emin) + emin;
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

void get_ipol_factor(const float value, const float* arr, const int n_arr, int *ind, double *ifac) {
  (*ind) = binary_search_float(arr, n_arr, (float) value);
  (*ifac) = (value - arr[*ind]) /
      (arr[*ind + 1] - arr[*ind]);
}



/* get RMS (ISCO) for the Kerr Case */
double kerr_rms(double a) {
  //	 accounts for negative spin
  double sign = 1.0;
  if (a < 0) {
    sign = -1.0;
  }

  double Z1 = 1.0 + pow(1.0 - a * a, 1.0 / 3.0) * (pow(1.0 + a, 1.0 / 3.0) + pow(1.0 - a, 1.0 / 3.0));
  double Z2 = sqrt((3.0 * a * a) + (Z1 * Z1));

  return 3.0 + Z2 - sign * sqrt((3.0 - Z1) * (3.0 + Z1 + (2 * Z2)));
}

/* get the rplus value (size if the black hole event horizon */
double kerr_rplus(double a) {
  return 1 + sqrt(1 - a * a);
}

/** calculate the doppler factor for a moving primary source **/
double doppler_factor(double del, double bet) {
  return sqrt(1.0 - bet * bet) / (1.0 + bet * cos(del));
}

/** calculates g = E/E_i in the lamp post geometry (see, e.g., 27 in Dauser et al., 2013, MNRAS) **/
double gi_potential_lp(double r, double a, double h, double bet, double del) {

  /** ! calculates g = E/E_i in the lamp post geometry
    ! (see, e.g., page 48, Diploma Thesis, Thomas Dauser) **/
  double ut_d = ((r * sqrt(r) + a) / (sqrt(r) * sqrt(r * r - 3 * r + 2 * a * sqrt(r))));
  double ut_h = sqrt((h * h + a * a) / (h * h - 2 * h + a * a));

  double gi = ut_d / ut_h;

  // check if we need to calculate the additional factor for the velocity
  if (fabs(bet) < 1e-6) {
    return gi;
  }

  double gam = 1.0 / sqrt(1.0 - bet * bet);

  // get the sign for the equation
  double sign = 1.0;
  if (del > M_PI / 2) {
    sign = -1.0;
  }

// ** Adam, 6.7.2021: **
// So given the expression for Carter’s q:
//q^2 = \sin^2\delta (h^2+a^2)^2/\Delta_h - a^2,
//it is clear that
//q^2 + a^2 = \sin^2\delta (h^2+a^2)^2/\Delta_h
//Therefore:
//(h^2+a^2)^2 - \Delta_h (q^2+a^2) = (h^2+a^2)^2 - \sin^2\delta (h^2+a^2)^2
//= (h^2+a^2)^2 \cos^2\delta
//
//Therefore:
//\sqrt{ (h^2+a^2)^2 - \Delta_h (q^2+a^2) } / (h^2+a^2) = \cos\delta
//
//Therefore your equation (27) becomes:
//glp = glp(beta=0) / { \gamma [ 1 -/+ \beta\cos\delta ] }

  double delta_eq = h * h - 2 * h + a * a;
  double q2 = (pow(sin(del), 2)) * (pow((h * h + a * a), 2) / delta_eq) - a * a;

  double beta_fac = sqrt(pow((h * h + a * a), 2) - delta_eq * (q2 + a * a));
  beta_fac = gam * (1.0 + sign * beta_fac / (h * h + a * a) * bet);

  return gi / beta_fac;
}

/* A simple implementation of the FFT taken from http://paulbourke.net/miscellaneous/dft/
   (uses the Radix-2 Cooley-Tukey algorithm)

   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform
*/
void FFT_R2CT(short int dir, long m, double *x, double *y) {

  long n, i, i1, j, k, i2, l, l1, l2;
  double c1, c2, tx, ty, t1, t2, u1, u2, z;

  /* Calculate the number of points */
  n = 1;
  for (i = 0; i < m; i++) {
    n *= 2;
  }

  /* Do the bit reversal */
  i2 = n >> 1;
  j = 0;
  for (i = 0; i < n - 1; i++) {
    if (i < j) {
      tx = x[i];
      ty = y[i];
      x[i] = x[j];
      y[i] = y[j];
      x[j] = tx;
      y[j] = ty;
    }
    k = i2;
    while (k <= j) {
      j -= k;
      k >>= 1;
    }
    j += k;
  }

  /* Compute the FFT */
  c1 = -1.0;
  c2 = 0.0;
  l2 = 1;
  for (l = 0; l < m; l++) {
    l1 = l2;
    l2 <<= 1;
    u1 = 1.0;
    u2 = 0.0;
    for (j = 0; j < l1; j++) {
      for (i = j; i < n; i += l2) {
        i1 = i + l1;
        t1 = u1 * x[i1] - u2 * y[i1];
        t2 = u1 * y[i1] + u2 * x[i1];
        x[i1] = x[i] - t1;
        y[i1] = y[i] - t2;
        x[i] += t1;
        y[i] += t2;
      }
      z = u1 * c1 - u2 * c2;
      u2 = u1 * c2 + u2 * c1;
      u1 = z;
    }
    c2 = sqrt((1.0 - c1) / 2.0);
    if (dir == 1)
      c2 = -c2;
    c1 = sqrt((1.0 + c1) / 2.0);
  }

  /* Scaling for forward transform */
  if (dir == 1) {
    for (i = 0; i < n; i++) {
      x[i] /= n;
      y[i] /= n;
    }
  }

}

/** get the number of zones on which we calculate the relline-spectrum **/
int get_num_zones(int model_type, int emis_type, int ion_grad_type) {

  char *env;
  env = getenv("RELXILL_NUM_RZONES");
  int env_n_zones = 0;
  if (env != NULL) {
    env_n_zones = atof(env);
  }


  // set the number of zones in radial direction (1 for relline/conv model, N_ZONES for xill models)
  if (is_iongrad_model(model_type, ion_grad_type)) {
    if (env != NULL) {
      if ((env_n_zones > 9) && (env_n_zones <= N_ZONES_MAX)) {
        return env_n_zones;
      } else {
        printf(" *** warning: value of %i for RELXILL_NUM_RZONES not within required interval of [10,%i] \n",
               env_n_zones,
               N_ZONES_MAX);
      }
    }
    return N_ZONES_MAX;
  } else if (is_relxill_model(model_type) && (emis_type == EMIS_TYPE_LP)) {

    if (env != NULL) {
      if ((env_n_zones > 0) && (env_n_zones <= N_ZONES_MAX)) {
        return env_n_zones;
      } else {
        printf(" *** warning: value of %i for RELXILL_NUM_RZONES not within required interval of [1,%i] \n",
               env_n_zones,
               N_ZONES_MAX);
      }
    }
    return N_ZONES;
  } else {
    return 1;
  }

}

/** rebin spectrum to a given energy grid
 *  length of ener is nbins+1       **/

void rebin_spectrum(const double *ener, double *flu, int nbins, const double *ener0, const double *flu0, int nbins0) {

  int ii;
  int jj;
  int imin = 0;
  int imax = 0;

  for (ii = 0; ii < nbins; ii++) {

    flu[ii] = 0.0;

    /* check of the bin is outside the given energy range */
    if ((ener0[0] <= ener[ii + 1]) && (ener0[nbins0] >= ener[ii])) {

      /* need to make sure we are in the correct bin */
      while (ener0[imin] <= ener[ii] && imin <= nbins0) {
        imin++;
      }
      // need to set it back, as we just crossed to the next bin
      if (imin > 0) {
        imin--;
      }
      while ((ener0[imax] <= ener[ii + 1] && imax < nbins0)) {
        imax++;
      }
      if (imax > 0) {
        imax--;
      }

      double elo = ener[ii];
      double ehi = ener[ii + 1];
      if (elo < ener0[imin]) elo = ener0[imin];
      if (ehi > ener0[imax + 1]) ehi = ener0[imax + 1];

      if (imax == imin) {
        flu[ii] = (ehi - elo) / (ener0[imin + 1] - ener0[imin]) * flu0[imin];
      } else {

        double dmin = (ener0[imin + 1] - elo) / (ener0[imin + 1] - ener0[imin]);
        double dmax = (ehi - ener0[imax]) / (ener0[imax + 1] - ener0[imax]);

        flu[ii] += flu0[imin] * dmin + flu0[imax] * dmax;

        for (jj = imin + 1; jj <= imax - 1; jj++) {
          flu[ii] += flu0[jj];
        }

      }

    }

  }
}

int do_renorm_model(relParam *rel_param) {

  int renorm;

  if (is_relxill_model(rel_param->model_type)) {
    if (rel_param->emis_type == EMIS_TYPE_LP || is_returnrad_model(rel_param->model_type)
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
  nthcomp_param[3] = 1.0;  // inp_type
  nthcomp_param[4] = z;

}


// for x0 descending and xn ascending, calculate the mean at xn
void inv_rebin_mean(double *x0, double *y0, int n0, double *xn, double *yn, int nn, int *status) {
  if (xn[0] > xn[nn - 1] || x0[nn - 1] > x0[0]) {
    RELXILL_ERROR(" *** grid in wrong order", status);
    return;
  }
  if (xn[0] < x0[n0 - 1] || xn[nn - 1] > x0[0]) {
    RELXILL_ERROR(" *** new grid is larger than the input grid", status);
    return;
  }

  int in = nn - 1; // traverse new array backwards
  int ii;
  for (ii = 0; ii < n0 - 1; ii++) {  // only go to the second to last bin

    if (x0[ii] > xn[in] && x0[ii + 1] <= xn[in]) {

      double ifac_r = (xn[in] - x0[ii + 1]) / (x0[ii] - x0[ii + 1]);
      yn[in] = interp_lin_1d(ifac_r, y0[ii + 1], y0[ii]);

      in--;
      if (in < 0) { // we can stop if once we are below zero as the array is filled
        break;
      }
    }

  }

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


double calc_g_inf(double height, double a) {
  return sqrt(1.0 - (2 * height /
      (height * height + a * a)));
}

void setArrayToZero(double *arr, int n) {
  for (int jj = 0; jj < n; jj++) {
    arr[jj] = 0.0;
  }
}

// for x0 and xn ascending, calculate the mean at xn  (x0 has n0+1 bins and xn has nn+1 bins)
void rebin_mean_flux(double *xn, double *yn, int nn, double *x0, double *y0, int n0, int *status) {
  if (xn[0] > xn[nn] || x0[0] > x0[n0]) {
    RELXILL_ERROR(" *** grid in wrong order", status);
    return;
  }
  /* if (xn[0]<x0[n0-1] || xn[nn-1]>x0[0]){
     RELXILL_ERROR(" *** new grid is larger than the input grid",status);
     return;
   } */

  int ii = 1; // we start at one, as we always need one bin lower than the xn[0] bin
  for (int in = 0; in < nn; in++) {  // only go to the second to last bin

    yn[in] = 0.0; // set it zero by default

    double xn_m = 0.5 * (xn[in] + xn[in + 1]);
    while (0.5 * (x0[ii - 1] + x0[ii]) < xn_m) {
      ii++;
      if (ii >= n0) {
        break;
      }
    }
    ii--;

    if (ii == 0) {  // at the lowest bin, we need to extrapolate
      ii = 1;
    }

    double x0_m_lo = 0.5 * (x0[ii - 1] + x0[ii]);
    double x0_m_hi = 0.5 * (x0[ii] + x0[ii + 1]);

    if (xn_m > x0_m_lo && xn_m <= x0_m_hi) {
      double ifac_r = (xn_m - x0_m_lo) / (x0_m_hi - x0_m_lo);
      yn[in] = interp_lin_1d(ifac_r, y0[ii - 1], y0[ii]);
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
double calcSumInEnergyBand(const double *array, int n_array, double *ener, double valLo, double valHi) {
  double testSum = 0.0;
  for (int jj = 0; jj < n_array; jj++) {
    if (ener[jj] >= valLo && ener[jj + 1] <= valHi) {
      testSum += array[jj];
    }
  }
  return testSum;
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



EnerGrid *new_EnerGrid(int *status) {

  EnerGrid *egrid = malloc(sizeof(EnerGrid));
  CHECK_MALLOC_RET_STATUS(egrid, status, egrid)
  return egrid;
}

Spectrum *new_Spectrum(int *status) {

  Spectrum *spec = malloc(sizeof(Spectrum));
  CHECK_MALLOC_RET_STATUS(spec, status, spec)
  return spec;
}

void free_Spectrum(Spectrum **spec) {

  if (*spec != NULL) {
    free((*spec)->ener);
    free((*spec)->flux);
    free(*spec);
  }

}

Spectrum *getNewSpec(double emin, double emax, int nbins, int *status) {

  double *ener = malloc(sizeof(double) * (nbins + 1));
  double *flux = malloc(sizeof(double) * nbins);

  Spectrum *spec = new_Spectrum(status);

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
