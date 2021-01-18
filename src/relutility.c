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

/** linear interpolation in 1 dimension **/
double interp_lin_1d(double ifac_r, double rlo, double rhi) {
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

/** calculate the gravitational redshift **/
double grav_redshift(const relParam *param) {
  if (param->emis_type == EMIS_TYPE_LP) {
    return 1.0 / sqrt(1.0 - 2 * param->height /
        (param->height * param->height + param->a * param->a)) - 1.0;
  } else {
    // important: without a geometrical assumption no grav. redshift can be calculated
    return 0.0;
  }
}

double relat_abberation(double del, double beta) {
  return acos((cos(del) - beta) / (1 - beta * cos(del)));
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
    int debug = atof(env);
    if (debug == 1) {
      return 1;
    }
  }
  return 0;
}

/** check if we are currently debugging the model **/
int shouldOutfilesBeWritten(void) {
  char *env;
  env = getenv("RELXILL_WRITE_OUTFILES");
  if (env != NULL) {
    int debug = atof(env);
    if (debug == 1) {
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

  double delta_eq = h * h - 2 * h + a * a;
  double q2 = (pow(sin(del), 2)) * (pow((h * h + a * a), 2) / delta_eq) - a * a;

  double beta_fac = sqrt(pow((h * h + a * a), 2) - delta_eq * (q2 + a * a));
  beta_fac = gam * (1.0 + sign * beta_fac / (h * h + a * 2) * bet);

  return gi / beta_fac;
}

/** print the xillver spectrum   **/
void save_xillver_spectrum(double *ener, double *flu, int n_ener, char *fname) {

  FILE *fp = fopen(fname, "w+");
  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    fprintf(fp, " %e \t %e \t %e \n", ener[ii], ener[ii + 1], flu[ii]);
  }
  if (fclose(fp)) exit(1);
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

void rebin_spectrum(double *ener, double *flu, int nbins, double *ener0, double *flu0, int nbins0) {

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

/*** we calculate the disk density from  Shakura & Sunyaev (1973)
 *    - for zone A as describe in their publication,  formula (Eq 2.11)
 *    - only the radial dependence is picked up here  (viscosity alpha=const.)
 *    - normalized such that dens(rms) = 1
 *                                                               ***/
double density_ss73_zone_a(double radius, double rms) {
  return pow((radius / rms), (3. / 2)) * pow((1 - sqrt(rms / radius)), -2);
}

// calculate the log(xi) for given density and emissivity
static double cal_lxi(double dens, double emis) {
  return log10(4.0 * M_PI * emis / dens);
}

// determine the radius of maximal ionization
static double cal_lxi_max_ss73(double *re, double *emis, int nr, double rin, int *status) {

  CHECK_STATUS_RET(*status, 0.0);

  double rad_max_lxi = pow((11. / 9.), 2)
      * rin;  // we use the same definition as Adam with r_peak = (11/9)^2 rin to be consistent (does not matter much)

  // radial AD grid is sorted descending (!)
  int kk = inv_binary_search(re, nr, rad_max_lxi);
  double interp = (rad_max_lxi - re[kk + 1]) / (re[kk] - re[kk + 1]);

  double emis_max_lxi = interp_lin_1d(interp, emis[kk + 1], emis[kk]);

  double lxi_max = cal_lxi(density_ss73_zone_a(rad_max_lxi, rin), emis_max_lxi);

  return lxi_max;
}

static void save_ion_profile(ion_grad *ion) {

  FILE *fp = fopen("test_ion_grad_relxill.dat", "w+");
  int ii;
  for (ii = 0; ii < ion->nbins; ii++) {
    fprintf(fp, " %e \t %e \t %e \n", ion->r[ii], ion->r[ii + 1], ion->lxi[ii]);
  }
  if (fclose(fp)) exit(1);

}

/** *** set log(xi) to obey the limits of the xillver table: TODO: check if we need to adjust the normalization as well  ***
 *  NOTE: with correctly set xpsec/isis limits, it is only possible to reach the lower boundary       **/
static void lxi_set_to_xillver_bounds(double *pt_lxi) {
  /**  TODO: Need to define this globally **/
  double xlxi_tab_min = 0.0;
  double xlxi_tab_max = 4.7;

  // #1: set the value of xi to the lowest value of the table
  if (*pt_lxi < xlxi_tab_min) {
    *pt_lxi = xlxi_tab_min;
  } else if (*pt_lxi > xlxi_tab_max) {
    //	#2: high ionization: we approximately assume such a highly ionized disk acts as a mirror
    *pt_lxi = xlxi_tab_max;
  }

}

ion_grad *calc_ion_gradient(relParam *rel_param,
                            double xlxi0,
                            double xindex,
                            int type,
                            double *rgrid,
                            int n,
                            int *status) {

  CHECK_STATUS_RET(*status, NULL);

  ion_grad *ion = new_ion_grad(rgrid, n, status);
  CHECK_STATUS_RET(*status, NULL);

  double rmean[n];
  double del_inc[n];
  int ii;
  for (ii = 0; ii < n; ii++) {
    rmean[ii] = 0.5 * (rgrid[ii] + rgrid[ii + 1]);
  }

  if (type == ION_GRAD_TYPE_PL) {
    for (ii = 0; ii < n; ii++) {
      ion->lxi[ii] = (exp(xlxi0))
          * pow((rmean[ii] / rmean[0]), -1.0 * xindex);  // TODO: check if we need to subtract xlxi_tab_min here
      ion->lxi[ii] = log(ion->lxi[ii]);

      lxi_set_to_xillver_bounds(&(ion->lxi[ii]));

    }

  } else if (type == ION_GRAD_TYPE_ALPHA) {
    double dens[n];
    double rin = rgrid[0];

    // TODO: use a better approach to not linearly interpolate but rather average over the profile?
    double emis_zones[n];

    // we need the emissivity profile (should be cached, so no extra effort required here)
    RelSysPar *sysPar = get_system_parameters(rel_param, status);
    emisProfile *emis_profile = sysPar->emis;

    assert(emis_profile->del_inc != NULL);
    inv_rebin_mean(emis_profile->re, emis_profile->emis, sysPar->nr, rmean, emis_zones, n, status);
    inv_rebin_mean(emis_profile->re, emis_profile->del_inc, sysPar->nr, rmean, del_inc, n, status);
    inv_rebin_mean(emis_profile->re, emis_profile->del_emit, sysPar->nr, rmean, ion->del_emit, n, status);

    // calculate the maximal ionization assuming r^-3 and SS73 alpha disk
    double lxi_max = cal_lxi_max_ss73(emis_profile->re, emis_profile->emis, emis_profile->nr, rin, status);

    // the maximal ionization is given as input parameter, so we need to normalize our calculation by this value
    double fac_lxi_norm = xlxi0 - lxi_max; // subtraction instead of division because of the log

    /** calculate the density for a  stress-free inner boundary condition, i.e., R0=rin in SS73)  **/
    for (ii = 0; ii < n; ii++) {
      dens[ii] = density_ss73_zone_a(rmean[ii], rin);

      // now we can use the emissivity to calculate the ionization
      ion->lxi[ii] = cal_lxi(dens[ii], emis_zones[ii]) + fac_lxi_norm;

      ion->lxi[ii] += log10(cos(M_PI / 4) / cos(del_inc[ii]));

      lxi_set_to_xillver_bounds(&(ion->lxi[ii]));
    }

  } else if (type
      == ION_GRAD_TYPE_CONST) {  // should not happen, as this will be approximated by 1 zone (but just in case we get here...)
    for (ii = 0; ii < n; ii++) {
      ion->lxi[ii] = xlxi0;
    }
  } else {
    printf(" *** ionization type with number %i not implemented \n", type);
    printf("     choose either %i for the PL, %i for the ALPHA-disk, or %i for constant\n",
           ION_GRAD_TYPE_PL, ION_GRAD_TYPE_ALPHA, ION_GRAD_TYPE_CONST);
    RELXILL_ERROR("unknown ionization gradient type", status);
  }

  if (is_debug_run()) {
    save_ion_profile(ion);
  }

  if (*status != EXIT_SUCCESS) {
    RELXILL_ERROR("calculating the ionization gradient failed due to previous error", status);
  }

  return ion;
}

ion_grad *new_ion_grad(double *r, int n, int *status) {

  ion_grad *ion = (ion_grad *) malloc(sizeof(ion_grad));
  CHECK_MALLOC_RET_STATUS(ion, status, NULL)

  ion->r = (double *) malloc((n + 1) * sizeof(double));
  CHECK_MALLOC_RET_STATUS(ion->r, status, NULL)
  ion->lxi = (double *) malloc((n) * sizeof(double));
  CHECK_MALLOC_RET_STATUS(ion->lxi, status, NULL)
  ion->fx = (double *) malloc((n) * sizeof(double));
  CHECK_MALLOC_RET_STATUS(ion->fx, status, NULL)
  ion->del_emit = (double *) malloc((n) * sizeof(double));
  CHECK_MALLOC_RET_STATUS(ion->del_emit, status, NULL)

  ion->nbins = n;

  int ii;
  for (ii = 0; ii < n; ii++) {
    ion->r[ii] = r[ii];
    ion->lxi[ii] = 0.0;
    ion->fx[ii] = 0.0;
    ion->del_emit[ii] = M_PI / 4.; // assume default 45 deg (xillver assumption), only used if beta>0
  }
  // radius goes to n+1
  ion->r[n] = r[n];

  return ion;
}

void free_ion_grad(ion_grad *ion) {

  if (ion != NULL) {
    if (ion->r != NULL) {
      free(ion->r);
    }
    if (ion->lxi != NULL) {
      free(ion->lxi);
    }
    if (ion->fx != NULL) {
      free(ion->fx);
    }
    free(ion);
  }
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
  for (ii = 0; ii < n0 - 2; ii++) {  // only go to the second to last bin

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

double calc_g_inf(double height, double a) {
  return sqrt(1.0 - (2 * height /
      (height * height + a * a)));
}

void setArrayToZero(double *arr, int n) {
  for (int jj = 0; jj < n; jj++) {
    arr[jj] = 0.0;
  }
}


void multiplyArray(double *arr, int n, double factor) {
  for (int jj = 0; jj < n; jj++) {
    arr[jj] *= factor;
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
