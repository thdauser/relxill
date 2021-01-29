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
#include "relbase.h"

#include "writeOutfiles.h"

// new CACHE routines
cnode *cache_relbase = NULL;

/** caching parameters **/
relParam *cached_rel_param = NULL;
xillParam *cached_xill_param = NULL;

int save_1eV_pos = 0;
// double cached_int_romb_rad = -1.0;

const double ener_xill_norm_lo = 0.1;
const double ener_xill_norm_hi = 1000;
enum {
  n_ener_xill = 3000
};
double *global_ener_xill = NULL;

const int n_ener_std = N_ENER_CONV;
double *global_ener_std = NULL;

specCache *spec_cache = NULL;


static specCache *new_specCache(int n_cache, int n_ener, int *status) {

  specCache *spec = (specCache *) malloc(sizeof(specCache));
  CHECK_MALLOC_RET_STATUS(spec, status, NULL)

  spec->n_cache = n_cache;
  spec->nzones = 0;
  spec->n_ener = n_ener;

  spec->fft_xill = (double ***) malloc(sizeof(double **) * n_cache);
  CHECK_MALLOC_RET_STATUS(spec->fft_xill, status, NULL)

  spec->fft_rel = (double ***) malloc(sizeof(double **) * n_cache);
  CHECK_MALLOC_RET_STATUS(spec->fft_rel, status, NULL)

  spec->xill_spec = (xillSpec **) malloc(sizeof(xillSpec *) * n_cache);
  CHECK_MALLOC_RET_STATUS(spec->xill_spec, status, NULL)

  int ii;
  int jj;
  int m = 2;
  for (ii = 0; ii < n_cache; ii++) {
    spec->fft_xill[ii] = (double **) malloc(sizeof(double *) * m);
    CHECK_MALLOC_RET_STATUS(spec->fft_xill[ii], status, NULL)
    spec->fft_rel[ii] = (double **) malloc(sizeof(double *) * m);
    CHECK_MALLOC_RET_STATUS(spec->fft_rel[ii], status, NULL)

    for (jj = 0; jj < m; jj++) {
      spec->fft_xill[ii][jj] = (double *) malloc(sizeof(double) * n_ener);
      CHECK_MALLOC_RET_STATUS(spec->fft_xill[ii][jj], status, NULL)
      spec->fft_rel[ii][jj] = (double *) malloc(sizeof(double) * n_ener);
      CHECK_MALLOC_RET_STATUS(spec->fft_rel[ii][jj], status, NULL)
    }

    spec->xill_spec[ii] = NULL;

  }

  spec->out_spec = NULL;

  return spec;
}

static void init_specCache(specCache **spec, int *status) {

  if ((*spec) == NULL) {
    (*spec) = new_specCache(N_ZONES_MAX, N_ENER_CONV, status);
  }
}

specCache *init_globalSpecCache(int *status) {
  init_specCache(&spec_cache, status);
  CHECK_RELXILL_ERROR("failed initializing glocal Spec Cache", status);
  return spec_cache;
}

/** convolve the (bin-integrated) spectra f1 and f2 (which need to have a certain binning)
 *  fout: gives the output
 *  f1 input (reflection) specrum
 *  f2 filter
 *  ener has length n+1 and is the energy array
 *  requirements: needs "spec_cache" to be set up
 * **/
static void fft_conv_spectrum(double *ener, const double *fxill, const double *frel, double *fout, int n,
                              int re_rel, int re_xill, int izone, specCache *cache, int *status) {

  long m = 0;
  switch (n) {
    case 512: m = 9;
      break;
    case 1024: m = 10;
      break;
    case 2048: m = 11;
      break;
    case 4096: m = 12;
      break;
    case 8192: m = 13;
      break;
    default: *status = EXIT_FAILURE;
      printf(" *** error: Number of Bins %i not allowed in Convolution!! \n", n);
      break;
  }
  CHECK_STATUS_VOID(*status);

  // needs spec cache to be set up
  assert(cache != NULL);

  /* need to find out where the 1keV for the filter is, which defines if energies are blue or redshifted*/
  if (save_1eV_pos == 0 ||
      (!((ener[save_1eV_pos] <= 1.0) &&
          (ener[save_1eV_pos + 1] > 1.0)))) {
    save_1eV_pos = binary_search(ener, n + 1, 1.0);
  }

  int ii;
  int irot;
  double xcomb[n];
  double ycomb[n];

  /**********************************************************************/
  /** cache either the relat. or the xillver part, as only one of the
   * two changes most of the time (reduce time by 1/3 for convolution) **/
  /**********************************************************************/

  /** #1: for the xillver part **/
  if (re_xill) {
    for (ii = 0; ii < n; ii++) {
      cache->fft_xill[izone][0][ii] = fxill[ii] / (ener[ii + 1] - ener[ii]);
      cache->fft_xill[izone][1][ii] = 0.0;
    }
    FFT_R2CT(1, m, cache->fft_xill[izone][0], cache->fft_xill[izone][1]);
  }
  double *x1 = cache->fft_xill[izone][0];
  double *y1 = cache->fft_xill[izone][1];

  /** #2: for the relat. part **/
  if (re_rel) {
    for (ii = 0; ii < n; ii++) {
      irot = (ii - save_1eV_pos + n) % n;
      cache->fft_rel[izone][0][irot] = frel[ii] / (ener[ii + 1] - ener[ii]);
      cache->fft_rel[izone][1][ii] = 0.0;
    }
    FFT_R2CT(1, m, cache->fft_rel[izone][0], cache->fft_rel[izone][1]);
  }
  double *x2 = cache->fft_rel[izone][0];
  double *y2 = cache->fft_rel[izone][1];

  /* complex multiplication
   * (we need the real part, so we already use the output variable here
   *  to save computing time */
  for (ii = 0; ii < n; ii++) {
    xcomb[ii] = x1[ii] * x2[ii] - y1[ii] * y2[ii];
    ycomb[ii] = y1[ii] * x2[ii] + x1[ii] * y2[ii];
  }

  FFT_R2CT(-1, m, xcomb, ycomb);

  for (ii = 0; ii < n; ii++) {
    fout[ii] = xcomb[ii] * (ener[ii + 1] - ener[ii]);
  }

}

static void print_angle_dist(RelCosne *spec, int izone) {

  FILE *fp = fopen("test_angle_dist.dat", "w+");
  int ii;
  for (ii = 0; ii < spec->n_cosne; ii++) {
    fprintf(fp, " %e %e \n", spec->cosne[ii],
            spec->dist[izone][ii]);
  }
  if (fclose(fp)) {
    exit(1);
  }

}

/**
 * @Function: calcFFTNormFactor
 * @Synopsis: calculate the normalization of the FFT, which is defined to keep the normalization of the
 *           input spectrum and the relat. smearing
 * Take the sum in the given energy band of interested, to avoid problems at the border of the FFT
 * convolution.
 */
double calcFFTNormFactor(const double *ener, const double *fxill, const double *frel, const double *fout, int n) {

  double sum_relline = 0.0;
  double sum_xillver = 0.0;
  double sum_conv = 0.0;
  for (int jj = 0; jj < n; jj++) {
    if (ener[jj] >= EMIN_XILLVER && ener[jj + 1] < EMAX_XILLVER) {
      sum_xillver += fxill[jj];
      sum_relline += frel[jj];
      sum_conv += fout[jj];
    }
  }

  return sum_relline * sum_xillver / sum_conv;
}

void normalizeFFTOutput(const double *ener, const double *fxill, const double *frel, double *fout, int n) {
  double norm_fac = calcFFTNormFactor(ener, fxill, frel, fout, n);

  for (int ii = 0; ii < n; ii++) {
    fout[ii] *= norm_fac;
  }

}
void convolveSpectrumFFTNormalized(double *ener, const double *fxill, const double *frel, double *fout, int n,
                                   int re_rel, int re_xill, int izone, specCache *spec_cache_ptr, int *status) {

  fft_conv_spectrum(ener, fxill, frel, fout, n, re_rel, re_xill, izone, spec_cache_ptr, status);

  normalizeFFTOutput(ener, fxill, frel, fout, n);

}

static void calc_xillver_angdep(double *xill_flux, xillSpec *xill_spec, const double *dist, const int *status) {

  CHECK_STATUS_VOID(*status);

  int ii;
  int jj;
  for (ii = 0; ii < xill_spec->n_ener; ii++) {
    xill_flux[ii] = 0.0;
  }

  for (ii = 0; ii < xill_spec->n_incl; ii++) {
    for (jj = 0; jj < xill_spec->n_ener; jj++) {
      xill_flux[jj] += dist[ii] * xill_spec->flu[ii][jj];
    }
  }

}

static double get_rzone_energ_shift(relParam *param, double rad, double del_emit) {

  // del_emit: any value (for beta!=0 is okay, as it does not matter for the energy shift!)
  return gi_potential_lp(rad, param->a, param->height, param->beta, del_emit);
}

/** check if the complete spectrum is cached and if the energy grid did not change
 *  -> additionally we adapte the energy grid to the new dimensions and energy values
 *     if it's not the case yet
 */
static int is_all_cached(specCache *cache, int n_ener_inp, double *ener_inp, int recompute_xill, int recompute_rel) {

  if ((recompute_xill + recompute_rel) == 0 && (cache->out_spec != NULL)) {
    /** first need to check if the energy grid did change (might happen) **/
    if (cache->out_spec->n_ener != n_ener_inp) {
      return 0;
    }
    int ii;
    // we know that n_ener
    for (ii = 0; ii < n_ener_inp; ii++) {
      if (fabs(cache->out_spec->ener[ii] - ener_inp[ii]) > 1e-4) {
        int jj;
        for (jj = 0; jj < n_ener_inp; jj++) {
          cache->out_spec->ener[jj] = ener_inp[jj];
        }
        return 0;
      }
    }
    return 1;
  } else {
    return 0;
  }

}

static void check_caching_relxill(relParam *rel_param, xillParam *xill_param, int *re_rel, int *re_xill) {


  /** always re-compute if
   *  - the number of zones changed
   *  - we are interested in some output files
   **/
  if ((cached_rel_param != NULL) && (shouldOutfilesBeWritten() == 0)) {

    if (rel_param->num_zones != cached_rel_param->num_zones) {
      if (is_debug_run()) {
        printf("  *** warning :  the number of radial zones was changed from %i to %i \n",
               rel_param->num_zones, cached_rel_param->num_zones);
      }
      *re_rel = 1;
      *re_xill = 1;
      return;
    }
  } else {
    *re_rel = 1;
    *re_xill = 1;
    return;
  }

  /** did any of the relat. parameters change?  **/
  *re_rel = redo_relbase_calc(rel_param, cached_rel_param);

  *re_xill = redo_xillver_calc(rel_param, xill_param, cached_rel_param, cached_xill_param);
}

/** renorm a model (flu) to have the same flux as another model (flu)
 *  (bin-integrated flux, same energy grid!) **/
static void renorm_model(const double *flu0, double *flu, int nbins) {

  double sum_inp = 0.0;
  double sum_out = 0.0;
  int ii;
  for (ii = 0; ii < nbins; ii++) {
    sum_inp += flu0[ii];
    sum_out += flu[ii];
  }
  for (ii = 0; ii < nbins; ii++) {
    flu[ii] *= sum_inp / sum_out;
  }

}

void renorm_xill_spec(float *spec, int n, double lxi, double dens) {
  int ii;
  for (ii = 0; ii < n; ii++) {
    spec[ii] /= pow(10, lxi);
    if (fabs(dens - 15) > 1e-6) {
      spec[ii] /= pow(10, dens - 15);
    }
  }
}

void get_std_relxill_energy_grid(int *n_ener, double **ener, int *status) {
  if (global_ener_std == NULL) {
    global_ener_std = (double *) malloc((N_ENER_CONV + 1) * sizeof(double));
    CHECK_MALLOC_VOID_STATUS(global_ener_std, status)
    get_log_grid(global_ener_std, (N_ENER_CONV + 1), EMIN_RELXILL, EMAX_RELXILL);
  }
  (*n_ener) = N_ENER_CONV;
  (*ener) = global_ener_std;

}

/** BASIC RELCONV FUNCTION : convole any input spectrum with the relbase kernel
 *  (ener has the length n_ener+1)
 *  **/
void relconv_kernel(double *ener_inp, double *spec_inp, int n_ener_inp, relParam *rel_param, int *status) {

  /* get the (fixed!) energy grid for a RELLINE for a convolution
   * -> as we do a simple FFT, we can now take into account that we
   *    need it to be number = 2^N */

  // always do the convolution on this grid
  int n_ener;
  double *ener;
  get_std_relxill_energy_grid(&n_ener, &ener, status);

  rel_spec *rel_profile = relbase(ener, n_ener, rel_param, NULL, status);

  // simple convolution only makes sense for 1 zone !
  assert(rel_profile->n_zones == 1);

  double rebin_flux[n_ener];
  double conv_out[n_ener];
  rebin_spectrum(ener, rebin_flux, n_ener,
                 ener_inp, spec_inp, n_ener_inp);

  spec_cache = init_globalSpecCache(status);
  CHECK_STATUS_VOID(*status);
  fft_conv_spectrum(ener, rebin_flux, rel_profile->flux[0], conv_out, n_ener,
                    1, 1, 0, spec_cache, status);
  CHECK_STATUS_VOID(*status);

  // need to renormalize the convolution? (not that only LP has a physical norm!!)
  if (!do_not_normalize_relline()) {
    renorm_model(rebin_flux, conv_out, n_ener);
  }

  // rebin to the output grid
  rebin_spectrum(ener_inp, spec_inp, n_ener_inp, ener, conv_out, n_ener);

}

/**
 * @Function: get_xillver_angdep_spec
 * @Synopsis: Calculate the Angle Weighted Xillver Spectrum on the Standard Relxill Spectrum
 * @Input:  ener[n_ener]
 *         rel_dist[n_incl]
 *         xill_spec
 * @Output: o_xill_flux  (needs to be allocated, will be overwritten)
 *  [reason for the required allocation is that this will be called in a large
 *   loop and otherwise we would need to allocate a 3000 bin array very often]
 **/
void get_xillver_angdep_spec(double *o_xill_flux,
                             int n_ener,
                             double *ener,
                             double *rel_dist,
                             xillSpec *xill_spec,
                             int *status) {

  double xill_angdist_inp[xill_spec->n_ener];

  calc_xillver_angdep(xill_angdist_inp, xill_spec, rel_dist, status);

  rebin_spectrum(ener, o_xill_flux, n_ener,
                 xill_spec->ener, xill_angdist_inp, xill_spec->n_ener);

}

/*
 * BASIC RELXILL KERNEL FUNCTION : convolve a xillver spectrum with the relbase kernel
 * (ener has the length n_ener+1)
*/
void relxill_kernel(double *ener_inp,
                    double *spec_inp,
                    int n_ener_inp,
                    xillParam *xill_param,
                    relParam *rel_param,
                    int *status) {

  /** only do the calculation once **/
  int n_ener;
  double *ener;
  get_std_relxill_energy_grid(&n_ener, &ener, status);
  assert(ener != NULL);

  xillTable *xill_tab = NULL;
  get_init_xillver_table(&xill_tab, xill_param, status);
  CHECK_STATUS_VOID(*status);

  // in case of an ionization gradient, we need to update the number of zones
  if (is_iongrad_model(rel_param->model_type, xill_param->ion_grad_type)) {
    rel_param->num_zones = get_num_zones(rel_param->model_type, rel_param->emis_type, xill_param->ion_grad_type);
  }

  // make sure the output array is set to 0
  int ii;
  for (ii = 0; ii < n_ener_inp; ii++) {
    spec_inp[ii] = 0.0;
  }

  /*** LOOP OVER THE RADIAL ZONES ***/
  double conv_out[n_ener];
  double single_spec_inp[n_ener_inp];
  for (ii = 0; ii < n_ener; ii++) {
    conv_out[ii] = 0.0;
  }

  /** be careful, as xill_param->ect can get over-written so save the following values**/
  double ecut0 = xill_param->ect;

  /** note that in case of the nthcomp model Ecut is in the frame of the primary source
      but for the bkn_powerlaw it is given in the observer frame */
  double ecut_primary = 0.0;
  if (xill_param->prim_type == PRIM_SPEC_ECUT) {
    ecut_primary = ecut0 * (1 + grav_redshift(rel_param));
  } else if (xill_param->prim_type == PRIM_SPEC_NTHCOMP) {
    ecut_primary = ecut0;
  }

  int recompute_xill = 1;
  int recompute_rel = 1;
  check_caching_relxill(rel_param, xill_param, &recompute_rel, &recompute_xill);

  init_specCache(&spec_cache, status);
  CHECK_STATUS_VOID(*status);

  /** is both already cached we can see if we can simply use the output flux value **/
  if (is_all_cached(spec_cache, n_ener_inp, ener_inp, recompute_xill, recompute_rel)) {
    CHECK_STATUS_VOID(*status);
    for (ii = 0; ii < n_ener_inp; ii++) {
      spec_inp[ii] = spec_cache->out_spec->flux[ii];
    }

    /** if NOT, we need to do a whole lot of COMPUTATIONS **/
  } else {
    CHECK_STATUS_VOID(*status);

    /* *** first, stored the parameters for which we are calculating **/
    set_cached_xill_param(xill_param, &cached_xill_param, status);
    set_cached_rel_param(rel_param, &cached_rel_param, status);

    /* calculate the relline profile **/

    rel_spec *rel_profile = relbase(ener, n_ener, rel_param, xill_tab, status);
    CHECK_STATUS_VOID(*status);

    /* init the xillver spectrum structure **/
    xillSpec *xill_spec_table = NULL;
    double xill_flux[n_ener];


    // currently only working for the LP version (as relxill always has just 1 zone)
    ion_grad *ion = NULL;
    if (is_iongrad_model(rel_param->model_type, xill_param->ion_grad_type)) {

      // make sure the number of zones is correctly set:
      assert(rel_param->num_zones
                 == get_num_zones(rel_param->model_type, rel_param->emis_type, xill_param->ion_grad_type));

      ion = calc_ion_gradient(rel_param, xill_param->lxi, xill_param->ion_grad_index, xill_param->ion_grad_type,
                              rel_profile->rgrid, rel_profile->n_zones, status);
      CHECK_STATUS_VOID(*status);
    }

    for (ii = 0; ii < rel_profile->n_zones; ii++) {
      assert(spec_cache != NULL);

      /** avoid problems where no relxill bin falls into an ionization bin **/
      if (calcSum(rel_profile->flux[ii], rel_profile->n_ener) < 1e-12) {
        continue;
      }

      // now calculate the reflection spectra for each zone (using the angular distribution)
      assert(rel_profile->rel_cosne != NULL);
      if ( shouldOutfilesBeWritten() ) {
        print_angle_dist(rel_profile->rel_cosne, ii);
      }

      // get the energy shift in order to calculate the proper Ecut value (if nzones>1)
      // (the latter part of the IF is a trick to get the same effect as NZONES=1 if during a running
      //  session the number of zones is reset)
      if (rel_profile->n_zones == 1) {
        xill_param->ect = ecut0;
      } else {
        // choose the (linear) middle of the radial zone
        double rzone = 0.5 * (rel_profile->rgrid[ii] + rel_profile->rgrid[ii + 1]);
        double del_emit = 0.0;  // only relevant if beta!=0 (doppler boosting)
        if (ion != NULL) {
          assert(ion->nbins == rel_profile->n_zones);
          del_emit = ion->del_emit[ii];
        }
        xill_param->ect = ecut_primary * get_rzone_energ_shift(rel_param, rzone, del_emit);
      }

      // if we have an ionization gradient defined, we need to set the xlxi to the value of the current zone
      if (ion != NULL) {
        xill_param->lxi = ion->lxi[ii];
      }

      // call the function which calculates the xillver spectrum
      //  - always need to re-compute if we have an ionization gradient, TODO: better caching here
      if (recompute_xill) {
        if (spec_cache->xill_spec[ii] != NULL) {
          free_xill_spec(spec_cache->xill_spec[ii]);
        }
        spec_cache->xill_spec[ii] = get_xillver_spectra(xill_param, status);
      }
      xill_spec_table = spec_cache->xill_spec[ii];

      get_xillver_angdep_spec(xill_flux, n_ener, ener, rel_profile->rel_cosne->dist[ii], xill_spec_table, status);


      // convolve the spectrum **
      //(important for the convolution: need to recompute fft for xillver
      //always if rel changes, as the angular distribution changes !!)
      fft_conv_spectrum(ener, xill_flux, rel_profile->flux[ii], conv_out, n_ener,
                        recompute_rel, 1, ii, spec_cache, status);
      CHECK_STATUS_VOID(*status);

      double normFacFFT = calcFFTNormFactor(ener, xill_flux, rel_profile->flux[ii], conv_out, n_ener);

      // rebin to the output grid
      rebin_spectrum(ener_inp, single_spec_inp, n_ener_inp, ener, conv_out, n_ener);

      // add it to the final output spectrum
      for (int jj = 0; jj < n_ener_inp; jj++) {
        spec_inp[jj] += single_spec_inp[jj] * normFacFFT;
      }

      if (shouldOutfilesBeWritten() && rel_profile->n_zones <= 10) {
        char vstr[200];
        double test_flu[n_ener_inp];
        for (int jj = 0; jj < n_ener_inp; jj++) {
          test_flu[jj] = single_spec_inp[jj] * normFacFFT;
        }
        if (sprintf(vstr, "test_relxill_spec_zones_%03i.dat", ii + 1) == -1) {
          RELXILL_ERROR("failed to get filename", status);
        }
        save_xillver_spectrum(ener_inp, test_flu, n_ener_inp, vstr);
      }

    } /**** END OF LOOP OVER RADIAL ZONES [ii] *****/

    /** important: set the cutoff energy value back to its original value **/
    xill_param->ect = ecut0;

    /** free the ionization parameter structure **/
    free_ion_grad(ion);

    /** initialize the cached output spec array **/
    if ((spec_cache->out_spec != NULL)) {
      if (spec_cache->out_spec->n_ener != n_ener_inp) {
        free_out_spec(spec_cache->out_spec);
        spec_cache->out_spec = init_out_spec(n_ener_inp, ener_inp, status);
        CHECK_STATUS_VOID(*status);
      }
    } else {
      spec_cache->out_spec = init_out_spec(n_ener_inp, ener_inp, status);
      CHECK_STATUS_VOID(*status);
    }

    for (ii = 0; ii < n_ener_inp; ii++) {
      spec_cache->out_spec->flux[ii] = spec_inp[ii];
    }
  } /************* END OF THE HUGE COMPUTATION ********************/

  /** add a primary spectral component and normalize according to the given refl_frac parameter**/
  add_primary_component(ener_inp, n_ener_inp, spec_inp, rel_param, xill_param, status);
}

void set_stdNormXillverEnerygrid(int *status) {
  if (global_ener_xill == NULL) {
    global_ener_xill = (double *) malloc((N_ENER_XILLVER + 1) * sizeof(double));
    CHECK_MALLOC_VOID_STATUS(global_ener_xill, status)
    get_log_grid(global_ener_xill, N_ENER_XILLVER + 1, EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION);
  }
}

EnerGrid *get_stdXillverEnergygrid(int *status) {
  CHECK_STATUS_RET(*status, NULL);

  set_stdNormXillverEnerygrid(status);
  CHECK_STATUS_RET(*status, NULL);

  EnerGrid *egrid = new_EnerGrid(status);
  egrid->ener = global_ener_xill;
  egrid->nbins = N_ENER_XILLVER;

  return egrid;
}

double calcNormWrtXillverTableSpec(const double *flux, const double *ener, const int n, int *status) {
  /* get the normalization of the spectrum with respect to xillver
 *  - everything is normalized using dens=10^15 cm^3
 *  - normalization defined, e.g., in Appendix of Dauser+2016
 *  - needs to be calculated on the specific energy grid (defined globally)
 */

  set_stdNormXillverEnerygrid(status);
  assert(global_ener_xill != NULL);

  double keV2erg = 1.602177e-09;

  // need this to make sure no floating point problems arise
  double floatCompFactor = 1e-6;

  if (ener[n] < (EMAX_XILLVER_NORMALIZATION - floatCompFactor)
      || ener[0] > (EMIN_XILLVER_NORMALIZATION + floatCompFactor)) {
    RELXILL_ERROR("can not calculate the primary spectrum normalization", status);
    printf("  the given energy grid from %e keV to %e keV does not cover the boundaries\n", ener[0], ener[n]);
    printf("  from [%e,%e] necessary for the calcualtion\n", EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION);
    return 0.0;
  }

  double sum_pl = 0.0;
  for (int ii = 0; ii < n; ii++) {
    if (ener[ii] >= EMIN_XILLVER_NORMALIZATION && ener[ii] <= EMAX_XILLVER_NORMALIZATION) {
      sum_pl += flux[ii] * 0.5 * (ener[ii] + ener[ii + 1]) * 1e20 * keV2erg;
    }
  }
  double norm_xillver_table = 1e15 / 4.0 / M_PI;
  return sum_pl / norm_xillver_table;
}

static void printReflectionStrengthInfo(double *ener,
                                        int n_ener,
                                        const double *flu,
                                        relParam *rel_param,
                                        xillParam *xill_param,
                                        const double *pl_flux,
                                        lpReflFrac *struct_refl_frac) {


  // todo: all this to be set by a ENV
  int imin = binary_search(ener, n_ener + 1, RSTRENGTH_EMIN);
  int imax = binary_search(ener, n_ener + 1, RSTRENGTH_EMAX);

  double sum_pl = 0.0;
  double sum = 0.0;
  for (int ii = imin; ii <= imax; ii++) {
    sum_pl += pl_flux[ii];
    sum += flu[ii];
  }

  printf("For a = %.3f, Rin = %.3f, and h = %.2f rg", rel_param->a, rel_param->rin, rel_param->height);
  if (is_iongrad_model(rel_param->model_type, xill_param->ion_grad_type) || rel_param->beta > 1e-6) {
    printf(" and beta=%.3f v/c", rel_param->beta);
  }
  printf(": \n - reflection fraction  %.3f \n - reflection strength is: %.3f \n",
         struct_refl_frac->refl_frac,
         sum / sum_pl);
  printf(" - photons falling into the black hole or plunging region: %.2f%%\n", struct_refl_frac->f_bh * 100);
  printf(" - gravitational redshift from the observer to the primary source is %.3f\n", grav_redshift(rel_param));
}

void calculatePrimarySpectrum(double *pl_flux_xill, double *ener, int n_ener,
                              const relParam *rel_param, const xillParam *xill_param, int *status) {

  CHECK_STATUS_VOID(*status);
  assert(global_ener_xill != NULL);

  if (xill_param->prim_type == PRIM_SPEC_ECUT) {
    /** note that in case of the nthcomp model Ecut is in the frame of the primary source
	    but for the bkn_powerlaw it is given in the observer frame */

    /** IMPORTANT: defintion of Ecut is ALWAYS in the frame of the observer by definition **/
    /**    (in case of the nthcomp primary continuum ect is actually kte ) **/
    double ecut_rest = xill_param->ect;

    for (int ii = 0; ii < n_ener; ii++) {
      pl_flux_xill[ii] = exp(1.0 / ecut_rest) *
          pow(0.5 * (ener[ii] + ener[ii + 1]), -xill_param->gam) *
          exp(-0.5 * (ener[ii] + ener[ii + 1]) / ecut_rest) *
          (ener[ii + 1] - ener[ii]);
    }

  } else if (xill_param->prim_type == PRIM_SPEC_NTHCOMP) {

    double nthcomp_param[5];
    /** important, kTe is given in the primary source frame, so we have to add the redshift here
		 *     however, only if the REL model **/
    double z = 0.0;
    if (rel_param != NULL && rel_param->emis_type == EMIS_TYPE_LP) {
      z = grav_redshift(rel_param);
    }
    get_nthcomp_param(nthcomp_param, xill_param->gam, xill_param->ect, z);
    c_donthcomp(ener, n_ener, nthcomp_param, pl_flux_xill);

  } else if (xill_param->prim_type == PRIM_SPEC_BB) {

    double en;
    for (int ii = 0; ii < n_ener; ii++) {
      en = 0.5 * (ener[ii] + ener[ii + 1]);
      pl_flux_xill[ii] = en * en / (pow(xill_param->kTbb, 4) * (exp(en / xill_param->kTbb) - 1));
      pl_flux_xill[ii] *= (ener[ii + 1] - ener[ii]);
    }
  } else {
    RELXILL_ERROR("trying to add a primary continuum to a model where this does not make sense (should not happen!)",
                  status);
  }
}

void add_primary_component(double *ener, int n_ener, double *flu, relParam *rel_param,
                           xillParam *xill_param, int *status) {

  double pl_flux[n_ener];

  /** need to create a spcific energy grid for the primary component to fulfill the XILLVER NORM condition (Dauser+2016) **/
  EnerGrid *egrid = get_stdXillverEnergygrid(status);
  CHECK_STATUS_VOID(*status);
  double pl_flux_xill[egrid->nbins]; // global energy grid
  calculatePrimarySpectrum(pl_flux_xill, egrid->ener, egrid->nbins, rel_param, xill_param, status);

  double primarySpecNormFactor = 1. / calcNormWrtXillverTableSpec(pl_flux_xill, egrid->ener, egrid->nbins, status);

  /** bin the primary continuum onto the Input grid **/
  rebin_spectrum(ener, pl_flux, n_ener, egrid->ener, pl_flux_xill, egrid->nbins); //TODO: bug, if E<0.1keV in ener grid

  free(egrid);

  for (int ii = 0; ii < n_ener; ii++) {
    pl_flux[ii] *= primarySpecNormFactor;
  }

  /** 2 **  decide if we need to do relat. calculations **/
  if (is_xill_model(xill_param->model_type)) {

    for (int ii = 0; ii < n_ener; ii++) {
      flu[ii] *= fabs(xill_param->refl_frac);
    }
  } else {

    assert(rel_param != NULL);

    // should be cached, as it has been calculated before
    RelSysPar *sysPar = get_system_parameters(rel_param, status);

    lpReflFrac *struct_refl_frac = sysPar->emis->returnFracs;

    if ( xill_param->fixReflFrac > 0 ) {
      /** set the reflection fraction calculated from the height and
       *  spin of the primary source, in this case for the physical
       *  value from Rin to Rout          						 */
      xill_param->refl_frac = struct_refl_frac->refl_frac;

      // special case, if set to "3", it will return only the reflected spectrum
      // with the normalization as predicted
      if (xill_param->fixReflFrac == 3){
        xill_param->refl_frac = - struct_refl_frac->refl_frac;
      }
    }

    /** 4 ** and apply it to primary and reflected spectra **/
    if (rel_param->emis_type == EMIS_TYPE_LP) {
      double g_inf = sqrt(1.0 - (2 * rel_param->height /
          (rel_param->height * rel_param->height + rel_param->a * rel_param->a)));


      /** if the user sets the refl_frac parameter manually, we need to calculate the ratio
       *  to end up with the correct normalization
       */
      double norm_fac_refl = (fabs(xill_param->refl_frac)) / struct_refl_frac->refl_frac;

      double prim_fac = struct_refl_frac->f_inf / 0.5 * pow(g_inf, xill_param->gam);

      for (int ii = 0; ii < n_ener; ii++) {
        pl_flux[ii] *= prim_fac;
        flu[ii] *= norm_fac_refl;
      }
    } else {
      for (int ii = 0; ii < n_ener; ii++) {
        flu[ii] *= fabs(xill_param->refl_frac);
      }
    }

    /** 5 ** if desired, we ouput the reflection fraction and strength (as defined in Dauser+2016) **/
    if ((xill_param->fixReflFrac == 2) && (rel_param->emis_type == EMIS_TYPE_LP)) {
      printReflectionStrengthInfo(ener, n_ener, flu, rel_param, xill_param, pl_flux, struct_refl_frac);
    }

  }

  /** 6 ** add power law component only if desired (i.e., refl_frac > 0)**/
  if (xill_param->refl_frac >= 0) {
    for (int ii = 0; ii < n_ener; ii++) {
      flu[ii] += pl_flux[ii];
    }
  }

}


int did_xill_param_change(xillParam *cpar, xillParam *par) {
  if (comp_single_param_val(par->afe, cpar->afe)) return 1;
  if (comp_single_param_val(par->dens, cpar->dens)) return 1;
  if (comp_single_param_val(par->ect, cpar->ect)) return 1;
  if (comp_single_param_val(par->gam, cpar->gam)) return 1;
  if (comp_single_param_val(par->lxi, cpar->lxi)) return 1;
  if (comp_single_param_val(par->kTbb, cpar->kTbb)) return 1;
  if (comp_single_param_val(par->frac_pl_bb, cpar->frac_pl_bb)) return 1;
  if (comp_single_param_val(par->z, cpar->z)) return 1;

  if (comp_single_param_val((double) par->prim_type, (double) cpar->prim_type)) return 1;
  if (comp_single_param_val((double) par->model_type, (double) cpar->model_type)) return 1;

  if (comp_single_param_val(par->ion_grad_index, cpar->ion_grad_index)) return 1;
  if (comp_single_param_val((double) par->ion_grad_type, (double) cpar->ion_grad_type)) return 1;

  return 0;
}

/* check if values, which need a re-computation of the relline profile, have changed */
int redo_xillver_calc(relParam *rel_param, xillParam *xill_param, relParam *ca_rel_param, xillParam *ca_xill_param) {

  int redo = 1;

  if ((ca_rel_param != NULL) && (ca_xill_param != NULL)) {

    redo = did_xill_param_change(ca_xill_param, xill_param);

    /** did spin or h change (means xillver needs to be re-computed as well, due to Ecut) **/
    if (comp_single_param_val(rel_param->a, ca_rel_param->a) ||
        comp_single_param_val(rel_param->height, ca_rel_param->height)) {
      redo = 1;
    }

  }

  return redo;
}

int redo_relbase_calc(relParam *rel_param, relParam *ca_rel_param) {

  if (did_rel_param_change(ca_rel_param, rel_param)) {
    return 1;
  } else {
    return 0;
  }

}


/* the relbase function calculating the basic relativistic line shape for a given parameter setup
 * (assuming a 1keV line, by a grid given in keV!)
 * input: ener(n_ener), param
 * optinal input: xillver grid
 * output: photar(n_ener)     */
rel_spec *relbase_multizone(double *ener,
                            const int n_ener,
                            relParam *param,
                            xillTable *xill_tab,
                            double *radialZones,
                            int nzones,
                            int *status) {

  CHECK_STATUS_RET(*status, NULL);

  inpar *inp = set_input(ener, n_ener, param, NULL, status);

  // check caching here and also re-set the cached parameter values
  cache_info *ca_info = cli_check_cache(cache_relbase, inp, check_cache_relpar, status);

  // set a pointer to the spectrum
  rel_spec *spec = NULL;

  // initialize parameter values (has an internal cache)
  RelSysPar *sysPar = get_system_parameters(param, status);
  assert(sysPar != NULL);

  if (is_relbase_cached(ca_info) == 0) {

    // init the spectra where we store the flux
    param->num_zones = nzones;
    init_rel_spec(&spec, param, xill_tab, radialZones, &ener, n_ener, status);

    // calculate line profile (returned units are 'cts/bin')
    relline_profile(spec, sysPar, status);

    // normalize it and calculate the angular distribution (if necessary)
    renorm_relline_profile(spec, param, status);

    // last step: store parameters and cached rel_spec (this prepends a new node to the cache)
    set_cache_relbase(&cache_relbase, param, spec, status);
    if (is_debug_run() && *status == EXIT_SUCCESS) {
      printf(" DEBUG:  Adding new RELBASE eval to cache; the count is %i \n", cli_count_elements(cache_relbase));
    }
  } else {
    if (is_debug_run()) {
      printf(" DEBUG:  RELBASE-Cache: re-using calculated values\n");
    }
    spec = ca_info->store->data->relbase_spec;
  }

  if (shouldOutfilesBeWritten()) {
    save_emis_profiles(sysPar);
    save_relline_profile(spec);
  }


  // free the input structure
  free(inp);
  free(ca_info);

  // CHECK_RELXILL_DEFAULT_ERROR(status);

  return spec;
}


rel_spec *relbase(double *ener, const int n_ener, relParam *param, xillTable *xill_tab, int *status) {

  double *rgrid = get_rzone_grid(param->rin, param->rout, param->num_zones, param->height, status);

  return relbase_multizone(ener, n_ener, param, xill_tab, rgrid, param->num_zones, status);
}

void free_rel_cosne(RelCosne *spec) {
  if (spec != NULL) {
    //	free(spec->ener);  we do not need this, as only a pointer for ener is assigned
    free(spec->cosne);
    if (spec->dist != NULL) {
      int ii;
      for (ii = 0; ii < spec->n_zones; ii++) {
        free(spec->dist[ii]);
      }
    }
    free(spec->dist);
    free(spec);
  }
}

void free_rel_spec(rel_spec *spec) {
  if (spec != NULL) {
    free(spec->ener);
    free(spec->rgrid);
    if (spec->flux != NULL) {
      int ii;
      for (ii = 0; ii < spec->n_zones; ii++) {
        if (spec->flux[ii] != NULL) {
          free(spec->flux[ii]);
        }
      }
    }
    free(spec->flux);
    if (spec->rel_cosne != NULL) {
      free_rel_cosne(spec->rel_cosne);
    }
    free(spec);
  }
}

void free_cached_tables(void) {

  free_relprofile_cache();

  free_cached_relTable();
  free_cached_lpTable();
  free_cached_xillTable();

  free(cached_rel_param);
  free(cached_xill_param);

  free_specCache();

  free(global_ener_std);
  free(global_ener_xill);

}

void free_fft_cache(double ***sp, int n1, int n2) {

  int ii;
  int jj;
  if (sp != NULL) {
    for (ii = 0; ii < n1; ii++) {
      if (sp[ii] != NULL) {
        for (jj = 0; jj < n2; jj++) {
          free(sp[ii][jj]);
        }
      }
      free(sp[ii]);
    }
    free(sp);
  }

}

OutSpec *init_out_spec(int n_ener, const double *ener, int *status) {

  OutSpec *spec = (OutSpec *) malloc(sizeof(OutSpec));
  CHECK_MALLOC_RET_STATUS(spec, status, NULL)

  spec->n_ener = n_ener;
  spec->ener = (double *) malloc(sizeof(double) * n_ener);
  CHECK_MALLOC_RET_STATUS(spec->ener, status, NULL)
  spec->flux = (double *) malloc(sizeof(double) * n_ener);
  CHECK_MALLOC_RET_STATUS(spec->flux, status, NULL)

  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    spec->ener[ii] = ener[ii];
    spec->flux[ii] = 0.0;
  }

  return spec;
}

void free_out_spec(OutSpec *spec) {
  if (spec != NULL) {
    free(spec->ener);
    free(spec->flux);
    free(spec);
  }
}

void free_specCache(void) {

  int ii;
  int m = 2;
  if (spec_cache != NULL) {
    if (spec_cache->xill_spec != NULL) {
      for (ii = 0; ii < spec_cache->n_cache; ii++) {
        if (spec_cache->xill_spec[ii] != NULL) {
          free_xill_spec(spec_cache->xill_spec[ii]);
        }
      }
      free(spec_cache->xill_spec);
    }

    if (spec_cache->fft_xill != NULL) {
      free_fft_cache(spec_cache->fft_xill, spec_cache->n_cache, m);
    }

    if (spec_cache->fft_rel != NULL) {
      free_fft_cache(spec_cache->fft_rel, spec_cache->n_cache, m);
    }

    free_out_spec(spec_cache->out_spec);

  }

  free(spec_cache);

}

/** free the CLI cache **/

void free_cache(void) {
  free_cache_syspar();
  cli_delete_list(&cache_relbase);
}
