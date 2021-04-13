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
#include "fftw/fftw3.h"   // assumes installation in heasoft

#include "writeOutfiles.h"

// new CACHE routines
cnode *cache_relbase = NULL;


int save_1eV_pos = 0;

double *global_ener_xill = NULL;

const int n_ener_std = N_ENER_CONV;
double *global_ener_std = NULL;

specCache *global_spec_cache = NULL;


static specCache *new_specCache(int n_cache, int n_ener, int *status) {

  specCache *spec = (specCache *) malloc(sizeof(specCache));
  CHECK_MALLOC_RET_STATUS(spec, status, NULL)

  spec->n_cache = n_cache;
  spec->nzones = 0;
  spec->n_ener = n_ener;

  spec->conversion_factor_energyflux = NULL;

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

static void init_specCache(specCache **spec, const int n_zones, int *status) {
  if ((*spec) == NULL) {
    (*spec) = new_specCache(n_zones, N_ENER_CONV, status);
  }
}


specCache *init_global_specCache(int *status) {
  init_specCache(&global_spec_cache, N_ZONES_MAX, status);
  CHECK_RELXILL_ERROR("failed initializing Relconv Spec Cache", status);
  return global_spec_cache;
}

static double* calculate_energyflux_conversion(const double* ener, int n_ener, int* status){

  double* factor = (double *) malloc(sizeof(double) * n_ener );
  CHECK_MALLOC_RET_STATUS(factor, status, NULL)

  for(int ii=0; ii<n_ener; ii++){
    factor[ii] = 0.5*(ener[ii]+ener[ii+1]) / (ener[ii+1] - ener[ii]);
  }

  return factor;
}

/** convolve the (bin-integrated) spectra f1 and f2 (which need to have a certain binning)
 *  fout: gives the output
 *  f1 input (reflection) specrum
 *  f2 filter
 *  ener has length n+1 and is the energy array
 *  requirements: needs "spec_cache" to be set up
 * **/
void fft_conv_spectrum(double *ener, const double *fxill, const double *frel, double *fout, int n,
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

  if (cache->conversion_factor_energyflux == NULL){
    cache->conversion_factor_energyflux = calculate_energyflux_conversion(ener, n, status);
  }


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
      cache->fft_xill[izone][0][ii] = fxill[ii] * cache->conversion_factor_energyflux[ii] ;
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
      cache->fft_rel[izone][0][irot] = frel[ii] * cache->conversion_factor_energyflux[ii]; /// (ener[ii + 1] - ener[ii]) * ener[ii];
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
    fout[ii] = xcomb[ii] /  cache->conversion_factor_energyflux[ii]; //* (ener[ii + 1] - ener[ii]) / ener[ii];
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
  for (int ii = 0; ii < n; ii++) {
    spec[ii] /=  pow(10, lxi);  // do not cast to float (fails refdata)
    if (fabs(dens - 15) > 1e-6) {
      spec[ii] /=  pow(10, dens - 15); // do not cast to float (fails refdata)
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

  specCache* spec_cache = init_global_specCache(status);
  CHECK_STATUS_VOID(*status);
  convolveSpectrumFFTNormalized(ener, rebin_flux, rel_profile->flux[0], conv_out, n_ener,
                    1, 1, 0, spec_cache, status);
  CHECK_STATUS_VOID(*status);

  // rebin to the output grid
  rebin_spectrum(ener_inp, spec_inp, n_ener_inp, ener, conv_out, n_ener);

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
 * output: photar(n_ener)  [photons/bin]   */
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
  CHECK_STATUS_RET(*status, NULL);
  assert(sysPar != NULL);

  if (is_relbase_cached(ca_info) == 0) {

    // init the spectra where we store the flux
    param->num_zones = nzones;
    init_rel_spec(&spec, param, xill_tab, radialZones, &ener, n_ener, status);

    // calculate line profile (returned units are 'photons/bin')
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
    free(radialZones); // necessary is it is allocated outside this function TODO: change it!
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

  // TODO: implement cache in a general way
 // free(cached_rel_param);
 // free(cached_xill_param);

  free_specCache(global_spec_cache);

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

void free_specCache(specCache* spec_cache) {

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

    if (spec_cache->conversion_factor_energyflux != NULL){
      free(spec_cache->conversion_factor_energyflux);
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
