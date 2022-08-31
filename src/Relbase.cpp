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
#include "Relbase.h"
#include "Xillspec.h"
#include "Relphysics.h"

extern "C" {
#include "fftw/fftw3.h"   // assumes installation in heasoft
#include "writeOutfiles.h"
}

// new CACHE routines
cnode *cache_relbase = nullptr;


int save_1eV_pos = 0;


double *global_ener_std = nullptr;

specCache *global_spec_cache = nullptr;


static specCache *new_specCache(int n_cache, int *status) {

  auto *spec = new specCache;

  spec->n_cache = n_cache;
  spec->nzones = 0;
  spec->n_ener = N_ENER_CONV;

  spec->conversion_factor_energyflux = nullptr;

  spec->fft_xill = new double**[n_cache];
  spec->fft_rel = new double**[n_cache];

  spec->fftw_xill = new fftw_complex*[n_cache];
  spec->fftw_rel = new fftw_complex*[n_cache];

  spec->fftw_backwards_input = new fftw_complex[spec->n_ener];
  spec->fftw_output = new double[spec->n_ener];

  spec->plan_c2r = fftw_plan_dft_c2r_1d(spec->n_ener, spec->fftw_backwards_input, spec->fftw_output,FFTW_ESTIMATE);

  spec->xill_spec = new xillSpec*[n_cache];

  int ii;
  int jj;
  int m = 2;
  for (ii = 0; ii < n_cache; ii++) {
    spec->fft_xill[ii] = new double*[m];
    spec->fft_rel[ii] = new double*[m];

    spec->fftw_xill[ii] = new fftw_complex[spec->n_ener];
    spec->fftw_rel[ii] = new fftw_complex[spec->n_ener];

   for (jj = 0; jj < m; jj++) {
      spec->fft_xill[ii][jj] = new double[spec->n_ener];
      spec->fft_rel[ii][jj] = new double[spec->n_ener];
    }
    spec->xill_spec[ii] = nullptr;
  }
  spec->out_spec = nullptr;

  return spec;
}

static void init_specCache(specCache **spec, const int n_zones, int *status) {
  if ((*spec) == nullptr) {
    (*spec) = new_specCache(n_zones, status);
  }
}


specCache *init_global_specCache(int *status) {
  init_specCache(&global_spec_cache, N_ZONES_MAX, status);
  CHECK_RELXILL_ERROR("failed initializing Relconv Spec Cache", status);
  return global_spec_cache;
}

static double* calculate_energyflux_conversion(const double* ener, int n_ener, int* status){

  auto* factor = new double[n_ener];
  CHECK_MALLOC_RET_STATUS(factor, status, nullptr)

  for(int ii=0; ii<n_ener; ii++){
    factor[ii] = 0.5*(ener[ii]+ener[ii+1]) / (ener[ii+1] - ener[ii]);
  }

  return factor;
}








/** @brief FFTW VERSION: convolve the (bin-integrated) spectra f1 and f2 (which need to have a certain binning)
 *  @details fout: gives the output
 *  f1 input (reflection) specrum
 *  f2 filter
 *  ener has length n+1 and is the energy array
 *  requirements: needs "specCache" to be set up
 * **/
void fftw_conv_spectrum(double *ener, const double *fxill, const double *frel, double *fout, int n,
                       int re_rel, int re_xill, int izone, specCache *cache, int *status) {

  CHECK_STATUS_VOID(*status);

  // needs spec cache to be set up
  assert(cache != nullptr);

  if (cache->conversion_factor_energyflux == nullptr){
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

  /**********************************************************************/
  /** cache either the relat. or the xillver part, as only one of the
   * two changes most of the time (reduce time by 1/3 for convolution) **/
  /**********************************************************************/

  /** #1: for the xillver part **/
  if (re_xill) {
    for (ii = 0; ii < n; ii++) {
      cache->fft_xill[izone][0][ii] = fxill[ii] * cache->conversion_factor_energyflux[ii] ;
    }

    fftw_plan plan_xill = fftw_plan_dft_r2c_1d(n, cache->fft_xill[izone][0], cache->fftw_xill[izone],FFTW_ESTIMATE);
    fftw_execute(plan_xill);
    fftw_destroy_plan(plan_xill);
  }

  /** #2: for the relat. part **/
  if (re_rel){
    for (ii = 0; ii < n; ii++) {
      irot = (ii - save_1eV_pos + n) % n;
      cache->fft_rel[izone][0][irot] = frel[ii] * cache->conversion_factor_energyflux[ii];
    }

    fftw_plan plan_rel = fftw_plan_dft_r2c_1d(n, cache->fft_rel[izone][0], cache->fftw_rel[izone],FFTW_ESTIMATE);
    fftw_execute(plan_rel);
    fftw_destroy_plan(plan_rel);
  }

  // complex multiplication (TODO: fix that complex multiplication is not by hand)
  for (ii = 0; ii < n; ii++) {
//    cache->fftw_backwards_input[ii] = cache->fftw_xill[izone][ii] * cache->fftw_rel[izone][ii];
    cache->fftw_backwards_input[ii][0] =
        cache->fftw_xill[izone][ii][0] * cache->fftw_rel[izone][ii][0] -
        cache->fftw_xill[izone][ii][1] * cache->fftw_rel[izone][ii][1];

    cache->fftw_backwards_input[ii][1] =
        cache->fftw_xill[izone][ii][0] * cache->fftw_rel[izone][ii][1] +
            cache->fftw_xill[izone][ii][1] * cache->fftw_rel[izone][ii][0];

  }

  fftw_execute(cache->plan_c2r);

  for (ii = 0; ii < n; ii++) {
    fout[ii] = cache->fftw_output[ii] /  cache->conversion_factor_energyflux[ii]; 
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

  fftw_conv_spectrum(ener, fxill, frel, fout, n, re_rel, re_xill, izone, spec_cache_ptr, status);

  normalizeFFTOutput(ener, fxill, frel, fout, n);

}

void get_relxill_conv_energy_grid(int *n_ener, double **ener, int *status) {
  if (global_ener_std == nullptr) {
    global_ener_std = (double *) malloc((N_ENER_CONV + 1) * sizeof(double));
    CHECK_MALLOC_VOID_STATUS(global_ener_std, status)
    get_log_grid(global_ener_std, (N_ENER_CONV + 1), EMIN_RELXILL_CONV, EMAX_RELXILL_CONV);
  }
  (*n_ener) = N_ENER_CONV;
  (*ener) = global_ener_std;

}


void set_flux_outside_defined_range_to_zero(const double* ener, double* spec, int n_ener, double emin, double emax){
  int warned = 0;
  for (int ii=0; ii<n_ener; ii++){
    if (ener[ii+1]<emin || ener[ii]>emax){
      if (is_debug_run() && warned==0){
        printf(" *** warning: relconv applied outside the allowed energy range %.2f-%.0f\n",
               RELCONV_EMIN, RELCONV_EMAX);
        printf("     -> values outside are set to zero\n\n");
        warned=1;
      }
      spec[ii] = 0;
    }
  }
}


/**
 * @brief basic relconv function: convolve any input spectrum with the relbase kernel
 * @description
 *   it is only defined the in energy range of 0.01-1000 keV (see RELCONV_EMIN, RELCONV_EMAX variables)
 *   and zero outside this range
 * @param double[n_ener_inp+1] ener_inp
 * @param double[n_ener_inp] spec_ener_inp
 *  **/
void relconv_kernel(double *ener_inp, double *spec_inp, int n_ener_inp, relParam *rel_param, int *status) {

  // get the (fixed!) energy grid for a RELLINE for a convolution
  // -> as we do a simple FFT, we can now take into account that we
  // need it to be number = 2^N */
  int n_ener; double *ener;
  get_relxill_conv_energy_grid(&n_ener, &ener, status);

  relline_spec_multizone *rel_profile = relbase(ener, n_ener, rel_param, status);

  // simple convolution only makes sense for 1 zone !
  assert(rel_profile->n_zones == 1);

  auto rebin_flux =  new double[n_ener];
  rebin_spectrum(ener, rebin_flux, n_ener, ener_inp, spec_inp, n_ener_inp);

  specCache* spec_cache = init_global_specCache(status);
  CHECK_STATUS_VOID(*status);
  auto conv_out = new double[n_ener];
  convolveSpectrumFFTNormalized(ener, rebin_flux, rel_profile->flux[0], conv_out, n_ener,
                    1, 1, 0, spec_cache, status);
  CHECK_STATUS_VOID(*status);

  // rebin to the output grid
  rebin_spectrum(ener_inp, spec_inp, n_ener_inp, ener, conv_out, n_ener);

  set_flux_outside_defined_range_to_zero(ener_inp, spec_inp, n_ener_inp, RELCONV_EMIN, RELCONV_EMAX);

  delete[] rebin_flux;
  delete[] conv_out;

}


static void print_reflection_strength(double *ener,
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
  if (is_iongrad_model(rel_param->ion_grad_type) || rel_param->beta > 1e-6) {
    printf(" and beta=%.3f v/c", rel_param->beta);
  }
  printf(" (using boost=1): \n - reflection fraction  %.3f \n - reflection strength is: %.3f \n",
         struct_refl_frac->refl_frac,
         sum / sum_pl);
  printf(" - photons falling into the black hole or plunging region: %.2f%%\n", struct_refl_frac->f_bh * 100);
  printf(" - energy shift from the primary source to the observer is %.3f\n", energy_shift_source_obs(rel_param));
}


void add_primary_component(double *ener, int n_ener, double *flu, relParam *rel_param, xillParam *xill_input_param,
                           RelSysPar *sys_par, int *status) {

  xillTableParam *xill_table_param = get_xilltab_param(xill_input_param, status);
  double *pl_flux = calc_normalized_primary_spectrum(ener, n_ener, rel_param, xill_table_param, status);
  free(xill_table_param);
  CHECK_STATUS_VOID(*status);

  // For the non-relativistic model and if not the LP geometry, we simply multiply by the reflection fraction
  if (is_xill_model(xill_input_param->model_type) || rel_param->emis_type != EMIS_TYPE_LP) {
    for (int ii = 0; ii < n_ener; ii++) {
      flu[ii] *= fabs(xill_input_param->refl_frac);
    }

  } else { // we are in the LP geometry

    assert(rel_param != nullptr);

    lpReflFrac *struct_refl_frac = sys_par->emis->photon_fate_fractions;

    if (xill_input_param->interpret_reflfrac_as_boost) {
      // if set, it is given as boost, wrt predicted refl_frac
      xill_input_param->refl_frac *= struct_refl_frac->refl_frac;
    }

    double g_inf = energy_shift_source_obs(rel_param);
    double prim_fac = struct_refl_frac->f_inf_rest / 0.5 * pow(g_inf, xill_input_param->gam);
    if (rel_param->beta
        > 1e-4) { // flux boost of primary radiation taking into account here (therfore we need f_inf_rest above)
      prim_fac *= pow(doppler_factor_source_obs(rel_param), 2);
    }

    // if the user sets the refl_frac parameter manually, we need to calculate the ratio
    // to end up with the correct normalization
    double norm_fac_refl = (fabs(xill_input_param->refl_frac)) / struct_refl_frac->refl_frac;

    for (int ii = 0; ii < n_ener; ii++) {
      pl_flux[ii] *= prim_fac;
      flu[ii] *= norm_fac_refl;
    }

    /** 5 ** if desired, we ouput the reflection fraction and strength (as defined in Dauser+2016) **/
    if (shouldAuxInfoGetPrinted()) {
      print_reflection_strength(ener, n_ener, flu, rel_param, xill_input_param, pl_flux, struct_refl_frac);
    }

  }

  // Finally, add the power law component if refl_frac >= 0
  if (xill_input_param->refl_frac >= 0) {
    for (int ii = 0; ii < n_ener; ii++) {
      flu[ii] += pl_flux[ii];
    }
  }

  delete[] pl_flux;

}

int did_xill_param_change(const xillParam *cpar, const xillParam *par) {
  if (comp_single_param_val(par->afe, cpar->afe)) {
    return 1;
  }
  if (comp_single_param_val(par->dens, cpar->dens)) {
    return 1;
  }
  if (comp_single_param_val(par->ect, cpar->ect)) {
    return 1;
  }
  if (comp_single_param_val(par->gam, cpar->gam)) {
    return 1;
  }
  if (comp_single_param_val(par->lxi, cpar->lxi)) {
    return 1;
  }
  if (comp_single_param_val(par->kTbb, cpar->kTbb)) {
    return 1;
  }
  if (comp_single_param_val(par->frac_pl_bb, cpar->frac_pl_bb)) {
    return 1;
  }
  if (comp_single_param_val(par->z, cpar->z)) {
    return 1;
  }

  if (comp_single_param_val((double) par->prim_type, (double) cpar->prim_type)) {
    return 1;
  }
  if (comp_single_param_val((double) par->model_type, (double) cpar->model_type)) return 1;

  if (comp_single_param_val(par->iongrad_index, cpar->iongrad_index)) return 1;

  return 0;
}

/* check if values, which need a re-computation of the relline profile, have changed */
int redo_xillver_calc(const relParam *rel_param, const xillParam *xill_param,
                      const relParam *ca_rel_param, const xillParam *ca_xill_param) {

  int redo = 1;

  if ((ca_rel_param != nullptr) && (ca_xill_param != nullptr)) {

    redo = did_xill_param_change(ca_xill_param, xill_param);

    // xillver needs to be re-computed, Ecut changes for the following parameters **/
    if (comp_single_param_val(rel_param->a, ca_rel_param->a) ||
        comp_single_param_val(rel_param->height, ca_rel_param->height) ||
        comp_single_param_val(rel_param->beta, ca_rel_param->beta)) {
      redo = 1;
    }

  }

  return redo;
}

int redo_relbase_calc(const relParam *rel_param, const relParam *ca_rel_param) {

  if (did_rel_param_change(ca_rel_param, rel_param)) {
    return 1;
  } else {
    return 0;
  }

}

/** @brief relbase function calculating the basic relativistic line shape for a given parameter setup
 *  @details
 *    - assuming a 1keV line, by a grid given in keV!
 *    - it is cached
 * input: ener(n_ener), param
 * input: RelSysPar
 * optional input: xillver grid
 * output: photar(n_ener)  [photons/bin]
**/
relline_spec_multizone* relbase_profile(double *ener, int n_ener, relParam *param,
                                       RelSysPar *sysPar,
                                       xillTable *xill_tab,
                                       const double *radialZones,
                                       int nzones,
                                       int *status) {


  inpar* inp = get_inputvals_struct(ener, n_ener, param, status);
  cache_info *ca_info = cli_check_cache(cache_relbase, inp, check_cache_relpar, status);
  relline_spec_multizone *spec = nullptr;

  // set a pointer to the spectrum
  if (is_relbase_cached((ca_info)) == 0) {

    // init the spectra where we store the flux
    param->num_zones = nzones;
    init_relline_spec_multizone(&spec, param, xill_tab, radialZones, &ener, n_ener, status);

    calc_relline_profile(spec, sysPar, status); // returned units are 'photons/bin'

    // normalize it and calculate the angular distribution (if necessary)
    renorm_relline_profile(spec, param, status);

    // last step: store parameters and cached relline_spec_multizone (this prepends a new node to the cache)
    add_relspec_to_cache(&cache_relbase, param, spec, status);
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

  free(inp);
  free(ca_info);

  return spec;
}



/** @brief relbase wrapper function, calculating the relat system params plus the relbase profile
 *  @details
 *    - for more details see relbase_profile function
 *    - uses only a single zone on the disk
 *    - not used for any relxill-case (xill_table=nullptr, as no angular dependency is taken into account)
 * input: ener(n_ener), param
 * optional input: xillver grid
 * output: photar(n_ener)  [photons/bin]
**/
relline_spec_multizone *relbase(double *ener, const int n_ener, relParam *param, int *status) {

  // initialize parameter values (has an internal cache)
  RelSysPar *sysPar = get_system_parameters(param, status);
  CHECK_STATUS_RET(*status, nullptr);
  assert(sysPar != nullptr);

  double *rgrid = get_rzone_grid(param->rin, param->rout, param->num_zones, param->height, status);

  relline_spec_multizone* rel_spec = relbase_profile(ener, n_ener, param, sysPar, nullptr, rgrid, param->num_zones, status);

  delete[] rgrid;

  return rel_spec;
}



void free_rel_cosne(RelCosne *spec) {
  if (spec != nullptr) {
    //	free(spec->ener);  we do not need this, as only a pointer for ener is assigned
    free(spec->cosne);
    if (spec->dist != nullptr) {
      int ii;
      for (ii = 0; ii < spec->n_zones; ii++) {
        free(spec->dist[ii]);
      }
    }
    free(spec->dist);
    free(spec);
  }
}

void free_rel_spec(relline_spec_multizone *spec) {
  if (spec != nullptr) {
    free(spec->ener);
    free(spec->rgrid);
    if (spec->flux != nullptr) {
      int ii;
      for (ii = 0; ii < spec->n_zones; ii++) {
        if (spec->flux[ii] != nullptr) {
          free(spec->flux[ii]);
        }
      }
    }
    free(spec->flux);
    if (spec->rel_cosne != nullptr) {
      free_rel_cosne(spec->rel_cosne);
    }
    free(spec);
  }
}

void free_cached_tables() {
  free_relprofile_cache();

  free_cached_relTable();
  free_cached_lpTable();
  free_cached_xillTable();

  // TODO: implement cache in a general way
 // free(cached_rel_param);
 // free(cached_xill_param);

  free_specCache(global_spec_cache);

  free(global_ener_std);
  // free(global_ener_xill); // TODO, implement free of this global energy grid

}

void free_fft_cache(double ***sp, int n1, int n2) {

  int ii;
  int jj;
  if (sp != nullptr) {
    for (ii = 0; ii < n1; ii++) {
      if (sp[ii] != nullptr) {
        for (jj = 0; jj < n2; jj++) {
          free(sp[ii][jj]);
        }
      }
      free(sp[ii]);
    }
    free(sp);
  }

}

spectrum *new_spectrum(int n_ener, const double *ener, int *status) {

  auto *spec = new spectrum;
  spec->n_ener = n_ener;
  spec->ener = new double[n_ener];
  spec->flux = new double[n_ener];

  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    spec->ener[ii] = ener[ii];
    spec->flux[ii] = 0.0;
  }

  return spec;
}

void free_spectrum(spectrum *spec) {
  if (spec != nullptr) {
    delete[]spec->ener;
    delete[] spec->flux;
    delete spec;
  }
}

void free_fftw_complex_cache(fftw_complex** val, int n){
  for(int ii=0; ii<n; ii++){
    fftw_free(val[ii]);
  }
}

void free_specCache(specCache* spec_cache) {

  int ii;
  int m = 2;
  if (spec_cache != nullptr) {
    if (spec_cache->xill_spec != nullptr) {
      for (ii = 0; ii < spec_cache->n_cache; ii++) {
        if (spec_cache->xill_spec[ii] != nullptr) {
          free_xill_spec(spec_cache->xill_spec[ii]);
        }
      }
      free(spec_cache->xill_spec);
    }

    if (spec_cache->fft_xill != nullptr) {
      free_fft_cache(spec_cache->fft_xill, spec_cache->n_cache, m);
    }

    if (spec_cache->fftw_rel != nullptr) {
      free_fft_cache(spec_cache->fft_rel, spec_cache->n_cache, m);
    }

    free_fftw_complex_cache(spec_cache->fftw_rel, spec_cache->n_cache);
    free_fftw_complex_cache(spec_cache->fftw_xill, spec_cache->n_cache);
    fftw_free(spec_cache->fftw_backwards_input);
    fftw_destroy_plan(spec_cache->plan_c2r);
    free(spec_cache->fftw_output);

    if (spec_cache->conversion_factor_energyflux != nullptr){
      free(spec_cache->conversion_factor_energyflux);
    }

    free_spectrum(spec_cache->out_spec);

  }

  free(spec_cache);

}

/** free the CLI cache **/

void free_cache() {
  free_cache_syspar();
  cli_delete_list(&cache_relbase);
}


