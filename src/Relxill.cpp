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
#include "Relxill.h"

#include "IonGradient.h"
#include "XspecSpectrum.h"
#include "Relreturn_Corona.h"

extern "C" {
#include "xilltable.h"
}

/** caching parameters **/
relParam *cached_rel_param = nullptr;
xillParam *cached_xill_param = nullptr;

///////////////////////////////////////
// Forward Definitions of Functions  //
///////////////////////////////////////

double calc_ecut_at_primary_source(const xillParam *xill_param,
                                   const relParam *rel_param,
                                   double ecut_input,
                                   int *status);
void free_arrays_relxill_kernel(int n_zones,
                                const double *conv_out,
                                const double *single_spec_inp,
                                double *const *xill_angledep_spec);



///////////////////////////////////////
// Function Definitions              //
///////////////////////////////////////

void relxill_convolution_multizone(const XspecSpectrum &spectrum,
                                   const relline_spec_multizone *rel_profile,
                                   specCache *spec_cache,
                                   const relParam *rel_param,
                                   const CachingStatus &caching_status,
                                   int *status);
/**
 * check if the complete spectrum is cached and if the energy grid did not change
 *  -> additionally we adapt the energy grid to the new dimensions and energy values
 *     if it's not the case yet
 */
static void check_caching_energy_grid(CachingStatus &caching_status, specCache *cache, const XspecSpectrum &spectrum) {

  caching_status.energy_grid = cached::no;

  if (cache->out_spec == nullptr) {
    return;
  } else {
    // first need to check if the energy grid has a different number of bins
    if (cache->out_spec->n_ener != spectrum.num_flux_bins()) {
      return;
    }

    for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) { // now check if the energies themselves changes
      if (fabs(cache->out_spec->ener[ii] - spectrum.energy[ii]) > 1e-4) {
        //        for (int jj = 0; jj < spectrum.num_flux_bins(); jj++) {
        //          cache->out_spec->ener[jj] = spectrum.energy[jj];
        //        }
        return;
      }
    }
    // no nullptr, same number of bins and energies are the same
    caching_status.energy_grid = cached::yes;
  }


}

/**
 * @brief test if the rel_param number of zones is different from the global cache
 * @param rel_param
 * @return bool
 */
static int did_number_of_zones_change(const relParam *rel_param) {

  if (cached_rel_param == nullptr) {
    return 1;
  }

  if (rel_param->num_zones != cached_rel_param->num_zones) {
    if (is_debug_run()) {
      printf("  *** warning :  the number of radial zones was changed from %i to %i \n",
             rel_param->num_zones, cached_rel_param->num_zones);
    }
    return 1;
  } else {
    return 0;
  }
}

static void check_caching_parameters(CachingStatus &caching_status,
                                     const relParam *rel_param,
                                     const xillParam *xill_param) {

  // special case: no caching if output files are to be written
  if (shouldOutfilesBeWritten() || did_number_of_zones_change(rel_param)) {
    caching_status.relat = cached::no;
    caching_status.xill = cached::no;
  } else {

    /** did any of the parameters change?  **/
    if (redo_relbase_calc(rel_param, cached_rel_param) == 0) {
      caching_status.relat = cached::yes;
    }
    if (redo_xillver_calc(rel_param, xill_param, cached_rel_param, cached_xill_param) == 0) {
      caching_status.xill = cached::yes;
    }

    // special case: as return_rad emissivity depends on xillver, if xillver not cached, relat also no cached
    if (rel_param->return_rad && caching_status.xill == cached::no){
      caching_status.relat = cached::no;
    }

  }

}

static void write_output_spec_zones(const XspecSpectrum &spectrum, double *spec_inp_single, int ii, int *status) {
  char vstr[200];
  if (sprintf(vstr, "test_relxill_spec_zones_%03i.dat", ii + 1) == -1) {
    RELXILL_ERROR("failed to get filename", status);
  }
  save_xillver_spectrum(spectrum.energy, spec_inp_single, spectrum.num_flux_bins(), vstr);
}

/** initialize the cached output spec array **/
void copy_spectrum_to_cache(const XspecSpectrum &spectrum,
                            specCache *spec_cache,
                            int *status) {

  CHECK_STATUS_VOID(*status);

  if ((spec_cache->out_spec != nullptr)) {
    if (spec_cache->out_spec->n_ener != spectrum.num_flux_bins()) {
      free_out_spec(spec_cache->out_spec);
      spec_cache->out_spec = init_out_spec(spectrum.num_flux_bins(), spectrum.energy, status);
    }
  } else {
    spec_cache->out_spec = init_out_spec(spectrum.num_flux_bins(), spectrum.energy, status);
  }

  for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) {
    spec_cache->out_spec->flux[ii] = spectrum.flux[ii];
  }
}


/**
 * @brief calculate the correction factors for the given radial grid
 * @details only calculated if the return_rad!=0 is set, otherwise a nullptr
 *  is returned, which means no correction is calculated
 * @param rel_param
 * @param xill_spec
 * @param rgrid
 * @param xill_table_param
 * @param status
 * @return
 */
rradCorrFactors* calc_rrad_corr_factors(xillSpec **xill_spec, const RadialGrid &rgrid,
                                       xillTableParam *const *xill_table_param, int *status) {

  rradCorrFactors* rrad_corr_factors = init_rrad_corr_factors(rgrid.radius, rgrid.num_zones);

  if (rrad_corr_factors == nullptr || *status != EXIT_SUCCESS) {
    assert(rrad_corr_factors != nullptr);
    return rrad_corr_factors;
  }

  for (int ii = 0; ii < rgrid.num_zones; ii++) {
    get_xillver_fluxcorrection_factors(xill_spec[ii],
                                       &(rrad_corr_factors->corrfac_flux[ii]),
                                       &(rrad_corr_factors->corrfac_gshift[ii]),
                                       xill_table_param[ii], status);
  }

  return rrad_corr_factors;
}

static void free_xill_table_param_array(const relParam *rel_param, xillTableParam *const *xill_table_param) {
  for (int ii = 0; ii < rel_param->num_zones; ii++) {
    free(xill_table_param[ii]);
  }
  delete[] xill_table_param;
}


///////////////////////////////////////
// MAIN: Relxill Kernel Function     //
///////////////////////////////////////


/** @brief convolve a xillver spectrum with the relbase kernel
 * @param spectrum
 * @param xill_param
 * @param rel_param
 * @param status
 */
void relxill_kernel(const XspecSpectrum &spectrum,
                    const ModelParams &params,
                    int *status) {

  relParam *rel_param = get_rel_params(params);
  xillParam *xill_param = get_xill_params(params);


  // in case of an ionization gradient, we need to update the number of zones
  assert(rel_param->num_zones == get_num_zones(rel_param->model_type, rel_param->emis_type, rel_param->ion_grad_type));


  specCache *spec_cache = init_global_specCache(status);
  assert(spec_cache != nullptr);
  CHECK_STATUS_VOID(*status);

  auto caching_status = CachingStatus();
  check_caching_parameters(caching_status, rel_param, xill_param);
  check_caching_energy_grid(caching_status, spec_cache, spectrum);

  RelSysPar *sys_par = get_system_parameters(rel_param, status);
  CHECK_STATUS_VOID(*status);

  if (caching_status.is_all_cached()) { // already cached, simply use the cached output flux value
    for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) {
      spectrum.flux[ii] = spec_cache->out_spec->flux[ii];
    }

  } else { // if NOT, we need to do a whole lot of COMPUTATIONS

    // stored the parameters for which we are calculating
    set_cached_xill_param(xill_param, &cached_xill_param, status);
    set_cached_rel_param(rel_param, &cached_rel_param, status);
    CHECK_STATUS_VOID(*status);

    // --- 0 ---
    auto radial_grid = RadialGrid(rel_param->rin, rel_param->rout, rel_param->num_zones, rel_param->height);

    // --- 1 --- calculate disk irradiation


    // --- 2 --- calculate ionization gradient
    IonGradient ion_gradient{radial_grid, rel_param->ion_grad_type};
    ion_gradient.calculate(*(sys_par->emis), xill_param);

    auto xill_param_zone = new xillTableParam *[rel_param->num_zones];

    // --- 3 --- calculate xillver reflection spectra  (for every zone)

    double ecut_primary = calc_ecut_at_primary_source(xill_param, rel_param, xill_param->ect, status);

    for (int ii = 0; ii < rel_param->num_zones; ii++) {
      xill_param_zone[ii] = get_xilltab_param(xill_param, status);

      // set xillver parameters for the given zone
      xill_param_zone[ii]->ect = ion_gradient.get_ecut_disk_zone(rel_param, ecut_primary, ii);
      xill_param_zone[ii]->lxi = ion_gradient.lxi[ii];
      xill_param_zone[ii]->dens = ion_gradient.dens[ii];

      // --- 3a: load xillver spectra
      // (always need to re-compute for an ionization gradient, TODO: can we do better caching?)
      if (caching_status.xill == cached::no) {
        if (spec_cache->xill_spec[ii] != nullptr) {
          free_xill_spec(spec_cache->xill_spec[ii]);
        }
        spec_cache->xill_spec[ii] = get_xillver_spectra_table(xill_param_zone[ii], status);
      }
    }

    // -- 4 -- returning radiation correction factors (only calculated if above a given threshold)
    rel_param->rrad_corr_factors =
        (rel_param->return_rad != 0 && rel_param->a > SPIN_MIN_RRAD_CALC_CORRFAC) ?
        calc_rrad_corr_factors(spec_cache->xill_spec, radial_grid, xill_param_zone, status) :
        nullptr;

    free_xill_table_param_array(rel_param, xill_param_zone);

    // --- 5 --- calculate multi-zone relline profile
    xillTable *xill_tab = nullptr; // needed for the calc_relline_profile call
    get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, status);

    // --- 5a --- calculate the emissivity including the rrad correction factors
    sys_par = get_system_parameters(rel_param, status); // no need to free this, is automatically done by the cache

    // --- 5b --- calculate the relbase profile
    int n_ener_conv; // energy grid for the convolution, only created
    double *ener_conv;
    get_relxill_conv_energy_grid(&n_ener_conv, &ener_conv, status);
    relline_spec_multizone *rel_profile =
        relbase_profile(ener_conv, n_ener_conv, rel_param, sys_par, xill_tab,
                        ion_gradient.radial_grid.radius, ion_gradient.nzones(), status);

    // --- 6 --- convolve the reflection with the relativistic kernel
    relxill_convolution_multizone(spectrum, rel_profile, spec_cache, rel_param, caching_status, status);

    copy_spectrum_to_cache(spectrum, spec_cache, status);
    free_rrad_corr_factors(&(rel_param->rrad_corr_factors));
  }

  add_primary_component(spectrum.energy, spectrum.num_flux_bins(), spectrum.flux, rel_param, xill_param,
                        sys_par, status);

  delete rel_param;
  delete xill_param;

}



void relxill_convolution_multizone(const XspecSpectrum &spectrum,
                                   const relline_spec_multizone *rel_profile,
                                   specCache *spec_cache,
                                   const relParam *rel_param,
                                   const CachingStatus &caching_status,
                                   int *status) {

  CHECK_STATUS_VOID(*status);

  const int n_ener_conv = rel_profile->n_ener;
  double* ener_conv = rel_profile->ener;

  auto conv_out = new double[n_ener_conv];
  auto single_spec_inp = new double[spectrum.num_flux_bins()];
  auto xill_angledep_spec = new double *[rel_param->num_zones];

  // make sure the output array is set to 0
  for (int ie = 0; ie < spectrum.num_flux_bins(); ie++) {
    spectrum.flux[ie] = 0.0;
  }

  for (int ii = 0; ii < rel_profile->n_zones; ii++) { /***** loop over ionization zones   ******/

    /** avoid problems where no relxill bin falls into an ionization bin **/
    if (calcSum(rel_profile->flux[ii], rel_profile->n_ener) < 1e-12) {
      continue;
    }

    // -- 6 -- get angle-dependent spectrum
    xill_angledep_spec[ii] = new double[n_ener_conv];
    get_xillver_angdep_spec(xill_angledep_spec[ii], ener_conv, n_ener_conv, rel_profile->rel_cosne->dist[ii],
                            spec_cache->xill_spec[ii], status);
     CHECK_STATUS_VOID(*status);

    // -- 7 -- convolve the spectrum **
    int recompute_xill = 1; // always recompute fft for xillver, as relat changes the angular distribution
    convolveSpectrumFFTNormalized(ener_conv, xill_angledep_spec[ii], rel_profile->flux[ii], conv_out, n_ener_conv,
                                  caching_status.recomput_relat(), recompute_xill, ii, spec_cache, status);
    CHECK_STATUS_VOID(*status);
    rebin_spectrum(spectrum.energy, single_spec_inp, spectrum.num_flux_bins(), ener_conv, conv_out, n_ener_conv);

    // add it to the final output spectrum
    for (int jj = 0; jj < spectrum.num_flux_bins(); jj++) {
      spectrum.flux[jj] += single_spec_inp[jj];
    }

    if (is_debug_run() && rel_profile->n_zones <= 10) {
      write_output_spec_zones(spectrum, single_spec_inp, ii, status);
    }

  } /**** END OF LOOP OVER RADIAL ZONES *****/

  free_arrays_relxill_kernel(rel_profile->n_zones, conv_out, single_spec_inp, xill_angledep_spec);
}

/**
 * @brief free arrays allocated and needed by relxill_kernel
 * @param n_zones
 * @param conv_out
 * @param single_spec_inp
 * @param xill_angledep_spec
 */
void free_arrays_relxill_kernel(int n_zones,
                                const double *conv_out,
                                const double *single_spec_inp,
                                double *const *xill_angledep_spec) {
  delete[] single_spec_inp;
  delete[] conv_out;
  for (int ii = 0; ii < n_zones; ii++) {
    delete[] xill_angledep_spec[ii];
  }
  delete[] xill_angledep_spec;
}

/** @brief calculate the ecut/kTe value at the primary source (from the xspec parameter input value)
 *  @details in case of the nthcomp model Ecut is already input in the frame of the primary source
 *   but for the bkn_powerlaw it is given in the observer frame. In case of a BB spectrum, it is also
 *   given in the source frame.
 * @param xill_param
 * @param rel_param
 * @param ecut_input
 * @param status
 * @return
 */
double calc_ecut_at_primary_source(const xillParam *xill_param,
                                   const relParam *rel_param,
                                   double ecut_input,
                                   int *status) {
  if (xill_param->prim_type == PRIM_SPEC_ECUT) {
    return ecut_input / energy_shift_source_obs(rel_param);
  } else if (xill_param->prim_type == PRIM_SPEC_NTHCOMP || xill_param->prim_type == PRIM_SPEC_BB) {
    return ecut_input;
  } else {
    RELXILL_ERROR("unknown primary spectrum type, can not calculate Ecut at the primary source", status);
    return -1;
  }
}


