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
#include "PrimarySource.h"

extern "C" {
#include "xilltable.h"
}

/** caching parameters **/
relParam *cached_rel_param = nullptr;
xillParam *cached_xill_param = nullptr;


///////////////////////////////////////
// Function Definitions              //
///////////////////////////////////////

void relxill_convolution_multizone(const XspecSpectrum &spectrum,
                                   const relline_spec_multizone *rel_profile,
                                   const SpectrumZones &xill_spec_zones,
                                   specCache *spec_cache,
                                   const relParam *rel_param,
                                   const CachingStatus &caching_status,
                                   int *status);


/**
 * check if the complete spectrum is cached and if the energy grid did not change
 *  -> additionally we adapt the energy grid to the new dimensions and energy values
 *     if it's not the case yet
 */
void CachingStatus::check_caching_energy_grid(specCache *cache, const XspecSpectrum &spectrum) {

  m_cache_energy_grid = cached::no;

  if (cache->out_spec == nullptr) {
    return;
  }

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
  m_cache_energy_grid = cached::yes;

}


/**
 * @brief test if the m_rel_param number of zones is different from the global cache
 * @param rel_param
 * @return bool
 */
auto CachingStatus::did_number_of_zones_change(const relParam *rel_param) -> int {

  if (cached_rel_param == nullptr) {
    return 1;
  }

  if (rel_param->num_zones != cached_rel_param->num_zones) {
    if (is_debug_run() != 0) {
      printf("  *** warning :  the number of radial zones was changed from %i to %i \n",
             rel_param->num_zones, cached_rel_param->num_zones);
    }
    return 1;
  }

  return 0;
}


void CachingStatus::check_caching_parameters(const relParam *rel_param, const xillParam *xill_param) {
  // special case: no caching if output files are to be written
  if ((shouldOutfilesBeWritten() != 0) || (did_number_of_zones_change(rel_param) != 0)) {
    m_cache_relat = cached::no;
    m_cache_xill = cached::no;
  } else {

    /** did any of the parameters change?  **/
    if (redo_relbase_calc(rel_param, cached_rel_param) == 0) {
      m_cache_relat = cached::yes;
    }
    if (redo_xillver_calc(rel_param, xill_param, cached_rel_param, cached_xill_param) == 0) {
      m_cache_xill = cached::yes;
    }

    // special case: as return_rad emissivity depends on xillver, if xill_param as well as rel_paramare not cached
    if ((rel_param->return_rad != 0) && m_cache_xill == cached::no) {
      m_cache_relat = cached::no;
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
      free_spectrum(spec_cache->out_spec);
      spec_cache->out_spec = new_spectrum(spectrum.num_flux_bins(), spectrum.energy, status);
    }
  } else {
    spec_cache->out_spec = new_spectrum(spectrum.num_flux_bins(), spectrum.energy, status);
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

  rradCorrFactors *rrad_corr_factors = init_rrad_corr_factors(rgrid);

  if (rrad_corr_factors == nullptr || *status != EXIT_SUCCESS) {
    assert(rrad_corr_factors != nullptr);
    return rrad_corr_factors;
  }

  for (int ii = 0; ii < rgrid.num_zones(); ii++) {
    get_xillver_fluxcorrection_factors(xill_spec[ii],
                                       &(rrad_corr_factors->corrfac_flux[ii]),
                                       &(rrad_corr_factors->corrfac_gshift[ii]),
                                       xill_table_param[ii], status);
  }

  return rrad_corr_factors;
}

static void free_xill_table_param_array(const relParam *rel_param, xillTableParam *const *xill_table_param) {
  for (int ii = 0; ii < rel_param->num_zones; ii++) {
    delete xill_table_param[ii];
  }
  delete[] xill_table_param;
}


void get_relxill_params(const ModelParams &params, relParam *&rel_param, xillParam *&xill_param) {
  rel_param = get_rel_params(params);
  xill_param = get_xill_params(params);

  // special case, for the LP model and the Ecut model, the cutoff energy is given in the observer frame
  // -> change it such that the parameter "ect" is always given in the source frame
  if (rel_param->emis_type == EMIS_TYPE_LP && xill_param->prim_type == PRIM_SPEC_ECUT) {
    xill_param->ect /= energy_shift_source_obs(rel_param);
  }

}

/*
 * @brief: calculate the xillver reflection spectra for the parameter array given as input
 *
 * @details:
 *  - the spec_cache structure has the spec_cache->xill_spec structure allocated with the maximal number of allowed zones
 *  - if the xillver parameters did not change, it will be re-used and not re-calculated
 */
xillSpec **get_xillver_reflection_spectra(specCache *spec_cache,
                                          xillTableParam **xill_param_zone,
                                          int nzones,
                                          cached caching_status_xill) {
  auto xill_spec = spec_cache->xill_spec;

  int status = EXIT_SUCCESS;
  // --- 3 --- calculate xillver reflection spectra  (for every zone)
  for (int ii = 0; ii < nzones; ii++) {
    if (caching_status_xill ==
        cached::no) { // note, also means if alpha_model or height, spin, beta changes it will be re-computed
      if (xill_spec[ii]
          != nullptr) { // as the spectrum is not cached, free the memory to be able to load a new spectrum
        free_xill_spec(xill_spec[ii]);
      }
      xill_spec[ii] = get_xillver_spectra_table(xill_param_zone[ii], &status);
    }
  }

  if (status != EXIT_SUCCESS) {
    throw std::exception();
  }

  return xill_spec;
}


///////////////////////////////////////
// MAIN: Relxill Kernel Function     //
///////////////////////////////////////

/** @brief convolve a xillver spectrum with the relbase kernel
 */
void relxill_kernel(const XspecSpectrum &spectrum,
                    const ModelParams &params,
                    int *status) {

  relParam *rel_param = nullptr;
  xillParam *xill_param = nullptr;
  get_relxill_params(params, rel_param, xill_param);

  // in case of an ionization gradient, we need to update the number of zones, make sure they are set correctly
  assert(rel_param->num_zones == get_num_zones(rel_param->model_type, rel_param->emis_type, rel_param->ion_grad_type));

  specCache *spec_cache = init_global_specCache(status);
  assert(spec_cache != nullptr);

  auto caching_status = CachingStatus(rel_param, xill_param, spec_cache, spectrum);

  RelSysPar *sys_par = get_system_parameters(rel_param, status);
  auto primary_source = PrimarySource(params, sys_par);

  if (caching_status.is_all_cached() == 1) { // if already cached, simply use the cached output flux value
    for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) {
      spectrum.flux[ii] = spec_cache->out_spec->flux[ii];
    }

  } else {
    // store the parameters for which we are calculating
    set_cached_xill_param(xill_param, &cached_xill_param, status);
    set_cached_rel_param(rel_param, &cached_rel_param, status);
    CHECK_STATUS_VOID(*status);

    // --- 1 --- calculate the accretion disk zones and set their parameters
    IonGradient ion_gradient{RadialGrid(rel_param->rin, rel_param->rout, rel_param->num_zones, rel_param->height),
                             rel_param->ion_grad_type, xill_param->iongrad_index};
    ion_gradient.calculate_gradient(*(sys_par->emis), primary_source.source_parameters);

    auto xill_param_zone =
        ion_gradient.get_xill_param_zone(primary_source.source_parameters.xilltab_param());

    // --- 2 --- get xillver reflection spectra (are internally stored in a general, cached structure "SpecCache")
    //           such that they are re-used of the caching_status.m_cache_xill==yes
    auto xill_refl_spectra_zone =
        get_xillver_reflection_spectra(spec_cache, xill_param_zone, ion_gradient.nzones(), caching_status.xill());

    // -- 3 -- returning radiation correction factors (only calculated if above a given threshold)
    rel_param->rrad_corr_factors =
        (rel_param->return_rad != 0 && rel_param->a > SPIN_MIN_RRAD_CALC_CORRFAC) ?
        calc_rrad_corr_factors(xill_refl_spectra_zone, ion_gradient.radial_grid, xill_param_zone, status) :
        nullptr;

    //  calculate the emissivity including the rrad correction factors (for those the disk parameters need to be known)
    sys_par = get_system_parameters(rel_param, status); // no need to free this, is automatically done by the cache

    // --- 4 --- calculate multi-zone relline profile
    xillTable *xill_tab = nullptr; // needed for the relbase_profile call
    get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, status);

    // calculate the relline profile
    int n_ener_conv; // energy grid for the convolution, only created
    double *ener_conv = nullptr;
    get_relxill_conv_energy_grid(&n_ener_conv, &ener_conv);
    relline_spec_multizone *rel_profile =
        relbase_profile(ener_conv, n_ener_conv, rel_param, sys_par, xill_tab,
                        ion_gradient.radial_grid.radius.data(),
                        static_cast<int>(ion_gradient.radial_grid.num_zones()),
                        status);

    // --- 5 --- calculate the xillver spectra depending on the angular distribution (stored in the rel_profile)
    auto xillver_spectra_zones =
        SpectrumZones(xill_refl_spectra_zone[0]->ener, xill_refl_spectra_zone[0]->n_ener, ion_gradient.nzones());
    for (int ii = 0; ii < ion_gradient.nzones(); ii++) {
      calc_xillver_angdep(xillver_spectra_zones.flux[ii],
                          xill_refl_spectra_zone[ii],
                          rel_profile->rel_cosne->dist[ii],
                          status);

    }

    // need to re-normalize the spectra due to the energy shift from the source to the disk
    // reason: xillver is defined on a fixed energy flux integrated from 0.1-1000keV (see Dauser+16, A1), therefore
    // shifting ecut/kTe in energy will change the normalization of the primary spectrum, which was used to calculate
    // the reflected spectrum. As the normalization of reflection is calculated for the normalized incident spectrum
    // on the disk, we need to correct for the change in normalization, in order for the incident spectrum matching
    // the normalization of the primary source spectrum

    // we need to calculate the normalization change from disk to source, therefore calculate from source to disk and take
    // the inverse
    auto norm_change_factors = calc_xillver_normalization_change_source_to_disk(
        ion_gradient.m_energy_shift_source_disk, ion_gradient.nzones(), primary_source.source_parameters.xilltab_param()
    );
    for (int ii = 0; ii < ion_gradient.nzones(); ii++) {
      for (int jj = 0; jj < xillver_spectra_zones.num_flux_bins; jj++) {
        xillver_spectra_zones.flux[ii][jj] /= norm_change_factors[ii];
      }

    }
    delete[] norm_change_factors;

    free_xill_table_param_array(rel_param, xill_param_zone);

    // --- 6 --- convolve the reflection with the relativistic kernel
    relxill_convolution_multizone(spectrum,
                                  rel_profile,
                                  xillver_spectra_zones,
                                  spec_cache,
                                  rel_param,
                                  caching_status,
                                  status);

    copy_spectrum_to_cache(spectrum, spec_cache, status);
    free_rrad_corr_factors(&(rel_param->rrad_corr_factors));
  }

  primary_source.add_primary_spectrum(spectrum);

  delete rel_param;
  delete xill_param;

}


void relxill_convolution_multizone(const XspecSpectrum &spectrum,
                                   const relline_spec_multizone *rel_profile,
                                   const SpectrumZones &xill_spec_zones,
                                   specCache *spec_cache,
                                   const relParam *rel_param,
                                   const CachingStatus &caching_status,
                                   int *status) {

  CHECK_STATUS_VOID(*status);

  const int n_ener_conv = rel_profile->n_ener;
  double* ener_conv = rel_profile->ener;

  auto *conv_out = new double[n_ener_conv];
  auto *xill_rebinned_spec = new double[n_ener_conv];
  auto single_spec_inp = new double[spectrum.num_flux_bins()];

  // make sure the output array is set to 0
  for (int ie = 0; ie < spectrum.num_flux_bins(); ie++) {
    spectrum.flux[ie] = 0.0;
  }

  for (int ii = 0; ii < rel_profile->n_zones; ii++) { /***** loop over ionization zones   ******/

    /** avoid problems where no relxill bin falls into an ionization bin **/
    if (calcSum(rel_profile->flux[ii], rel_profile->n_ener) < 1e-12) {
      continue;
    }

    rebin_spectrum(ener_conv, xill_rebinned_spec, n_ener_conv,
                   xill_spec_zones.energy() , xill_spec_zones.flux[ii], xill_spec_zones.num_flux_bins);

    // --2-- convolve the spectrum on the energy grid "ener_conv" **
    // always recompute fft for xillver if any parameter changedrelat changed, as m_cache_relat changes the angular distribution
    int recompute_xill = caching_status.any_parameter_changed();
    convolveSpectrumFFTNormalized(ener_conv, xill_rebinned_spec, rel_profile->flux[ii], conv_out, n_ener_conv,
                                  caching_status.recomput_relat(), recompute_xill, ii, spec_cache, status);
    CHECK_STATUS_VOID(*status);
    rebin_spectrum(spectrum.energy, single_spec_inp, spectrum.num_flux_bins(), ener_conv, conv_out, n_ener_conv);

    // --3-- add it to the final output spectrum
    for (int jj = 0; jj < spectrum.num_flux_bins(); jj++) {
      spectrum.flux[jj] += single_spec_inp[jj];
    }

    if (is_debug_run() && rel_profile->n_zones <= 10) {
      write_output_spec_zones(spectrum, single_spec_inp, ii, status);
    }

  } /**** END OF LOOP OVER RADIAL ZONES *****/

  // free arrays allocated and needed by relxill_kernel
  delete[] single_spec_inp;
  delete[] conv_out;
  delete[] xill_rebinned_spec;
}
