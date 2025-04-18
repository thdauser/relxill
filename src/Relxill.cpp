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
#include "ModelDefinition.h"

#include <chrono>

extern "C" {
#include "xilltable.h"
}

/** caching parameters **/
relParam *cached_rel_param = nullptr;
xillParam *cached_xill_param = nullptr;


///////////////////////////////////////
// Function Definitions              //
///////////////////////////////////////

void relxill_convolution_multizone(const RelxillSpec &relxill_spec,
                                   const relline_spec_multizone *rel_profile,
                                   const SpectrumZones &xill_spec_zones,
                                   specCache *spec_cache,
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

static void write_output_spec_zones(const RelxillSpec &spectrum, int ii, int *status) {
  char vstr[200];
  if (sprintf(vstr, "test_relxill_spec_zones_%03i.dat", ii + 1) == -1) {
    RELXILL_ERROR("failed to get filename", status);
  }
  save_xillver_spectrum(spectrum.energy(), spectrum.flux, spectrum.num_flux_bins, vstr);
}

/*
/** initialize the cached output spec array
 *  DEPRECATED: need to remove this part of the cache
 **/
void copy_spectrum_to_cache(const Spectrum &spectrum,
                            specCache *spec_cache,
                            int *status) {

  CHECK_STATUS_VOID(*status);

  spec_cache->out_spec = nullptr;


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

  for (size_t ii = 0; ii < rgrid.num_zones(); ii++) {
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


void get_relxill_params(const ModelDefinition &params, relParam *&rel_param, xillParam *&xill_param) {
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





bool do_renorm_relxill() {
  if (is_env_set("RELXILL_RENORMALIZE", 0) == 1){
    return true;
  } else {
    return false;
  }
}

void renorm_relxill_spectrum_1keV(const Spectrum &spectrum) {

  // normalize it to have 1 cts/sec/keV/cm2 at ref_value
  double ref_value = 3.0;
  int ind = binary_search(spectrum.energy(), spectrum.num_flux_bins, ref_value);

  double dE = spectrum.energy()[ind+1] - spectrum.energy()[ind];
  double norm_factor = 1.0 / (spectrum.flux[ind] / dE) ;  // flux is in cts/bin, need to convert it to cts/keV

  spectrum.multiply_flux_by(norm_factor);
}

void rebin_spectrum(const XspecSpectrum &spectrum, const RelxillSpec &relxill_spec) {
  _rebin_spectrum(spectrum.energy, spectrum.flux, spectrum.num_flux_bins(),
                  relxill_spec.energy(), relxill_spec.flux, relxill_spec.num_flux_bins);
}

// rebin to input energy grid (as given by Xspec)
void rebin_and_normalize_relxill_for_xspec(const XspecSpectrum &spectrum, const RelxillSpec &relxill_spec) {

  if (do_renorm_relxill()){
    auto relxill_renorm_spec = RelxillSpec(relxill_spec.flux);
    renorm_relxill_spectrum_1keV(relxill_renorm_spec);

    rebin_spectrum(spectrum, relxill_renorm_spec);

  } else {
    rebin_spectrum(spectrum, relxill_spec);
  }
}



///////////////////////////////////////
// MAIN: Relxill Kernel Function     //
///////////////////////////////////////


/** @brief convolve a xillver xspec_spectrum with the relbase kernel
 */
void relxill_kernel(const XspecSpectrum &xspec_spectrum,
                    const ModelDefinition &params,
                    int *status) {

  auto tstart = std::chrono::steady_clock::now();

  // check if we find the evaluation in the cache
  auto cached_elem = RelxillCache::instance().find_spec_pair(params);
  auto relxill_spec = RelxillCache::get_spec(cached_elem);

  // calculated, if it is not cached
  if (!RelxillCache::is_cached(cached_elem)) {

    relParam *rel_param = nullptr;
    xillParam *xill_param = nullptr;
    get_relxill_params(params, rel_param, xill_param);

    // in case of an ionization gradient, we need to update the number of zones, make sure they are set correctly
    assert(
        rel_param->num_zones == get_num_zones(rel_param->model_type, rel_param->emis_type, rel_param->ion_grad_type));

    specCache *spec_cache = init_global_specCache(status);
    assert(spec_cache != nullptr);

    auto caching_status = CachingStatus(rel_param, xill_param, spec_cache, xspec_spectrum);

    RelSysPar *sys_par = get_system_parameters(rel_param, status);
    auto primary_source = PrimarySource(params, sys_par);


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
    relline_spec_multizone *rel_profile =
        relbase_profile(relxill_spec.energy(), relxill_spec.num_flux_bins, rel_param, sys_par, xill_tab,
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
    // shifting ecut/kTe in energy will change the normalization of the primary xspec_spectrum, which was used to calculate
    // the reflected xspec_spectrum. As the normalization of reflection is calculated for the normalized incident xspec_spectrum
    // on the disk, we need to correct for the change in normalization, in order for the incident xspec_spectrum matching
    // the normalization of the primary source xspec_spectrum

    // we need to calculate the normalization change from disk to source, therefore calculate from source to disk and take
    // the inverse
    auto *norm_change_factors = calc_xillver_normalization_change_source_to_disk(
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
    relxill_convolution_multizone(relxill_spec, rel_profile, xillver_spectra_zones, spec_cache, caching_status, status);

    //   copy_spectrum_to_cache(relxill_spec, spec_cache, status);
    free_rrad_corr_factors(&(rel_param->rrad_corr_factors));


    // add the primary source xspec_spectrum
    primary_source.add_primary_spectrum(relxill_spec);

    RelxillCache::instance().add(params, relxill_spec);
    delete rel_param;
    delete xill_param;
  }

  rebin_and_normalize_relxill_for_xspec(xspec_spectrum, relxill_spec);

  if (is_debug_run()) {
    auto time_elapsed_msec = std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::steady_clock::now() - tstart).count();
    printf("Time elapsed for model evaluation:      %.3fmsec\n", static_cast<double>(time_elapsed_msec));
  }

}


void assert_energy_grids_identical(const RelxillSpec &spectrum, const relline_spec_multizone *rel_profile) {
  assert(rel_profile->n_ener == static_cast<int>(spectrum.num_flux_bins));
  assert(rel_profile->ener[0] == spectrum.energy()[0]);
  assert(rel_profile->ener[1] == spectrum.energy()[1]);
}

/***
 * @brief main convolution routine for relxill, on a standard, global energy grid
 * as defined my RelxillSpec
 * @param relxill_spec
 * @param rel_profile
 * @param xill_spec_zones
 * @param spec_cache
 * @param caching_status
 * @param status
 */
void relxill_convolution_multizone(const RelxillSpec &relxill_spec,
                                   const relline_spec_multizone *rel_profile,
                                   const SpectrumZones &xill_spec_zones,
                                   specCache *spec_cache,
                                   const CachingStatus &caching_status,
                                   int *status) {

  CHECK_STATUS_VOID(*status);

  assert_energy_grids_identical(relxill_spec, rel_profile);

  // by this approach, all 3 spectra have the same energy grid
  auto xill_rebinned_spec = RelxillSpec();
  auto conv_out_spec = RelxillSpec();

  // make sure the output array is set to 0
  for (size_t ie = 0; ie < relxill_spec.num_flux_bins; ie++) {
    relxill_spec.flux[ie] = 0.0;
  }

  for (int izone = 0; izone < rel_profile->n_zones; izone++) { /***** loop over ionization zones   ******/

    /** avoid problems where no relxill bin falls into an ionization bin **/
    if (calcSum(rel_profile->flux[izone], rel_profile->n_ener) < 1e-12) {
      continue;
    }

    _rebin_spectrum(xill_rebinned_spec.energy(), xill_rebinned_spec.flux, xill_rebinned_spec.num_flux_bins,
                    xill_spec_zones.energy(), xill_spec_zones.flux[izone], xill_spec_zones.num_flux_bins);

    // --2-- convolve the relxill_spec on the energy grid "ener_conv" **
    //       note, always recompute fft for xillver if any parameter changed, as m_cache_relat changes the angular distribution
    int recompute_xill = caching_status.any_parameter_changed();
    convolveSpectrumFFTNormalized(xill_rebinned_spec.energy(), xill_rebinned_spec.flux, rel_profile->flux[izone],
                                  conv_out_spec.flux, conv_out_spec.num_flux_bins,
                                  caching_status.recomput_relat(), recompute_xill, izone, spec_cache, status);
    CHECK_STATUS_VOID(*status);

    // --3-- add it to the final output relxill_spec
    for (size_t ie = 0; ie < relxill_spec.num_flux_bins; ie++) {
      relxill_spec.flux[ie] += conv_out_spec.flux[ie];
    }

    if (is_debug_run() && rel_profile->n_zones <= 10) {
      write_output_spec_zones(relxill_spec, izone, status);
    }


  } /**** END OF LOOP OVER RADIAL ZONES *****/

}
