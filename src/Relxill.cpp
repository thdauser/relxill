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
#include "Relxill.h"
#include "IonGradient.h"
#include "XspecSpectrum.h"

/** caching parameters **/
relParam *cached_rel_param = nullptr;
xillParam *cached_xill_param = nullptr;

///////////////////////////////////////
// Forward Definitions of Functions  //
///////////////////////////////////////
double calculate_ecut_on_disk(const relParam *rel_param, double ecut_primary,
                              const double *rgrid, int num_zones, const IonGradient &ion_gradient, int ii);

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

void copy_spectrum_to_cache(const XspecSpectrum &spectrum,
                            specCache *spec_cache,
                            int *status) {/** initialize the cached output spec array **/
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
                    xillParam *xill_param,
                    relParam *rel_param,
                    int *status) {


  // TODO: this should be done when loading the input parameters
  // in case of an ionization gradient, we need to update the number of zones
  if (is_iongrad_model(xill_param->ion_grad_type)) {
    rel_param->num_zones = get_num_zones(rel_param->model_type, rel_param->emis_type, xill_param->ion_grad_type);
  }

  double ecut_primary = calc_ecut_at_primary_source(xill_param, rel_param, xill_param->ect, status);
  CHECK_STATUS_VOID(*status);

  specCache *spec_cache = init_global_specCache(status);
  assert(spec_cache != nullptr);
  CHECK_STATUS_VOID(*status);

  auto caching_status = CachingStatus();
  check_caching_parameters(caching_status, rel_param, xill_param);
  check_caching_energy_grid(caching_status, spec_cache, spectrum);

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
    double *rgrid = get_rzone_grid(rel_param->rin, rel_param->rout, rel_param->num_zones, rel_param->height, status);
    CHECK_STATUS_VOID(*status);

    // --- 1 --- calculate disk irradiation
    RelSysPar *sysPar = get_system_parameters(rel_param, status);
    CHECK_STATUS_VOID(*status);


    // --- 2 --- calculate ionization gradient
    IonGradient ion_gradient{rgrid, rel_param->num_zones, xill_param->ion_grad_type};
    ion_gradient.calculate(rel_param, xill_param);

    xillTableParam* xill_table_param[rel_param->num_zones];

    // --- 3 --- calculate xillver reflection spectra  (for every zone)
    for (int ii = 0; ii < rel_param->num_zones; ii++) {

      xill_table_param[ii] = get_xilltab_param(xill_param, status);

      // set xillver parameters for the given zone
      xill_table_param[ii]->ect =
          calculate_ecut_on_disk(rel_param, ecut_primary, rgrid, rel_param->num_zones, ion_gradient, ii);
      xill_table_param[ii]->lxi = ion_gradient.lxi[ii];
      xill_table_param[ii]->dens = ion_gradient.dens[ii];


      // --- 3a: load xillver spectra
      if (caching_status.xill == cached::no) {
        //  - always need to re-compute for an ionization gradient, TODO: can we do better caching?
        if (spec_cache->xill_spec[ii] != nullptr) {
          free_xill_spec(spec_cache->xill_spec[ii]);
        }
        spec_cache->xill_spec[ii] = get_xillver_spectra_table(xill_table_param[ii], status);
      }

    }

    // -- 4 -- returning radiation
    for (int ii = 0; ii < rel_param->num_zones; ii++) {
       if (rel_param->return_rad > 0) {
        get_xillver_fluxcorrection_factors(spec_cache->xill_spec[ii], &(rel_param->return_rad_flux_correction_factor),
                                           &(rel_param->xillver_gshift_corr_fac),
                                           xill_table_param[ii], status);
        CHECK_STATUS_VOID(*status);
      }

    }

    // --- 5 --- calculate multi-zone relline profile
    xillTable *xill_tab = nullptr; // needed for the relline_profile call
    get_init_xillver_table(&xill_tab, get_xilltab_param(xill_param, status), status);
//    CHECK_STATUS_VOID(*status);

    int n_ener_conv; // energy grid for the convolution, only created
    double *ener_conv;
    get_relxill_conv_energy_grid(&n_ener_conv, &ener_conv, status);
    // depends on returning radiation (i.e., the xillver correction factors)
    relline_spec_multizone *rel_profile = relbase_multizone(ener_conv, n_ener_conv, rel_param, xill_tab, rgrid,
                                                            rel_param->num_zones, status);
 //   CHECK_STATUS_VOID(*status);


    relxill_convolution_multizone(spectrum, rel_profile, spec_cache, rel_param, caching_status, status);

    copy_spectrum_to_cache(spectrum, spec_cache, status);
    CHECK_STATUS_VOID(*status);

  }

  add_primary_component(spectrum.energy, spectrum.num_flux_bins(), spectrum.flux, rel_param, xill_param, status);
}



void relxill_convolution_multizone(const XspecSpectrum &spectrum,
                                   const relline_spec_multizone *rel_profile,
                                   specCache *spec_cache,
                                   const relParam *rel_param,
                                   const CachingStatus &caching_status,
                                   int *status) {

  CHECK_STATUS_VOID(*status);

  int n_ener_conv = rel_profile->n_ener;
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
    return ecut_input * (1 + grav_redshift(rel_param));
  } else if (xill_param->prim_type == PRIM_SPEC_NTHCOMP || xill_param->prim_type == PRIM_SPEC_BB) {
    return ecut_input;
  } else {
    RELXILL_ERROR("unknown primary spectrum type, can not calculate Ecut at the primary source", status);
    return -1;
  }
}

/** @brief get the energy shift in order to calculate the proper Ecut value on the disk (if nzones>1)
 *        - the cutoff is calculated for the (linear) middle of the radial zone
 *        - for nzones=1 the value at the primary source is returned
 * @param rel_param
 * @param ecut0
 * @param ecut_primary
 * @param rgrid
 * @param num_zones
 * @param ion_gradient
 * @param ii
 * @return
 */
double calculate_ecut_on_disk(const relParam *rel_param, double ecut_primary,
                              const double *rgrid, int num_zones, const IonGradient &ion_gradient, int ii) {

  if (num_zones == 1) {
    return ecut_primary; // TODO: not obvious what "ecut0" is (it is the input value, from the model fitting)
  } else {
    double rzone = 0.5 * (rgrid[ii] + rgrid[ii + 1]);
    double del_emit = (ion_gradient.del_emit == nullptr) ? 0.0
                                                         : ion_gradient.del_emit[ii];  // only relevant if beta!=0 (doppler boosting)
    return ecut_primary * gi_potential_lp(rzone, rel_param->a, rel_param->height, rel_param->beta, del_emit);
  }
}

