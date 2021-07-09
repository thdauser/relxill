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


/** caching parameters **/
relParam *cached_rel_param = NULL;
xillParam *cached_xill_param = NULL;




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
  if ((cached_rel_param != nullptr) && (shouldOutfilesBeWritten() == 0)) {

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

  /** did any of the parameters change?  **/
  *re_rel = redo_relbase_calc(rel_param, cached_rel_param);
  *re_xill = redo_xillver_calc(rel_param, xill_param, cached_rel_param, cached_xill_param);
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

  xillTable *xill_tab = nullptr;
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

  specCache* spec_cache = init_global_specCache(status);
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


    // set Return Rad correction factor before checking the caching
    // (xillver parameters could influence it)
    if (rel_param->return_rad > 0) {
      get_xillver_fluxcorrection_factors(&(rel_param->return_rad_flux_correction_factor),
                                         &(rel_param->xillver_gshift_corr_fac),
                                         xill_param, status);
      if (*status != EXIT_SUCCESS) {
        RELXILL_ERROR("failed to calculate the flux correction factor", status);
      }
    }


    // stored the parameters for which we are calculating
    set_cached_xill_param(xill_param, &cached_xill_param, status);
    set_cached_rel_param(rel_param, &cached_rel_param, status);

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
        xill_param->ect = ecut_primary * gi_potential_lp(rzone, rel_param->a, rel_param->height, rel_param->beta, del_emit);
        //
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
      convolveSpectrumFFTNormalized(ener, xill_flux, rel_profile->flux[ii], conv_out, n_ener,
                        recompute_rel, 1, ii, spec_cache, status);
      CHECK_STATUS_VOID(*status);

      // rebin to the output grid
      rebin_spectrum(ener_inp, single_spec_inp, n_ener_inp, ener, conv_out, n_ener);

      // add it to the final output spectrum
      for (int jj = 0; jj < n_ener_inp; jj++) {
        spec_inp[jj] += single_spec_inp[jj] ;
      }

      if (shouldOutfilesBeWritten() && rel_profile->n_zones <= 10) {
        char vstr[200];
        double test_flu[n_ener_inp];
        for (int jj = 0; jj < n_ener_inp; jj++) {
          test_flu[jj] = single_spec_inp[jj];
        }
        if (sprintf(vstr, "test_relxill_spec_zones_%03i.dat", ii + 1) == -1) {
          RELXILL_ERROR("failed to get filename", status);
        }
        save_xillver_spectrum(ener_inp, test_flu, n_ener_inp, vstr);

//        save_xillver_spectrum(ener, xill_flux, n_ener, "test_fft_xill.dat" );
//        save_xillver_spectrum(ener,rel_profile->flux[ii], n_ener, "test_fft_rel.dat");
//        save_xillver_spectrum(ener, conv_out, n_ener, "test_fft_conv.dat");

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
