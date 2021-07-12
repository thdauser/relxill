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

#include "xillspec.h"

#include "xilltable.h"
#include "writeOutfiles.h"
#include "relbase.h"

/**
 * @description: the main routine for the xillver table: returns a spectrum for the given parameters
 *  - decides if the table needs to be initialized and/or more data loaded
 *  - automatically normalizes  the spectra to logN=1e15 and logXi=0
 *  - for relxill type models it returns spectra for all available inclinations
 * */
xillSpec *get_xillver_spectra(xillParam *param, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  xillTable *tab = NULL;
  const char *fname = get_init_xillver_table(&tab, param, status);

  CHECK_STATUS_RET(*status, NULL);
  assert(fname != NULL);

  // =1=  get the inidices
  int *indparam = get_xilltab_indices_for_paramvals(param, tab, status);

  // =2=  check if the necessary spectra for interpolation are loaded
  check_xilltab_cache(fname, param, tab, indparam, status);

  // =3= interpolate values
  xillSpec *spec = interp_xill_table(tab, param, indparam, status);

  CHECK_RELXILL_DEFAULT_ERROR(status);

  free(indparam);
  return spec;
}

/**
 * @brief calculates the angle averaged xillver spectrum from a xillSpec structure
 * @param xill_spec
 * @param status
 * @return (double*) spectrum
 */
double *calc_angle_averaged_xill_spec(xillSpec *xill_spec, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  double *output_spec = malloc(sizeof(double) * xill_spec->n_ener);
  CHECK_MALLOC_RET_STATUS(output_spec, status, NULL)

  memset(output_spec, 0, xill_spec->n_ener * sizeof(output_spec[0]));
  for (int ii = 0; ii < xill_spec->n_incl; ii++) {

    double incl_factor = norm_factor_semi_infinite_slab(xill_spec->incl[ii] * 180 / M_PI) / xill_spec->n_incl;
    for (int jj = 0; jj < xill_spec->n_ener; jj++) {
      output_spec[jj] += incl_factor * xill_spec->flu[ii][jj];
    }
  }

  return output_spec;
}

/**
 * @brief calculate the ENERGY flux in a given band from a given bin-integrated photon flux (cts/bin), the standard unit to store
 * all spectra in the relxill code
 **/
static double get_energy_flux_band(const double *ener,
                                   const double *photon_flux,
                                   int n_ener,
                                   double emin,
                                   double emax) {
  double sum = 0.0;
  for (int ii = 0; ii < n_ener; ii++) {
    if (ener[ii] >= emin && ener[ii + 1] <= emax) {
      sum += photon_flux[ii] * 0.5 * (ener[ii] + ener[ii + 1]);
    }
  }
  return sum;
}

/**
 * @brief calculate the PHOTON flux in a given band from a given bin-integrated photon flux (cts/bin), the standard unit to store
 * all spectra in the relxill code
 **/
static double get_photon_flux_band(const double *ener,
                                   const double *photon_flux,
                                   int n_ener,
                                   double emin,
                                   double emax) {
  double sum = 0.0;
  for (int ii = 0; ii < n_ener; ii++) {
    if (ener[ii] >= emin && ener[ii + 1] <= emax) {
      sum += photon_flux[ii];
    }
  }
  return sum;
}

double get_xillver_fluxcorr(double *flu, const double *ener, int n_ener,
                            const xillParam *xill_param, int *status) {
  double *direct_spec =
      calc_normalized_xillver_primary_spectrum(ener, n_ener, NULL, xill_param, status);

  if (shouldOutfilesBeWritten()) {
    save_xillver_spectrum(ener, direct_spec, n_ener, "test-debug-xillver-direct.dat");
    save_xillver_spectrum(ener, flu, n_ener, "test-debug-xillver-refl.dat");
  }

  double emin = 0.1;
  double emax = 1000;
  double ratio_refl_direct =
      get_energy_flux_band(ener, flu, n_ener, emin, emax) /
          get_energy_flux_band(ener, direct_spec, n_ener, emin, emax);

  free(direct_spec);
  return ratio_refl_direct;
}

double get_xillver_gshift_fluxcorr(double *flu, const double *ener, int n_ener, double gamma) {

  const double gshift_refvalue = 2. / 3.;  // shift it 1.5 to lower energies

  double ener_z[n_ener + 1];
  for (int ii = 0; ii <= n_ener; ii++) {  // ener array has n_ener+1 entries
    ener_z[ii] = ener[ii] / gshift_refvalue;
  }

  double flu_z[n_ener];

  rebin_spectrum(ener_z, flu_z, n_ener, ener, flu, n_ener);
  for (int ii = 0; ii < n_ener; ii++) {
    flu_z[ii] *= gshift_refvalue;  // take time dillation into account (dE already taken into account as bin-integ)
  }

  if (shouldOutfilesBeWritten()) {
    save_xillver_spectrum(ener, flu, n_ener, "test-debug-xillver-gshift-0.dat");
    save_xillver_spectrum(ener, flu_z, n_ener, "test-debug-xillver-gshift-z.dat");
  }

  double emin = 0.15;
  double emax = 500.0;

  // ratio of no-shift wrt to a shift if 1.5 to lower energies
  double gshift_ratio =
      get_photon_flux_band(ener, flu, n_ener, emin, emax) /
          get_photon_flux_band(ener, flu_z, n_ener, emin, emax);

  return gshift_ratio / pow(1.5, gamma);
}

/**
 * @brief calculate the (energy) flux correction factors for a xillver spectrum.
 * @detail Two correction factors are calculate, which are both the ratio of the energy flux of two spectra:
 * (1) fac_fluxcorr: ratio of reflected xillver spectrum to its input spectrum (below 1 for low ionization and above 1
 * for high ionization)
 * (2) fac_gshift_fluxcorr: ratio of xillver spectrum shifted to lower energy by a factor 1.5 (g=2/3) to a
 * non-shifted xillver spetctrum
 * @param (output) fac_fluxcorr (double*)
 * @param (output) fac_gshift_fluxcorr (double*)
 * @param xill_param
 */
void get_xillver_fluxcorrection_factors(double *fac_fluxcorr, double *fac_gshift_fluxcorr,
                                        xillParam *xill_param, int *status) {
  CHECK_STATUS_VOID(*status);

  xillSpec *xill_spec = get_xillver_spectra(xill_param, status);
  assert(xill_spec->n_incl > 1);  // has to be the case for the inclination interpolated xillver model

  double *angle_averaged_xill_spec = calc_angle_averaged_xill_spec(xill_spec, status);
  CHECK_STATUS_VOID(*status);

  if (fac_fluxcorr != NULL) {
    *fac_fluxcorr =
        get_xillver_fluxcorr(angle_averaged_xill_spec, xill_spec->ener, xill_spec->n_ener, xill_param, status);
  }

  if (fac_gshift_fluxcorr != NULL) {
    *fac_gshift_fluxcorr =
        get_xillver_gshift_fluxcorr(angle_averaged_xill_spec, xill_spec->ener, xill_spec->n_ener, xill_param->gam);
  }

  free_xill_spec(xill_spec);
  free(angle_averaged_xill_spec);

}

double norm_factor_semi_infinite_slab(double incl_deg) {
  return 0.5 * cos(incl_deg * M_PI / 180);
}

/**
 * @description: adds the proper flux normalization for a semi-infinate slab
 *  under inclination angle incl
 **/
void norm_xillver_spec(xillSpec *spec, double incl) {
  double norm_factor = norm_factor_semi_infinite_slab(incl);
  for (int ii = 0; ii < spec->n_ener; ii++) {
    spec->flu[0][ii] *= norm_factor;
  }

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



/**
 * calculate a xillver spectrum and combine it according to the angle distribution
 * (stored in rel_cosne_dist)
 * @param [output] xill_flux  (needs to be allocated)
 * @param ener
 * @param n_ener
 * @param xill_param
 * @param rel_cosne_dist
 * @param status
 */
void getNormalizedXillverSpec(double* xill_flux, double* ener, int n_ener, xillParam* xill_param,
                              double* rel_cosne_dist, int* status ){

  CHECK_STATUS_VOID(*status);

  xillSpec* xillSpecTable = get_xillver_spectra(xill_param, status);
  get_xillver_angdep_spec(xill_flux, n_ener, ener, rel_cosne_dist, xillSpecTable, status);

}

