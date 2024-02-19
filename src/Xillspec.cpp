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

#include "Xillspec.h"
#include "Relphysics.h"

extern "C" {
#include "xilltable.h"
#include "writeOutfiles.h"
#include "common.h"
}

EnerGrid *global_xill_egrid_coarse = nullptr;
EnerGrid *global_energy_grid_relxill = nullptr;


/** @brief main routine for the xillver table: returns a spectrum for the given parameters
 *  @details
 *   - decides if the table needs to be initialized and/or more data loaded
 *   - automatically normalizes  the spectra to logN=1e15 and logXi=0
 *   - for relxill type models it returns spectra for all available inclinations
 * @param param
 * @param status
 * @return xillSpec
 */



EnerGrid *get_relxill_conv_energy_grid() {
  if (global_energy_grid_relxill == nullptr) {
    global_energy_grid_relxill = new EnerGrid;
    global_energy_grid_relxill->nbins = N_ENER_CONV;
    global_energy_grid_relxill->ener = new double[global_energy_grid_relxill->nbins + 1];
    get_log_grid(global_energy_grid_relxill->ener, global_energy_grid_relxill->nbins + 1,
                 EMIN_RELXILL_CONV, EMAX_RELXILL_CONV);
  }

  return global_energy_grid_relxill;
}



xillSpec *get_xillver_spectra_table(xillTableParam *param, int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  xillTable *tab = nullptr;
  const char *fname = get_init_xillver_table(&tab, param->model_type, param->prim_type, status);

  CHECK_STATUS_RET(*status, nullptr);
  assert(fname != nullptr);

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


/** @brief similar to get_xillver_spectra_table, but uses the full xill_param input (from xpsecc)
 *  @details see get_xillver_spectra_table for more details, this function is just a wrapper
 * @param param
 * @param status
 * @return xillSpec
 */
xillSpec *get_xillver_spectra(xillParam *param, int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  xillTableParam *param_table = get_xilltab_param(param, status);
  xillSpec *xill_spec = get_xillver_spectra_table(param_table, status);
  CHECK_STATUS_RET(*status, nullptr);

  free(param_table);

  return xill_spec;
}

/**
 * @brief calculates the angle averaged xillver spectrum from a xillSpec structure
 * @param xill_spec
 * @param status
 * @return (double*) spectrum
 */
double *calc_angle_averaged_xill_spec(const xillSpec *xill_spec, int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  auto *output_spec = new double[xill_spec->n_ener];
  CHECK_MALLOC_RET_STATUS(output_spec, status, nullptr)

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



EnerGrid *get_coarse_xillver_energrid() {

  if (global_xill_egrid_coarse == nullptr) {
    global_xill_egrid_coarse = new EnerGrid ; // new_EnerGrid(status);
    global_xill_egrid_coarse->nbins = N_ENER_COARSE;
    global_xill_egrid_coarse->ener = new double[global_xill_egrid_coarse->nbins + 1];
    get_log_grid(global_xill_egrid_coarse->ener, global_xill_egrid_coarse->nbins + 1,
                 EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION);
  }

  return global_xill_egrid_coarse;
}

double calcNormWrtXillverTableSpec(const double *flux, const double *ener, const int n, int *status) {
  /* get the normalization of the spectrum with respect to xillver
 *  - everything is normalized using dens=10^15 cm^3
 *  - normalization defined, e.g., in Appendix of Dauser+2016
 *  - needs to be calculated on the specific energy grid (defined globally)
 */
  double keV2erg = 1.602177e-09;

  // need this to make sure no floating point problems arise
  double floatCompFactor = 1e-6;
  if (ener[n] < (EMAX_XILLVER_NORMALIZATION - floatCompFactor)
      || ener[0] > (EMIN_XILLVER_NORMALIZATION + floatCompFactor)) {
    RELXILL_ERROR("can not calculate the primary spectrum normalization", status);
    printf("  the given energy grid from %e keV to %e keV does not cover the boundaries\n", ener[0], ener[n]);
    printf("  from [%e,%e] necessary for the calculation\n", EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION);
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

/**
 * @brief cutoff pl model in xpsec units
 * @param pl_flux_xill [return] photon flux in counts/bin
 * @param ener
 * @param n_ener
 * @param xill_param
 * @param ener_shift: spectrum shifted in energy by this value, normalization is kept constant (i.e, norm not shifted by g^Gamma)
 */
void spec_cutoffpl(double *pl_flux_xill,
                   const double *ener,
                   int n_ener,
                   double ect,
                   double gamma,
                   double ener_shift = 1.0) {
  /** note that in case of the nthcomp model Ecut is in the frame of the primary source
      but for the bkn_powerlaw it is given in the observer frame */

  double ecut = ect * ener_shift; // Important: Ecut is given in the frame of the observer

  for (int ii = 0; ii < n_ener; ii++) {
    double en = 0.5 * (ener[ii] + ener[ii + 1]);
    pl_flux_xill[ii] = exp(1.0 / ecut) *
        pow(en, -gamma) *
        exp(-en / ecut) *
        (ener[ii + 1] - ener[ii]);
  }
}

/**
* @brief nthcomp model in xpsec units
* @param pl_flux_xill [return] photon flux in counts/bin
* @param ener
* @param n_ener
* @param xill_param
* @param ener_shift  spectrum shifted in energy by this value, normalization is kept constant! (i.e, norm not shifted by g^Gamma)
*/
void spec_nthcomp(double *pl_flux_xill,
                  const double *ener,
                  int n_ener,
                  const xillTableParam *xill_param,
                  double ener_shift = 1.0) {
  double nthcomp_param[5];
  double kTe = xill_param->ect; // Important: kTe is given in the frame of the source
  double z = 1 / ener_shift - 1; // convert energy shift to redshift
  get_nthcomp_param(nthcomp_param, xill_param->gam, kTe, z);
  c_donthcomp(ener, n_ener, nthcomp_param, pl_flux_xill);
}

/**
* @brief blackbody model in xpsec units
* @param pl_flux_xill [return] photon flux in counts/bin
* @param ener
* @param n_ener
* @param xill_param
*/
void spec_blackbody(double *pl_flux_xill, const double *ener, int n_ener, const xillTableParam *xill_param) {
  double en;
  for (int ii = 0; ii < n_ener; ii++) {
    en = 0.5 * (ener[ii] + ener[ii + 1]);
    pl_flux_xill[ii] = en * en / (pow(xill_param->kTbb, 4) * (exp(en / xill_param->kTbb) - 1));
    pl_flux_xill[ii] *= (ener[ii + 1] - ener[ii]);
  }
}

/**
 * @brief calculate the primary spectrum, which can be shifted by the factor ener_shift_source_obs
 * @param pl_flux_xill: photon flux in xspec units (cts/bin)
 * @param ener_shift_source_obs: energy shift applied the spectrum
 */
void calc_primary_spectrum(double *pl_flux_xill,
                           const double *ener,
                           int n_ener,
                           const xillTableParam *xill_param,
                           int *status,
                           double ener_shift_source_obs) {

  CHECK_STATUS_VOID(*status);
  assert(ener_shift_source_obs >= 0.0);

  // ecut and kTe are given in the frame of the source
  if (xill_param->prim_type == PRIM_SPEC_ECUT) {
    spec_cutoffpl(pl_flux_xill, ener, n_ener, xill_param->ect, xill_param->gam, ener_shift_source_obs);

  } else if (xill_param->prim_type == PRIM_SPEC_NTHCOMP) {
    spec_nthcomp(pl_flux_xill, ener, n_ener, xill_param, ener_shift_source_obs);

  } else if (xill_param->prim_type == PRIM_SPEC_BB) {
    spec_blackbody(pl_flux_xill, ener, n_ener, xill_param);

  } else {
    RELXILL_ERROR("trying to calculate a primary continuum for an undefined primary spectral type (should not happen!)",
                  status);
  }
}

double *calc_normalized_primary_spectrum(const double *ener, int n_ener,
                                         const relParam *rel_param, const xillTableParam *xill_param, int *status) {

  // need to create a specific energy grid for the primary component normalization to
  // fulfill the XILLVER NORM condition (Dauser+2016)
  EnerGrid *egrid = get_coarse_xillver_energrid();
  auto prim_spec_source = new double[egrid->nbins];
  calc_primary_spectrum(prim_spec_source, egrid->ener, egrid->nbins, xill_param, status);

  // calculate the normalization at the primary source (i.e., not shifted by the energy shift)
  double const
      norm_factor_prim_spec = 1. / calcNormWrtXillverTableSpec(prim_spec_source, egrid->ener, egrid->nbins, status);
  delete[] prim_spec_source;

  // if we have the LP geometry defined, the spectrum will be different between the primary source and the observer
  // due to energy shift: need to calculate the normalization for the source, but the spectrum at the observer
  double energy_shift = 1.0;
  if (rel_param != nullptr && rel_param->emis_type == EMIS_TYPE_LP) {
    energy_shift = energy_shift_source_obs(rel_param);
  }
  // output flux
  auto prim_spec_observer = new double[n_ener];
  calc_primary_spectrum(prim_spec_observer, ener, n_ener, xill_param, status, energy_shift);

  // take the normalization factor into account
  for (int ii = 0; ii < n_ener; ii++) {
    prim_spec_observer[ii] *= norm_factor_prim_spec;
  }

  return prim_spec_observer;
}

/**
 * @brief calculate the flux ration between the xillver reflection spectrum and its
 * input spectrum in the 0.1-1000 keV energy range
 * @param flu
 * @param ener
 * @param n_ener
 * @param xill_table_param
 * @param status
 * @return
 */
double get_xillver_fluxcorr(double *flu, const double *ener, int n_ener,
                            const xillTableParam *xill_table_param, int *status) {

  // use a coarser energy grid to evaluate the direct spectrum, as we are only interested
  // in the total flux in the band
  auto egrid_coarse = get_coarse_xillver_energrid();
  auto direct_spec =
      calc_normalized_primary_spectrum(egrid_coarse->ener, egrid_coarse->nbins,
                                       nullptr, xill_table_param, status);

  double const emin = 0.1;
  double const emax = 1000;
  double const ratio_refl_direct =
      get_energy_flux_band(ener, flu, n_ener, emin, emax) /
          get_energy_flux_band(egrid_coarse->ener, direct_spec, egrid_coarse->nbins, emin, emax);

  delete[] direct_spec;
  return ratio_refl_direct;
}

/**
 * @brief calculate the normalization factor to convert the initial xillver primary spectrum
 * (given by xill_param_0) under the given energy shift, by taking the fixed xillver incident
 * spectrum normalization in the 0.1-1000keV range into account (Dauser+16, A1). Note that the
 * normalization of the primary spectrum therefore scales as 1/norm_flux_cgs.
 * @param energy_shift_source_disk [energy shift from the source to the disk]
 * @param xill_param_0 [xillver parameters for the inital spectrum]
 * @return norm_flux_cgs
 */
double calc_xillver_normalization_change(double energy_shift, const xillTableParam *xill_param_0) {

  int status = EXIT_SUCCESS;
  EnerGrid *egrid = get_coarse_xillver_energrid();
  auto prime_spec = new double[egrid->nbins];

  calc_primary_spectrum(prime_spec, egrid->ener, egrid->nbins, xill_param_0, &status);
  double source_spec_norm_factor = 1. / calcNormWrtXillverTableSpec(prime_spec, egrid->ener, egrid->nbins, &status);
  // normalized source spec: prime_source_source*source_spec_norm_factor

  // get parameters for the source spectrum (ecut shifted in energy)
  auto xill_param_source = new xillTableParam;
  *xill_param_source = *xill_param_0; // shallow copy
  xill_param_source->ect = xill_param_0->ect * energy_shift;

  // calculate the primary spectrum at the source
  calc_primary_spectrum(prime_spec, egrid->ener, egrid->nbins, xill_param_source, &status);
  double disk_spec_norm_factor = 1. / calcNormWrtXillverTableSpec(prime_spec, egrid->ener, egrid->nbins, &status);
  // normalized prime spec: prime_spec_disk*disk_spec_norm_factor

  delete xill_param_source;
  delete[] prime_spec;

  return disk_spec_norm_factor / source_spec_norm_factor;
}

/**
 * @brief calculate the normalization factor to convert the initial xillver primary spectrum
 * (given by xill_param_0) under the given energy shift, by taking the fixed xillver incident
 * spectrum normalization in the 0.1-1000keV range into account (Dauser+16, A1).
 * @param energy_shift_source_disk [energy shift from the source to the disk for each zone]
 * @param n_zones [number of zones]
 * @param xill_param_0 [xillver parameters for the inital spectrum]
 * @return norm_flux_cgs
 */
double *calc_xillver_normalization_change_source_to_disk(const double *energy_shift,
                                                         int n_zones,
                                                         const xillTableParam *xill_param_0) {

  int status = EXIT_SUCCESS;
  EnerGrid *egrid = get_coarse_xillver_energrid();
  auto prime_spec = new double[egrid->nbins];

  calc_primary_spectrum(prime_spec, egrid->ener, egrid->nbins, xill_param_0, &status);
  double source_spec_norm_factor = 1. / calcNormWrtXillverTableSpec(prime_spec, egrid->ener, egrid->nbins, &status);
  // normalized disk spec: prime_spec_source*source_spec_norm_factor

  auto norm_factors = new double[n_zones];

  // get parameters for the source spectrum (ecut shifted in energy)
  auto xill_param_disk = new xillTableParam;
  *xill_param_disk = *xill_param_0; // shallow copy
  for (int ii = 0; ii < n_zones; ii++) {
    xill_param_disk->ect = xill_param_0->ect * energy_shift[ii];

    // calculate the primary spectrum at the source
    calc_primary_spectrum(prime_spec, egrid->ener, egrid->nbins, xill_param_disk, &status);
    double disk_spec_norm_factor = 1. / calcNormWrtXillverTableSpec(prime_spec, egrid->ener, egrid->nbins, &status);
    // normalized prime spec: prime_spec_source*source_spec_norm_factor

    norm_factors[ii] = disk_spec_norm_factor / source_spec_norm_factor;
  }

  delete xill_param_disk;
  delete[] prime_spec;

  return norm_factors;
}

/**
 * @brief calculate the gshift flux correction factor
 * @detail The gshift-flux-correction factor is defined as the ratio of the flux-change of the xillver spectrum
 * by shifting it to a factor 1.5 lower energies wrt to the flux change of a powerlaw-spectrum (with the same
 * gamma), shifted by the same amoount.
 *
 * C_g = ( F(g=g_r)/F(g=g_0)  )  /  ( (g_r)^Gamma / (g_0)^Gamma )
 *     = ( F(g=1.5)/F(g=1) ) /  ( (1.5)^Gamma
 * @param  flux (double*): xillver spectrum
 * @param  energy (double*): energy grid on which the spectrum is defined
 * @param  n_energy (int)
 * @param  gamma (photon index)
 * @return fac_gshift_fluxcorr: ratio of xillver spectrum shifted to lower energy by a factor 1.5 (g=2/3) to a
 * shifted powerlaw with gamma (gamma is the input to the xillver spec)
 */
double get_xillver_gshift_fluxcorr(double *flu, const double *ener, int n_ener, double gamma) {

  const double gshift_refvalue = 2. / 3.;  // shift it 1.5 to lower energies

  auto ener_z = new double[n_ener + 1];
  for (int ii = 0; ii <= n_ener; ii++) {  // ener array has n_ener+1 entries
    ener_z[ii] = ener[ii] / gshift_refvalue;
  }

  auto flu_z = new double[n_ener];

  rebin_spectrum(ener_z, flu_z, n_ener, ener, flu, n_ener);
  for (int ii = 0; ii < n_ener; ii++) {
    flu_z[ii] *= gshift_refvalue;  // take time dilation into account (dE already taken into account as bin-integ)
  }

  delete[] ener_z;

  double emin = 0.15;
  double emax = 500.0;

  // ratio of no-shift wrt to a shift if 1.5 to lower energies, i.e.
  // this is equivalent to F(g=1.5)/F(g=1)
  double gshift_ratio =
      get_photon_flux_band(ener, flu, n_ener, emin, emax) /
          get_photon_flux_band(ener, flu_z, n_ener, emin, emax);

  delete[] flu_z;

  return gshift_ratio / pow(1.5, gamma);
}

/**
 * @brief calculate the (energy) flux correction factors for a xillver spectrum.
 * @detail Two correction factors are calculate, which are both the ratio of the energy flux of two spectra:
 * (1) fac_fluxcorr: ratio of reflected xillver spectrum to its input spectrum (below 1 for low ionization and above 1
 * for high ionization)
 * (2) fac_gshift_fluxcorr: ratio of xillver spectrum shifted to lower energy by a factor 1.5 (g=2/3) to a
 * shifted powerlaw with gamma (gamma is the input to the xillver spec)
 * @param (output) fac_fluxcorr (double*)
 * @param (output) fac_gshift_fluxcorr (double*)
 * @param xill_param
 */
void get_xillver_fluxcorrection_factors(const xillSpec *xill_spec,
                                        double *fac_fluxcorr,
                                        double *fac_gshift_fluxcorr,
                                        xillTableParam *xill_table_param,
                                        int *status) {
  CHECK_STATUS_VOID(*status);

  double *angle_averaged_xill_spec = calc_angle_averaged_xill_spec(xill_spec, status);
  CHECK_STATUS_VOID(*status);

  if (fac_fluxcorr != nullptr) {
    *fac_fluxcorr =
        get_xillver_fluxcorr(angle_averaged_xill_spec, xill_spec->ener, xill_spec->n_ener, xill_table_param, status);
  }

  if (fac_gshift_fluxcorr != nullptr) {
    *fac_gshift_fluxcorr =
        get_xillver_gshift_fluxcorr(angle_averaged_xill_spec, xill_spec->ener, xill_spec->n_ener, xill_table_param->gam);
  }

  delete[] angle_averaged_xill_spec;

  if (*status != EXIT_SUCCESS) {
    RELXILL_ERROR("failed to calculate the flux correction factor", status);
  }

}

double norm_factor_semi_infinite_slab(double incl_deg) {
  return 0.5 * cos(incl_deg * M_PI / 180);
}

/**
 * @description: adds the proper flux normalization for a semi-infinate slab
 *  under inclination angle incl
 **/
void norm_xillver_spec(xillSpec *spec, double incl) {

  assert(spec->n_incl == 1); // make sure there is only one spectrum given (for the chosen inclination)

  double norm_factor = norm_factor_semi_infinite_slab(incl);
  for (int ii = 0; ii < spec->n_ener; ii++) {
    spec->flu[0][ii] *= norm_factor;
  }

}




/**
 * @Function: calc_xillver_angdep
 * @Synopsis: Calculate the Angle Weighted Xillver Spectrum on the Standard Relxill Spectrum
 * @Input: xill_spec
 *         rel_dist[n_incl]
 * @Output: xill_flux  (needs to be allocated, will be overwritten)
 **/
void calc_xillver_angdep(double *xill_flux, xillSpec *xill_spec, const double *dist, const int *status) {

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
                             double *ener,
                             int n_ener,
                             double *rel_dist,
                             xillSpec *xill_spec,
                             int *status) {

  auto xill_angdist_inp = new double[xill_spec->n_ener];

  calc_xillver_angdep(xill_angdist_inp, xill_spec, rel_dist, status);

  rebin_spectrum(ener, o_xill_flux, n_ener,
                 xill_spec->ener, xill_angdist_inp, xill_spec->n_ener);

  delete[] xill_angdist_inp;

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
  get_xillver_angdep_spec(xill_flux, ener, n_ener, rel_cosne_dist, xillSpecTable, status);

}

