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

#include "Relreturn_Corona.h"
#include "Relbase.h"

extern "C" {
#include "relreturn_datastruct.h"
#include "relutility.h"
}

/**
 * @brief From the known ratio between xillver and power law flux boost, calculate the expected flux
 * correction / factor for any g value (certainly not very accurate for large values of g). The function
 * is assuring that we are not over-correcting and does ensure that for g>1 we never have a flux reduction
 * and similarly for g<1 never a flux boost.
 * @param xill_gshift_fac  : factor between flux boost from xillver for F(g=1.5)/F(g=1) wrt to g^Gamma for g=1.5
 * @param g : energy shift
 * @param gamma: power law index gamma
 * @return flux boost factor (would be g^gamma for a powerlaw instead of the xillver reflection)
 */
double corrected_gshift_fluxboost_factor(double xill_gshift_fac, double g, double gamma) {

  if (fabs(g - 1) <= 1e-3) { // without any significant energy shift, there is no need for a flux correction
    return 1.0;
  }

  if (fabs(xill_gshift_fac - 1) < 1e-3) { // if the correction is not significant, we don't need any corrections
    return pow(g, gamma);
  }

  double g0 = 2. / 3;  // 1/g, where g is used to calculate xill_gshift_fac

  double alin; double blin;
  double gshift_corr_factor;
  if (xill_gshift_fac < 1) { // parabola which has f(g=0)=0
    double a = (xill_gshift_fac / g0 - 1) / (g0 - 1);
    double b = 1 - a;

    gshift_corr_factor = (g >= 1)
                         ? 1. / g * (1. / g * a + b)
                         : g * (g * a + b);

  } else { // linear interpolation
    alin = (xill_gshift_fac - 1) / (g0 - 1);
    blin = 1 - alin;

    gshift_corr_factor = (g >= 1)
                         ? (1. / g * alin + blin)
                         : (g * alin + blin);
  }

  double fluxboost_factor = pow(g, gamma) * gshift_corr_factor;

  // ensure that we are not over-correcting (meaning that for g>1 we do not allow a flux reduction)

  if (g > 1)
    assert(fluxboost_factor >= 0.9);
  if (g < 1){
    if (fluxboost_factor>1){
      if (fluxboost_factor > 1.1 ) { // print warning for a deviation of >10%
        if (is_debug_run()) {
          printf(" *** warning: for g=%.4f the gshift-fluxboost factor = %.4f > 1 -> resetting to 1 \n",
                 g, fluxboost_factor);
        }
      }
      fluxboost_factor = 1;
    }
  }
  assert(fluxboost_factor > 0);

  return fluxboost_factor;
}


emisProfile* calc_rrad_emis_corona(const returningFractions *ret_fractions, rradCorrFactors* corr_factors,
                                   const emisProfile* emis_input, double gamma, int* status) {

  // make very rough sanity checks here
  // assert(gshift_corr_factor>1e-3);
  // assert(gshift_corr_factor<1e3);

  const int ng = ret_fractions->tabData->ng;
  const int nrad = ret_fractions->nrad;

  // need to have emis_input on the SAME radial zone grid, ASCENDING (as tables are ascending in radius)
  //  -> require the grid of the emissivity profile to be identical, means they share the same pointer
  assert(emis_input->re == ret_fractions->rad);
  assert(emis_input->nr == ret_fractions->nrad);
  assert(emis_input->re[0] < emis_input->re[1]);

  auto emis_single_zone = new double[nrad];
  auto gfac  = new double[ng];

  emisProfile* emis_return = new_emisProfile(ret_fractions->rad, ret_fractions->nrad, status); // ret_fractions->rad is not owned by emisReturn

  for (int i_rad_incident = 0; i_rad_incident < nrad; i_rad_incident++) {
    int itab_rad_incident = ret_fractions->irad[i_rad_incident];

    for (int i_rad_emitted = 0; i_rad_emitted < nrad; i_rad_emitted++) {
      int itab_rad_emitted = ret_fractions->irad[i_rad_emitted];

      double ratio_ut_obs2emit =
          ut_disk(ret_fractions->rad[i_rad_incident], ret_fractions->a)
              / ut_disk(ret_fractions->rad[i_rad_emitted], ret_fractions->a);

      get_gfac_grid(gfac, ret_fractions->tabData->gmin[itab_rad_incident][itab_rad_emitted],
                    ret_fractions->tabData->gmax[itab_rad_incident][itab_rad_emitted], ng);
      emis_single_zone[i_rad_emitted] = 0.0;
      for (int jj = 0; jj < ng; jj++) {
        // TODO: switch to new table, means removing this factor
        double corr_fac_new_deriv = ratio_ut_obs2emit / gfac[jj];
        assert(corr_fac_new_deriv > 0);

        double emis_g_zone =
            ret_fractions->tabData->frac_g[itab_rad_incident][itab_rad_emitted][jj] * corr_fac_new_deriv;

        if (corr_factors != nullptr) {  // TODO: fix this ugly nullptr-if-case
          emis_g_zone *= corrected_gshift_fluxboost_factor(corr_factors->corrfac_gshift[i_rad_emitted], gfac[jj], gamma);
        } else {
          emis_g_zone *= corrected_gshift_fluxboost_factor(1.0, gfac[jj], gamma);
        }

        emis_single_zone[i_rad_emitted] += emis_g_zone;

      }

      emis_single_zone[i_rad_emitted] *=
          ret_fractions->frac_i[i_rad_incident][i_rad_emitted]
              * emis_input->emis[i_rad_emitted]
              * ret_fractions->tabData->f_ret[ret_fractions->irad[i_rad_emitted]];

    }

    emis_return->emis[i_rad_incident] = calcSum(emis_single_zone, nrad);
   if (corr_factors != nullptr){
      emis_return->emis[i_rad_incident] *= corr_factors->corrfac_flux[i_rad_incident];
    }
  }

  delete[] emis_single_zone;
  delete[] gfac;

  return emis_return;
}

/**
 * @brief determine the lower and upper radius value of the emissivity profile
 * @param emisInput
 * @param rlo_emis [output]
 * @param rhi_emis [output]
 */
void determine_rlo_rhi(const emisProfile *emisInput, double *rlo_emis, double *rhi_emis) {
  if (is_emis_grid_ascending(emisInput)){
    (*rlo_emis) = emisInput->re[0];
    (*rhi_emis) = emisInput->re[emisInput->nr-1];
  } else {
    (*rlo_emis) = emisInput->re[emisInput->nr-1];
    (*rhi_emis) = emisInput->re[0];
  }
  assert((*rlo_emis) < (*rhi_emis));
}

rradCorrFactors *init_rrad_corr_factors(const double *rlo, const double *rhi, int n_zones) {
  auto* corr_factors = new rradCorrFactors;

  corr_factors->rgrid = new double[n_zones+1];
  for (int ii=0; ii<n_zones; ii++){
    corr_factors->rgrid[ii] = rlo[ii];
  }
  corr_factors->rgrid[n_zones] = rhi[n_zones-1];

  corr_factors->n_zones= n_zones;

  corr_factors->corrfac_flux = new double[n_zones];
  corr_factors->corrfac_gshift = new double[n_zones];

  return corr_factors;
}

void free_rrad_corr_factors(rradCorrFactors** p_corr_factors){
  if (*p_corr_factors != nullptr ) {
    delete[] (*p_corr_factors)->rgrid;
    delete[] (*p_corr_factors)->corrfac_flux;
    delete[] (*p_corr_factors)->corrfac_gshift;
    delete (*p_corr_factors);
    *p_corr_factors = nullptr;
  }
}

/** @brief assign binned values to new grid
 *  - input is the x (rgrid0) and y (value0), which is the bin_lo/bin_hi-grid with n0+1 values
 *  - use the radius (rmean), which is the mean value to determine in which bin of the new
 *    grid the value falls
 *  - no interpolation is used
 *  @param value[nzones]
 *  @param rmean[nzones]
 *  @param rgrid[n0+1]
 *  @param value0[n0-1]
 *
 */
static void rebin_to_grid(double* value, const double* rmean, const double nzones,
                          const double*rgrid0, const double* value0, const int n0){

  for (int ii=0; ii<nzones; ii++){
    int ind = binary_search(rgrid0, n0+1, rmean[ii]);
    if (rmean[ii] < rgrid0[0]) {
      ind = 0;
    } else if (rmean[ii] > rgrid0[n0]) {
      ind = n0;
    }
    value[ii] = value0[ind];
  }

}

rradCorrFactors* rebin_corrfactors_to_rradtable_grid
    (rradCorrFactors* input_corr_factors, returningFractions* ret_fractions, int* status) {

  if (input_corr_factors == nullptr){
    return nullptr;
  } else {

    rradCorrFactors* rtable_corr_factors =
        init_rrad_corr_factors(ret_fractions->rlo, ret_fractions->rhi, ret_fractions->nrad);

    rebin_to_grid(rtable_corr_factors->corrfac_flux, ret_fractions->rad, ret_fractions->nrad,
                  input_corr_factors->rgrid, input_corr_factors->corrfac_flux, input_corr_factors->n_zones);

    rebin_to_grid(rtable_corr_factors->corrfac_gshift, ret_fractions->rad, ret_fractions->nrad,
                  input_corr_factors->rgrid, input_corr_factors->corrfac_gshift, input_corr_factors->n_zones);

    return rtable_corr_factors;
  }


}

/**
 * main function to calculate the returning radiation emissivity profile.
 *
 * Caveat: The radial grid for the relline profile, which is typically used for the emissivity
 * profile, is defined with the radius DESCENDING
 */
emisProfile *get_rrad_emis_corona(const emisProfile* emis_input, const relParam* param, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  double rlo_emis, rhi_emis;
  determine_rlo_rhi(emis_input, &rlo_emis, &rhi_emis);

  returningFractions *ret_fractions = get_rrad_fractions(param->a, rlo_emis, rhi_emis, status);

  emisProfile* emis_input_rebinned = new_emisProfile(ret_fractions->rad, ret_fractions->nrad, status); // ret_fractions->rad is not owned by emis_return
  inv_rebin_mean(emis_input->re, emis_input->emis, emis_input->nr,
                 emis_input_rebinned->re, emis_input_rebinned->emis, emis_input_rebinned->nr, status);


  rradCorrFactors* rrad_corr_factors =
      rebin_corrfactors_to_rradtable_grid(param->rrad_corr_factors, ret_fractions, status);

  emisProfile *emis_return = calc_rrad_emis_corona(ret_fractions, rrad_corr_factors,
                                                   emis_input_rebinned, param->gamma, status);
  CHECK_STATUS_RET(*status, NULL);

  emisProfile *emis_return_rebinned = new_emisProfile(emis_input->re, emis_input->nr, status);
  rebin_emisprofile_on_radial_grid(emis_return_rebinned, emis_return, status);

  free_emisProfile(emis_return);
  free_emisProfile(emis_input_rebinned);
  free_returningFractions(&ret_fractions);
  free_rrad_corr_factors(&rrad_corr_factors);

  return emis_return_rebinned;
}

