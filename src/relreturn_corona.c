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

    Copyright 2020 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "relreturn_corona.h"

#include "relreturn_datastruct.h"
#include "relutility.h"
#include "relbase.h"
#include "writeOutfiles.h"

/**
 *
 * @param ret_fractions
 * @param emis_input
 * @param gamma
 * @param irad
 * @return
 */

static double calc_rrad_emis_corona_singleZone(const returningFractions *ret_fractions, const emisProfile* emis_input, double gamma, int irad) {


  int ng = ret_fractions->tabData->ng;
  int nrad = ret_fractions->nrad;

  // need to have emis_input on the SAME radial zone grid, ASCENDING (as tables are ascending in radius)
  //  -> require the grid of the emissivity profile to be identical, means they share the same pointer
  assert(emis_input->re == ret_fractions->rad);
  assert(emis_input->nr == ret_fractions->nrad);
  assert(emis_input->re[0] < emis_input->re[1]);

  double emis_r[nrad];
  double gfac[ng];

  for (int ii = 0; ii < nrad; ii++) {
    get_gfac_grid(gfac, ret_fractions->tabData->gmin[irad][ii], ret_fractions->tabData->gmax[irad][ii], ng);

    emis_r[ii] = 0.0;
    for (int jj = 0; jj < ng; jj++) {
      emis_r[ii] += pow(gfac[jj], gamma) * ret_fractions->tabData->frac_g[irad][ii][jj];
    }

    // TODO:: f_ret needs to be interpolated for spin
    emis_r[ii] *= ret_fractions->frac_i[irad][ii] * emis_input->emis[ii] * ret_fractions->tabData->f_ret[ii];
  }

  return calcSum(emis_r, nrad);
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


/**
 * @brief: main function to calculate the returning radiation emissivity profile
 *
 * @detail Caveat: The radial grid for the relline profile, which is typically used for the emissivity
 * profile, is defined with the radius DESCENDING
 */
emisProfile *get_rrad_emis_corona(const emisProfile* emisInput, const relParam* param, int *status) {

  CHECK_STATUS_RET(*status, NULL);


  double rlo_emis, rhi_emis;
  determine_rlo_rhi(emisInput, &rlo_emis, &rhi_emis);

  returningFractions *fractions = get_rr_fractions(param->a, rlo_emis, rhi_emis, status);

  emisProfile* emis_input_rebinned = new_emisProfile(fractions->rad, fractions->nrad, status); // fractions->rad is not owned by emisReturn
  inv_rebin_mean(emisInput->re,emisInput->emis, emisInput->nr,
                 emis_input_rebinned->re, emis_input_rebinned->emis, emis_input_rebinned->nr, status);

  emisProfile* emisReturn = new_emisProfile(fractions->rad, fractions->nrad, status); // fractions->rad is not owned by emisReturn
  for (int irad = 0; irad < fractions->nrad; irad++) {
    emisReturn->emis[irad] = calc_rrad_emis_corona_singleZone(fractions, emis_input_rebinned, param->gamma, irad);
  }


  emisProfile *emisReturnRebinned = new_emisProfile(emisInput->re, emisInput->nr, status);
  rebin_emisprofile_on_radial_grid(emisReturnRebinned, emisReturn, status);

  free_emisProfile(emisReturn);

  return emisReturnRebinned;
}

