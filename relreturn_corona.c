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
#include "relphysics.h"
#include "relutility.h"
#include "relbase.h"
#include "common.h"
#include "relreturn.h"

static double calc_rrad_emis_corona_singleZone(returnFracIpol *dat,
                                               emisProfile *emisInput,
                                               double gamma,
                                               int irad,
                                               int *status) {

  // loop over all zones to get the emissivity incident on this one zone irad

  int ng = dat->tabData->ng;
  int nrad = dat->nrad;

  // need to have emisInput on the radial zone grid
  // [poor test, as radial grid could still be different]
  assert(emisInput->nr == dat->nrad);

  double emis_r[nrad];
  double gfac[ng];

  for (int ii = 0; ii < nrad; ii++) {

    get_gfac_grid(gfac, dat->tabData->gmin[irad][ii], dat->tabData->gmax[irad][ii], ng);

    emis_r[ii] = 0.0;
    for (int jj = 0; jj < ng; jj++) {
      emis_r[ii] += pow(gfac[jj], gamma) * dat->tabData->frac_g[irad][ii][jj];
    }

    emis_r[ii] *= dat->frac_i[irad][ii] * emisInput->emis[ii] * dat->tabData->f_ret[ii];
  }

  double emisZone = calcSum(emis_r, nrad);

  return emisZone;
}
/* function: interpolEmisProfile
 * uses an emissivity profile (with radius ascending) and interpolates it (logarithmically)
 * onto a grid of descending radial values (as required by the relconv/relxill kernel)
 */
void interpolEmisProfile(emisProfile *emisReb, emisProfile *emis0, int *status) {

  double *re = emisReb->re;

  assert(re[0] > re[1]); // output requires a descending grid :-(

  // verify that the input-grid is ascending
  assert(emis0->re[1] > emis0->re[0]);
  double *r0 = emis0->re;

  int ind_kmin = binary_search(r0, emis0->nr, re[emisReb->nr - 1]);
  assert(ind_kmin >= 0);

  if (ind_kmin >= emis0->nr - 2) {
    RELXILL_ERROR("interpolating emissivity profile failed", status);
    printf("  -> smallest radial value of new grid [%e] is above the second largest of the original radial grid [%e]\n",
           re[0], r0[emis0->nr - 2]);
    printf("     Rin and Rout have to be further apart \n");
    return;
  }

  int kk = ind_kmin;
  for (int ii = emisReb->nr - 1; ii >= 0; ii--) {
    while ((re[ii] >= r0[kk + 1])) {
      kk++;
      if (kk == emis0->nr - 1) {
        if (re[ii] - RELTABLE_MAX_R <= 1e-6) {
          kk = emis0->nr - 2;
          break;
        } else {
          RELXILL_ERROR("interpolation of rel_table on fine radial grid failed due to corrupted grid", status);
          printf("   --> radius %.4e ABOVE the maximal possible radius of %.4e \n",
                 re[ii], RELTABLE_MAX_R);
          CHECK_STATUS_VOID(*status);
        }
      }
    }

    // for larger angles logarithmic interpolation works slightly better
    double inter_r;
    if (emis0->del_emit[kk] / M_PI * 180.0 <= 75.0) {
      inter_r = (re[ii] - r0[kk]) / (r0[kk + 1] - r0[kk]);
    } else {
      inter_r = (log(re[ii]) - log(r0[kk])) /
          (log(r0[kk + 1]) - log(r0[kk]));
    }

    //  log grid for the intensity (due to the function profile)
    emisReb->emis[ii] = interp_log_1d(inter_r, emis0->emis[kk], emis0->emis[kk + 1]);

    emisReb->del_emit[ii] = interp_lin_1d(inter_r, emis0->del_emit[kk], emis0->del_emit[kk + 1]);
    emisReb->del_inc[ii] = interp_lin_1d(inter_r, emis0->del_inc[kk], emis0->del_inc[kk + 1]);

  }
}
emisProfile *get_rrad_emis_corona(double *re, int nr, relParam *param, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  returnFracIpol *dat = get_rr_fractions(param->a, param->rin, param->rout, status);

  double rmeanZone[dat->nrad];
  for (int ii = 0; ii < dat->nrad; ii++) {
    rmeanZone[ii] = 0.5 * (dat->rlo[ii] + dat->rhi[ii]);
  }
  emisProfile *emisInput = new_emisProfile(rmeanZone, dat->nrad, status); // think about free memeory etc
  get_emis_jet(emisInput, param, status);
  CHECK_STATUS_RET(*status, NULL);

  emisProfile *emisReturn = new_emisProfile(rmeanZone, dat->nrad, status); // think about free memeory etc
  for (int ii = 0; ii < dat->nrad; ii++) {
    emisReturn->emis[ii] = calc_rrad_emis_corona_singleZone(dat, emisInput, param->gamma, ii, status);
  }

  emisProfile *emisReturnRebinned = new_emisProfile(re, nr, status);


  // rebin to (typically finer) emissivity grid
  interpolEmisProfile(emisReturnRebinned, emisReturn, status);

  free_emisProfile(emisReturn);
  free_emisProfile(emisInput);

  CHECK_STATUS_RET(*status, emisReturnRebinned);

  return emisReturnRebinned;
}