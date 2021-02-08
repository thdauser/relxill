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
#include "relmodels.h"
#include "relutility.h"
#include "relbase.h"
#include "writeOutfiles.h"


static double calc_rrad_emis_corona_singleZone(returnFracIpol *dat,
                                               emisProfile *emisInput,
                                               double gamma,
                                               int irad,
                                               int *status) {

  // loop over all zones to get the emissivity incident on this one zone irad

  int ng = dat->tabData->ng;
  int nrad = dat->nrad;

  // need to have emisInput on the radial zone grid, ASCENDING (as tables are ascending in radius)
  // [poor test, as radial grid could still be different]
  assert(emisInput->nr == dat->nrad);
  assert((dat->rlo[irad] < emisInput->re[irad]) && (emisInput->re[irad] < dat->rhi[irad]));
  assert(emisInput->re[0] < emisInput->re[1]);

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
 * uses an emissivity profile (with radius descending) and interpolates it (logarithmically)
 * onto a grid of descending radial values (as required by the relconv/relxill kernel)
 * warning: values will be extrapolated onto the finer grid
 */
void interpolEmisProfile(emisProfile *emisReb, emisProfile *emis0, int *status) {

  double *re = emisReb->re;
  double *r0 = emis0->re;

  assert(re[0] > re[1]);                // output requires a descending grid :-(
  assert(r0[0] > r0[1]);  // verify that the input-grid is DESCENDING as well


  int ind_rin = inv_binary_search(r0, emis0->nr, re[emisReb->nr - 1]);
  assert(ind_rin >= 0);

  if (ind_rin <= 1) {
    RELXILL_ERROR("interpolating returning radiation emissivity profile failed", status);
    printf("  -> smallest radial value of new grid [%e] is above the second largest of the original radial grid [%e]\n",
           re[emis0->nr - 1], r0[1]);
    printf("     Rin and Rout have to be further apart \n");
    return;
  }

  int kk = ind_rin;
  for (int ii = emisReb->nr - 1; ii >= 0; ii--) {
    while ((re[ii] >= r0[kk]) && kk >= 0) {
      kk--;
      if (kk == 0) {
        if (re[ii] - RELTABLE_MAX_R <= 1e-6) {
          kk = 0;
          break;
        } else {
          RELXILL_ERROR("interpolation of rel_table on fine radial grid failed due to corrupted grid", status);
          printf("   --> radius %.4e ABOVE the maximal possible radius of %.4e \n",
                 re[ii], RELTABLE_MAX_R);
          CHECK_STATUS_VOID(*status);
        }
      } else if (kk < 0) {  // if we are already below the lowest value, we will extrapolate and therefore reset to kk=0
        kk = 0;
        break;
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

/* main function to calculate the returning radation emissivitiy profile
 *
 * Caveat: Emissivity Profiles are Define with radius DESCENDING
 */
emisProfile *get_rrad_emis_corona(double *re, int nr, const relParam *param, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  returnFracIpol *dat = get_rr_fractions(param->a, param->rin, param->rout, status);

  double rmeanZone[dat->nrad];
  for (int ii = 0; ii < dat->nrad; ii++) {
    rmeanZone[dat->nrad - ii - 1] = 0.5 * (dat->rlo[ii] + dat->rhi[ii]);
  }
  emisProfile *emisInput = new_emisProfile(rmeanZone, dat->nrad, status); // think about free memeory etc
  get_emis_jet(emisInput, param, status);
  write_data_to_file("test_emis_profile_input.dat", emisInput->re, emisInput->emis, emisInput->nr);
  CHECK_STATUS_RET(*status, NULL);

  emisProfile *emisReturn = new_emisProfile(rmeanZone, dat->nrad, status); // think about free memeory etc

  invertArray(emisReturn->re, emisReturn->nr);
  invertArray(emisReturn->emis, emisReturn->nr);

  for (int ii = 0; ii < dat->nrad; ii++) {
    emisReturn->emis[ii] = calc_rrad_emis_corona_singleZone(dat, emisInput, param->gamma, ii, status);
  }

  invertArray(emisReturn->re, emisReturn->nr);
  invertArray(emisReturn->emis, emisReturn->nr);

  emisProfile *emisReturnRebinned = new_emisProfile(re, nr, status);


  // rebin to (typically finer) emissivity grid
  interpolEmisProfile(emisReturnRebinned, emisReturn, status);

  free_emisProfile(emisReturn);
  free_emisProfile(emisInput);

  CHECK_STATUS_RET(*status, emisReturnRebinned);

  return emisReturnRebinned;
}

void relxill_returnKernel(double *ener_inp,
                          double *spec_inp,
                          int n_ener_inp,
                          xillParam *xill_param,
                          relParam *rel_param,
                          int *status) {

  CHECK_STATUS_VOID(*status);
  assert(xill_param->model_type == MOD_TYPE_RELXILLLPRET);

}

void init_par_relxill_ret(relParam **rel_param,
                          xillParam **xill_param,
                          const double *inp_par,
                          const int n_parameter,
                          int *status) {

  relParam *param = new_relParam(MOD_TYPE_RELXILLLPRET, EMIS_TYPE_LP, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLLPRET, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILLLPRET);

  param->height = inp_par[0];
  param->a = inp_par[1];
  param->incl = inp_par[2] * M_PI / 180;
  param->rin = inp_par[3];
  param->rout = inp_par[4];
  param->z = inp_par[5];
  xparam->z = inp_par[5];

  param->gamma = inp_par[6];
  xparam->gam = inp_par[6];
  xparam->lxi = inp_par[7];
  xparam->afe = inp_par[8];
  xparam->ect = inp_par[9];
  xparam->dens = 15.0;

  xparam->refl_frac = inp_par[10];
  xparam->fixReflFrac = (int) (inp_par[11] + 0.5); // make sure there is nor problem with integer conversion

  param->beta = 0.0;
  param->return_rad = (int) (inp_par[12] + 0.5); // make sure there is nor problem with integer conversion

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

}

/** RELXILL MODEL FUNCTION for the BB returning radiation **/
void tdrelxilllpret(const double *ener0,
                    int n_ener0,
                    double *photar,
                    const double *parameter,
                    int n_parameter,
                    int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxill_ret(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  double *ener = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);

  relxill_returnKernel(ener, photar, n_ener0, xill_param, rel_param, status);
  CHECK_STATUS_VOID(*status);

  free(ener);
  free(xill_param);
  free(rel_param);

}

