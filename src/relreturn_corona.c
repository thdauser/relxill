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


static double calc_rrad_emis_corona_singleZone(returningFractions *dat, const emisProfile* emis_input, double gamma, int irad) {

  // loop over all zones to get the emissivity incident on this one zone irad

  int ng = dat->tabData->ng;
  int nrad = dat->nrad;

  // need to have emisInput on the radial zone grid, ASCENDING (as tables are ascending in radius)
  // [poor test, as radial grid could still be different]
  assert(emis_input->nr == dat->nrad);
//  assert((dat->rlo[irad] < emisInput->re[irad]) && (emisInput->re[irad] < dat->rhi[irad]));
//  assert(emisInput->re[0] < emisInput->re[1]);

  double emis_r[nrad];
  double gfac[ng];

  for (int ii = 0; ii < nrad; ii++) {
    get_gfac_grid(gfac, dat->tabData->gmin[irad][ii], dat->tabData->gmax[irad][ii], ng);

    emis_r[ii] = 0.0;
    for (int jj = 0; jj < ng; jj++) {
      emis_r[ii] += pow(gfac[jj], gamma) * dat->tabData->frac_g[irad][ii][jj];
    }

    // TODO:: f_ret needs to be interpolated for spin
    emis_r[ii] *= dat->frac_i[irad][ii] * emis_input->emis[ii] * dat->tabData->f_ret[ii];
    // emis_r[ii] *= dat->frac_i[irad][ii]  * dat->tabData->f_ret[ii];
  }

  double emisZone = calcSum(emis_r, nrad);

  char buf[100];
  sprintf(buf, "test_rrad_single_zone_%02i.dat",irad);
  write_data_to_file(buf, dat->rad, emis_r, nrad);

  return emisZone;
}


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

  write_data_to_file("test_emis_profile_input.dat", emisInput->re, emisInput->emis, emisInput->nr);

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

