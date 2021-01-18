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

#include "rellp.h"

lpTable *cached_lp_table = NULL;

/** calculate the reflection fraction as defined in Dauser+2016 **/
static lpReflFrac *calc_refl_frac(emisProfile *emis_profile, double del_emit_ad_max, relParam *param, int *status) {

  // in case there is no relativity information, the refl_frac is 1
  if (param == NULL) {
    printf(" *** Warning: can not calculate reflection fraction as no relat. parameters are given \n");
    return NULL;
  }

  /** important: get the radial values for which the RELLINE is calculated
   *             should be Rin=r_isco & Rout=1000rg  **/

  // get the angle emitted in the rest-frame of the primary source, which hits the inner and outer edge of the disk
  double del_bh = emis_profile->del_emit[inv_binary_search(emis_profile->re, emis_profile->nr, param->rin)];
  double del_ad = emis_profile->del_emit[inv_binary_search(emis_profile->re, emis_profile->nr, param->rout)];

  /** calculate the coordinate transformation / relat abberation
   *   - an observer on the accretion disk sees the rays from
   *     del_bh up to del_ad
   *   - for the reflection fraction we therefore need to convert from
   *     the moving source (which the disk observer sees) into the
   *     local frame
   *   -> therefore we need to calculate the abberation of -beta
   */
  if (param->beta > 1e-6) {
    del_bh = relat_abberation(del_bh, -1. * param->beta);
    del_ad = relat_abberation(del_ad, -1. * param->beta);
  }

  lpReflFrac *str = new_lpReflFrac(status);
  CHECK_STATUS_RET(*status, str);

  str->f_bh = 0.5 * (1.0 - cos(del_bh));
  str->f_ad = 0.5 * (cos(del_bh) - cos(del_ad));
  /** photons are not allowed to cross the disk
   *  (so they only reach infinity if they don't hit the disk plane) */
  str->f_inf = 0.5 * (1.0 + cos(del_emit_ad_max));

  // photons are not allowed to cross the disk plane
  if (str->f_inf > 0.5) {
    str->f_inf = 0.5;
  }

  str->refl_frac = str->f_ad / str->f_inf;

  return str;
}

static void norm_emis_profile(const double *re, int nr, double *emis) {

  double integ_area = 0.0;
  for (int ii = 0; ii < nr; ii++) {
    integ_area += emis[ii] * trapez_integ_single(re, ii, nr) * 2;
  }

  for (int ii = 0; ii < nr; ii++) {
    emis[ii] /= integ_area;
  }

}

/** routine for the broken power law emissivity **/
static void get_emis_bkn(double *emis, const double *re, int nr,
                         double index1, double index2, double rbr) {

  double alpha;

  int ii;
  for (ii = 0; ii < nr; ii++) {
    alpha = index1;
    if (re[ii] > rbr) {
      alpha = index2;
    }
    emis[ii] = pow(re[ii] / rbr, -alpha);
  }

  norm_emis_profile(re, nr, emis);

}

double getPrimarySpecScalingFactor(double ginf, double gamma, double f_inf) {
  double prim_fac = f_inf / 0.5 * pow(ginf, gamma);
  return prim_fac;
}

/*
 * ignores param->height and param->beta
 */
static void calc_emis_jet_point_source(emisProfile *emisProf, relParam *param, double height, double beta,
                                       lpTable *tab, int ind_a, double ifac_a, int *status) {

  lpDat *dat[2];
  int ind_h[2];
  double ifac_h[2];

  for (int ii = 0; ii < 2; ii++) {
    dat[ii] = tab->dat[ind_a + ii];
    ind_h[ii] = binary_search_float(dat[ii]->h, tab->n_h, (float) height);
    ifac_h[ii] = (height - dat[ii]->h[ind_h[ii]]) /
        (dat[ii]->h[ind_h[ii] + 1] - dat[ii]->h[ind_h[ii]]);

    // make sure the incident angle is defined as positive value (otherwise the interpolation
    // will create problems / jumps )
    int jj;
    int kk;
    for (jj = 0; jj < tab->n_h; jj++) {
      for (kk = 0; kk < tab->n_rad; kk++) {
        dat[ii]->del_inc[jj][kk] = fabs(dat[ii]->del_inc[jj][kk]);
        dat[ii]->del[jj][kk] = fabs(dat[ii]->del[jj][kk]);
      }
    }

  }

  double jet_rad[tab->n_rad];
  double jet_emis[tab->n_rad];
  double jet_del[tab->n_rad];
  double jet_del_inc[tab->n_rad];

  // interpolate everything for the given a-h values on the original grid
  for (int ii = 0; ii < tab->n_rad; ii++) {
    // #1: intensity
    jet_emis[ii] =
        (1.0 - ifac_a) * (1.0 - ifac_h[0]) * dat[0]->intens[ind_h[0]][ii]
            + (1.0 - ifac_a) * (ifac_h[0]) * dat[0]->intens[ind_h[0] + 1][ii]
            + (ifac_a) * (1.0 - ifac_h[1]) * dat[1]->intens[ind_h[1]][ii]
            + (ifac_a) * (ifac_h[1]) * dat[1]->intens[ind_h[1] + 1][ii];

    jet_del[ii] =
        (1.0 - ifac_a) * (1.0 - ifac_h[0]) * dat[0]->del[ind_h[0]][ii]
            + (1.0 - ifac_a) * (ifac_h[0]) * dat[0]->del[ind_h[0] + 1][ii]
            + (ifac_a) * (1.0 - ifac_h[1]) * dat[1]->del[ind_h[1]][ii]
            + (ifac_a) * (ifac_h[1]) * dat[1]->del[ind_h[1] + 1][ii];

    jet_del_inc[ii] =
        (1.0 - ifac_a) * (1.0 - ifac_h[0]) * dat[0]->del_inc[ind_h[0]][ii]
            + (1.0 - ifac_a) * (ifac_h[0]) * dat[0]->del_inc[ind_h[0] + 1][ii]
            + (ifac_a) * (1.0 - ifac_h[1]) * dat[1]->del_inc[ind_h[1]][ii]
            + (ifac_a) * (ifac_h[1]) * dat[1]->del_inc[ind_h[1] + 1][ii];

    // #2: r-grid
    jet_rad[ii] = interp_lin_1d(ifac_a, dat[0]->rad[ii], dat[1]->rad[ii]);
  }

  // and now rebin it to the given radial grid
  double inter_r;

  double *re = emisProf->re;
  int n_r = emisProf->nr;

  // get the extent of the disk (indices are defined such that tab->r[ind+1] <= r < tab->r[ind]
  int ind_rmin = binary_search(jet_rad, tab->n_rad, re[n_r - 1]);
  assert(ind_rmin > 0);
  int kk = ind_rmin;
  for (int ii = n_r - 1; ii >= 0; ii--) {
    while ((re[ii] >= jet_rad[kk + 1])) {
      kk++;
      if (kk >= tab->n_rad - 1) { //TODO: construct table such that we don't need this?
        if (re[ii] - RELTABLE_MAX_R <= 1e-6) {
          kk = tab->n_rad - 2;
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
    if (jet_del[kk] / M_PI * 180.0 <= 75.0) {
      inter_r = (re[ii] - jet_rad[kk]) / (jet_rad[kk + 1] - jet_rad[kk]);
    } else {
      inter_r = (log(re[ii]) - log(jet_rad[kk])) /
          (log(jet_rad[kk + 1]) - log(jet_rad[kk]));
    }

    //  log grid for the intensity (due to the function profile)
    emisProf->emis[ii] = interp_log_1d(inter_r, jet_emis[kk], jet_emis[kk + 1]);

    emisProf->del_emit[ii] = interp_lin_1d(inter_r, jet_del[kk], jet_del[kk + 1]);
    emisProf->del_inc[ii] = interp_lin_1d(inter_r, jet_del_inc[kk], jet_del_inc[kk + 1]);

    /** multiply by the additional factor gi^gamma (see Dauser et al., 2013) **/
    emisProf->emis[ii] *= pow(gi_potential_lp(re[ii], param->a, height, beta, emisProf->del_emit[ii]), param->gamma);

    // take the beaming of the jet into account (see Dauser et al., 2013)
    if (param->beta > 1e-6) {
      emisProf->emis[ii] *= pow(doppler_factor(emisProf->del_emit[ii], beta), 2);
    }
  }

  // del_emit for the largest radius of the table (need for refl_frac)
  double del_emit_ad_max = jet_del[tab->n_rad - 1];
  emisProf->returnFracs = calc_refl_frac(emisProf, del_emit_ad_max, param, status);

  assert(emisProf->returnFracs != NULL);
  double g_inf = calc_g_inf(height, param->a);
  emisProf->normFactorPrimSpec = getPrimarySpecScalingFactor(g_inf, param->gamma, emisProf->returnFracs->f_inf);

}

int modelLampPostPointsource(relParam *param) {
  double htopPrecLimit = 1e-3;
  if ((fabs(param->htop) <= 1.0) || (param->htop - htopPrecLimit <= param->height)) {
    return 1;
  } else {
    return 0;
  }
}

/*
 * Calculate the velocity (in units of c) for a given height, hbase and a given
 * velocity (beta) at 100Rg. Specfici definitions are:
 *  - beta is defined as the velocity at 100r_g above the jet base at hbase
 *  - we assume constant acceleration such that currentBeta=bet100
 */
double jetSpeedConstantAccel(double beta100, double height, double hbase) {

  const double HEIGHT_REF = 100.0;  // in units of Rg

  double rel_gamma = 1.0 / sqrt(1 - pow(beta100, 2));
  double acceleration = 1.0 / HEIGHT_REF * (sqrt(1 + pow(rel_gamma * beta100, 2)) - 1);

  double x = height - hbase;
  double t = sqrt(pow(x, 2) + 2 * x / acceleration);

  // speed for constant accel in SRT
  double beta = acceleration * t / sqrt(1 + (pow(acceleration * t, 2)));

  return beta;
}

// get the extended source geometry in height and allocate necessary parameters
extPrimSource *getExtendedJetGeom(const relParam *param, int *status) {

  extPrimSource *source = new_extendedPrimarySource(NHBINS_VERTICALLY_EXTENDED_SOURCE, status);
  CHECK_MALLOC_RET_STATUS(source, status, source);

  get_log_grid(source->heightArr, source->nh + 1, param->height, param->htop);

  // check and set the parameters as defined for the extended jet
  assert(param->height < param->htop);
  double beta100Rg = param->beta;
  for (int ii = 0; ii < source->nh; ii++) {
    source->heightMean[ii] = 0.5 * (source->heightArr[ii] + source->heightArr[ii + 1]);
    if (beta100Rg > 1e-6) {
      source->beta[ii] = jetSpeedConstantAccel(beta100Rg, source->heightMean[ii], param->height);;
    } else {
      source->beta[ii] = 0.0;
    }
  }

  return source;
}

static void addSingleReturnFractions(lpReflFrac *reflFracAvg, lpReflFrac *singleReflFrac, double fraction) {

  reflFracAvg->refl_frac += singleReflFrac->refl_frac * fraction;
  reflFracAvg->f_ad += singleReflFrac->f_ad * fraction;
  reflFracAvg->f_inf += singleReflFrac->f_inf * fraction;
  reflFracAvg->f_bh += singleReflFrac->f_bh * fraction;

}

/*
 *  EXTENDED LAMP POST:
 *  - if htop <= heigh=hbase we assume it's a point-like jet
 *  - the meaning of beta for the extended jet is the velocity at 100Rg, in case the
 *    profile is of interest, it will be output in the debug mode
 */
void calc_emis_jet_extended(emisProfile *emisProf,
                            relParam *param,
                            lpTable *tab,
                            int ind_a,
                            double ifac_a,
                            int *status) {

  extPrimSource *source = getExtendedJetGeom(param, status);
  CHECK_STATUS_VOID(*status);

  emisProfile *emisProfSingle = new_emisProfile(emisProf->re, emisProf->nr, status);

  setArrayToZero(emisProf->emis, emisProf->nr);
  emisProf->returnFracs = new_lpReflFrac(status);
  emisProf->normFactorPrimSpec = 0.0;

  for (int ii = 0; ii < source->nh; ii++) {

    calc_emis_jet_point_source(emisProfSingle,
                               param,
                               source->heightMean[ii],
                               source->beta[ii],
                               tab,
                               ind_a,
                               ifac_a,
                               status);

    // assuming a constant luminosity in the frame of the jet
    double
        heightIntegrationFactor = (source->heightArr[ii + 1] - source->heightArr[ii]) / (param->htop - param->height);

    for (int jj = 0; jj < emisProf->nr; jj++) {
      emisProf->emis[jj] += emisProfSingle->emis[jj] * heightIntegrationFactor;

      emisProf->del_inc[jj] += emisProfSingle->del_inc[jj] * heightIntegrationFactor;
      emisProf->del_emit[jj] += emisProfSingle->del_emit[jj] * heightIntegrationFactor;

    }

    addSingleReturnFractions(emisProf->returnFracs, emisProfSingle->returnFracs, heightIntegrationFactor);

    emisProf->normFactorPrimSpec += emisProfSingle->normFactorPrimSpec * heightIntegrationFactor;

    free_lpReflFrac(&(emisProfSingle->returnFracs));
  }

  if (shouldOutfilesBeWritten()) {
    save_radial_profile("test_rellxill_heightVelocityProfile.txt", source->heightMean, source->beta, source->nh);
  }

  free_extendedPrimarySource(source);
  free_emisProfile(emisProfSingle);

}


static void get_emis_alphadisk(double *emis, double *re, int n) {

  for (int ii = 0; ii < n; ii++) {
    emis[ii] = 1. / pow(re[ii], 3) * (1 - 1. / sqrt(re[ii] / re[0]));
  }

  // normalized to 1?

}

static void get_emis_constant(double *emis, int n) {

  for (int ii = 0; ii < n; ii++) {
    emis[ii] = 1.0;
  }

  // normalized to 1?

}


/*
 * LAMP POST GEOMETRY  --- MAIN ROUTINE
 */
void get_emis_jet(emisProfile *emis_profile, relParam *param, int *status) {

  CHECK_STATUS_VOID(*status);

  if (cached_lp_table == NULL) {
    read_lp_table(LPTABLE_FILENAME, &cached_lp_table, status);
    CHECK_STATUS_VOID(*status);
  }

  int ind_a = binary_search_float(cached_lp_table->a, cached_lp_table->n_a, (float) param->a);
  double ifac_a = (param->a - cached_lp_table->a[ind_a]) /
      (cached_lp_table->a[ind_a + 1] - cached_lp_table->a[ind_a]);

  if (modelLampPostPointsource(param)) {
    calc_emis_jet_point_source(emis_profile, param, param->height, param->beta, cached_lp_table, ind_a, ifac_a, status);
  } else {
    calc_emis_jet_extended(emis_profile, param, cached_lp_table, ind_a, ifac_a, status);
  }

}

/*
 *  MAIN ROUTINE to calculate the EMISSIVITY PROFILE
 */
emisProfile *calc_emis_profile(double *re, int nr, relParam *param, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  emisProfile *emis = new_emisProfile(re, nr, status);

  double invalid_angle = -1.0;

  /**  *** Broken Power Law Emissivity ***  **/
  if (param->emis_type == EMIS_TYPE_BKN) {

    get_emis_bkn(emis->emis, emis->re, emis->nr,
                 param->emis1, param->emis2, param->rbr);
    // set the angles in this case to a default value
    int ii;
    for (ii = 0; ii < nr; ii++) {
      emis->del_emit[ii] = invalid_angle;
      emis->del_inc[ii] = invalid_angle;
    }

  } else if (param->emis_type == EMIS_TYPE_ALPHA) {
    get_emis_alphadisk(emis->emis, re, nr);

  } else if (param->emis_type == EMIS_TYPE_CONST) {
    get_emis_constant(emis->emis, nr);

    /**  *** Lamp Post Emissivity ***  **/
  } else if (param->emis_type == EMIS_TYPE_LP) {

    get_emis_jet(emis, param, status);
  } else {

    RELXILL_ERROR(" calculation of emissivity profile not possible \n", status);
    printf("   -> emis_type=%i not known \n", param->emis_type);
    return NULL;
  }

  CHECK_RELXILL_ERROR("calculating the emissivity profile failed due to previous error", status);

  return emis;
}

void free_cached_lpTable(void) {
  free_lpTable(cached_lp_table);
}



// constructors and destructors //

extPrimSource *new_extendedPrimarySource(int nh, int *status) {

  extPrimSource *source = malloc(sizeof(extPrimSource));
  CHECK_MALLOC_RET_STATUS(source, status, source);

  source->nh = nh;
  source->heightArr = malloc(sizeof(double) * nh + 1);
  source->heightMean = malloc(sizeof(double) * nh);
  source->beta = malloc(sizeof(double) * nh);

  return source;
}

void free_extendedPrimarySource(extPrimSource *source) {
  if (source != NULL) {
    free(source->heightArr);
    free(source->heightMean);
    free(source->beta);
    free(source);
  }

}

lpReflFrac *new_lpReflFrac(int *status) {

  lpReflFrac *str = (lpReflFrac *) malloc(sizeof(lpReflFrac));
  CHECK_MALLOC_RET_STATUS(str, status, NULL)
  str->refl_frac = 0.0;
  str->f_inf = 0.0;
  str->f_ad = 0.0;
  str->f_bh = 0.0;

  return str;
}

void free_lpReflFrac(lpReflFrac **str) {
  if (*str != NULL) {
    free(*str);
    *str = NULL;
  }
}

emisProfile *new_emisProfile(double *re, int nr, int *status) {

  emisProfile *emis = (emisProfile *) malloc(sizeof(emisProfile));
  CHECK_MALLOC_RET_STATUS(emis, status, NULL)

  emis->re = re;
  emis->nr = nr;

  emis->emis = (double *) malloc(nr * sizeof(double));
  CHECK_MALLOC_RET_STATUS(emis->emis, status, emis)
  emis->del_emit = (double *) malloc(nr * sizeof(double));
  CHECK_MALLOC_RET_STATUS(emis->del_emit, status, emis)
  emis->del_inc = (double *) malloc(nr * sizeof(double));
  CHECK_MALLOC_RET_STATUS(emis->del_inc, status, emis)

  for (int ii = 0; ii < nr; ii++) {
    emis->emis[ii] = 0.0;
  }

  emis->normFactorPrimSpec = 0.0;

  emis->returnFracs = NULL;

  return emis;
}

void free_emisProfile(emisProfile *emis_profile) {
  if (emis_profile != NULL) {
    free(emis_profile->emis);
    free(emis_profile->del_emit);
    free(emis_profile->del_inc);

    free_lpReflFrac(&(emis_profile->returnFracs));

    free(emis_profile);
  }
}

