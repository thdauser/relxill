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

#include "Rellp.h"
#include "Relphysics.h"

extern "C" {
#include "writeOutfiles.h"
}

lpTable *cached_lp_table = nullptr;


/*
 * @brief: calculate the reflection fraction as defined in Dauser+2016 *
 * @params:
 *  - emis_profile: calculated for the given Rin/Rout of the model
 *  - del_emit_ad_max: emission angle hitting the outer edge of the fully, simulated accretion disk at
 *    1000Rg, regardless of Rout (as photons are not allowed to cross the disk plane)
 **/
static lpReflFrac *calc_refl_frac(emisProfile *emis_profile,
                                  double rin,
                                  double rout,
                                  double del_emit_ad_max,
                                  double beta,
                                  int *status) {

  // get the angle emitted in the rest-frame of the primary source, which hits the inner and outer edge of the disk
  //  -> important: get the radial values for which the RELLINE is calculated
  double del_bh = emis_profile->del_emit[inv_binary_search(emis_profile->re, emis_profile->nr, rin)];
  double del_ad = emis_profile->del_emit[inv_binary_search(emis_profile->re, emis_profile->nr, rout)];

  // photons are not allowed to cross the disk plane, so no photon (for beta=0) emitted
  // at del<pi/2 is allowed to be counted to reaching infinity
  // (this can happen for a large height, as Rout_simulation=1000Rg)
  if (del_emit_ad_max < M_PI / 2.0) {
    del_emit_ad_max = M_PI / 2.0;
  }

  /** calculate the coordinate transformation / relat abberation
   *   - an observer on the accretion disk sees the rays from
   *     del_bh up to del_ad
   *   - for the reflection fraction we therefore need to convert from
   *     the moving source (which the disk observer sees) into the
   *     local frame
   *   -> therefore we need to calculate the abberation of -beta
   */
  if (beta > 1e-6) {
    del_bh = relat_abberation(del_bh, -1. * beta);
    del_ad = relat_abberation(del_ad, -1. * beta);
  }

  lpReflFrac *str = new_lpReflFrac(status);
  CHECK_STATUS_RET(*status, str);

  str->f_bh = 0.5 * (1.0 - cos(del_bh));
  str->f_ad = 0.5 * (cos(del_bh) - cos(del_ad));

  // photons are not allowed to cross the disk plane, so we need the angle del_emit_ad_max
  str->f_inf_rest = 0.5 * (1.0 + cos(del_emit_ad_max));

  // in case of a moving source, we also have f_inf different from f_inf_rest (need for calculating the refl_frac)
  if (beta > 1e-6) {
    str->f_inf = 0.5 * (1.0 + cos(relat_abberation(del_emit_ad_max, -1. * beta)));
  } else {
    str->f_inf = str->f_inf_rest;
  }

  assert(str->f_inf + str->f_ad < 1.0);

  str->refl_frac = str->f_ad / str->f_inf;

  return str;
}

static void norm_emis_profile(const double *re, const int nr, double *emis) {

  double integ_area = 0.0;
  auto delta_area = new double[nr];

  for (int ii = 0; ii < nr; ii++) {
    if (re[1] < re[0]) {  // desencing grid
      delta_area[ii] = trapez_integ_single(re, ii, nr) * 2;
    } else {
      delta_area[ii] = trapez_integ_single_rad_ascending(re, ii, nr) * 2;
    }
    integ_area += emis[ii] * delta_area[ii];
  }
  delete[] delta_area;

  for (int ii = 0; ii < nr; ii++) {
    emis[ii] /= integ_area;
  }
}


/**
 *
 * @param emis_prof (required to be descending in radius)
 * @param emis_prof_tab (required to be ascending in radius)
 * @param status
 *
 * @detail function "invert_emis_profile" can be used to convert
 */
void rebin_emisprofile_on_radial_grid(emisProfile *emis_prof, const emisProfile* emis_prof_tab, int *status) {

  if (is_emis_grid_ascending(emis_prof)==1){
    RELXILL_ERROR("rebinning emissivity profile failed (require output radial grid of emissivity to be descending with radius",
                  status);
    assert(emis_prof->re[0]>emis_prof->re[1]);
    return;
  }

  if (is_emis_grid_ascending(emis_prof_tab)==0){
    RELXILL_ERROR("rebinning emissivity profile failed (require input emissivity to be ASCENDING with radius",
                  status);
    assert(emis_prof_tab->re[1]>emis_prof_tab->re[0]);
    return;
  }


  double* re = emis_prof->re;
  int nr = emis_prof->nr;

  double* re_tab = emis_prof_tab->re;
  int nr_tab = emis_prof_tab->nr;

  assert(emis_prof->re[0]>emis_prof->re[1]); // decreasing radius in input emis profile
  assert(emis_prof_tab->re[0]<emis_prof_tab->re[1]); // increasing radius in the table

  // get the extent of the disk (indices are defined such that tab->r[ind] <= r < tab->r[ind+1]
  int ind_rmin = binary_search(re_tab, nr_tab, re[ nr - 1]);

  assert(ind_rmin >= 0);
  assert(ind_rmin < nr_tab - 1);
  int kk = ind_rmin;
  for (int ii = nr - 1; ii >= 0; ii--) {
    while ((re[ii] >= re_tab[kk + 1])) {
      kk++;
      if (kk >= nr_tab - 1) { //TODO: construct table such that we don't need this?
        if (re[ii] - RELTABLE_MAX_R <= 1e-6) {
          kk = nr_tab - 2;
          break;
        } else {
          RELXILL_ERROR("interpolation of rel_table on fine radial grid failed due to corrupted grid", status);
          printf("   --> radius %.4e ABOVE the maximal possible radius of %.4e \n",
                 re[ii], RELTABLE_MAX_R);
          CHECK_STATUS_VOID(*status);
        }
      }
    }

    const double inter_r = get_ipol_factor_radius(re_tab[kk], re_tab[kk + 1], emis_prof_tab->del_emit[kk], re[ii]);

    //  log grid for the intensity (due to the function profile)
    emis_prof->emis[ii] = interp_log_1d(inter_r, emis_prof_tab->emis[kk], emis_prof_tab->emis[kk + 1]);
    emis_prof->del_emit[ii] = interp_lin_1d(inter_r, emis_prof_tab->del_emit[kk], emis_prof_tab->del_emit[kk + 1]);
    emis_prof->del_inc[ii] = interp_lin_1d(inter_r, emis_prof_tab->del_inc[kk], emis_prof_tab->del_inc[kk + 1]);
  }
}



int is_emis_grid_ascending(const emisProfile* emis){
  if (emis->re[0]<emis->re[emis->nr-1]){
    return 1;
  } else {
    return 0;
  }
}


// interpolate everything for the given a-h values on the original grid
// as h(a), first need to interpolate in h and then for a
static void interpol_emisprofile_spin_height(emisProfile *emis_profile_table,
                                             double ifac_a,
                                             const double *ifac_h,
                                             const int *ind_h,
                                             lpDat *const *dat_a,
                                             const int nr) {
  for (int ii = 0; ii < nr; ii++) {
    emis_profile_table->emis[ii] =
        (1.0 - ifac_a) * interp_lin_1d(ifac_h[0],dat_a[0]->intens[ind_h[0]][ii],dat_a[0]->intens[ind_h[0] + 1][ii]  )
            + (ifac_a) * interp_lin_1d(ifac_h[1],dat_a[1]->intens[ind_h[1]][ii],dat_a[1]->intens[ind_h[1] + 1][ii]  );

    emis_profile_table->del_emit[ii] =
        (1.0 - ifac_a) * interp_lin_1d(ifac_h[0],dat_a[0]->del[ind_h[0]][ii],dat_a[0]->del[ind_h[0] + 1][ii]  )
            + (ifac_a) * interp_lin_1d(ifac_h[1],dat_a[1]->del[ind_h[1]][ii],dat_a[1]->del[ind_h[1] + 1][ii]  );

    emis_profile_table->del_inc[ii] =
        (1.0 - ifac_a) * interp_lin_1d(ifac_h[0],dat_a[0]->del_inc[ind_h[0]][ii],dat_a[0]->del_inc[ind_h[0] + 1][ii]  )
            + (ifac_a) * interp_lin_1d(ifac_h[1],dat_a[1]->del_inc[ind_h[1]][ii],dat_a[1]->del_inc[ind_h[1] + 1][ii]  );

  }
}

/*
 * @function:  interpol_lptable
 * @synoposis: interpolate the LP table for a given spin and height
 * @output: emisProfile
 *   - the "emis" here is the photon flux, as stored in the table
 *   - radial grid is the one of the table
 */
emisProfile* interpol_lptable(double a, double height, lpTable* tab, int* status){

  CHECK_STATUS_RET(*status, nullptr);

  int ind_a;
  double ifac_a;
  get_ipol_factor((float) a, tab->a, tab->n_a, &ind_a, &ifac_a);


  lpDat *dat_ind_a[2] = {tab->dat[ind_a], tab->dat[ind_a+1]};

  auto* jet_rad = (double*) malloc(sizeof(double)*tab->n_rad);
  CHECK_MALLOC_RET_STATUS(jet_rad,status,nullptr)
  for (int ii = 0; ii < tab->n_rad; ii++) {
    jet_rad[ii] = interp_lin_1d(ifac_a, dat_ind_a[0]->rad[ii], dat_ind_a[1]->rad[ii]);
  }

  emisProfile* emis_profile_table = new_emisProfile(jet_rad, tab->n_rad, status);

  int ind_h[2];
  double ifac_h[2];
  for (int ii = 0; ii < 2; ii++) {
    get_ipol_factor((float) height, dat_ind_a[ii]->h, tab->n_h, &(ind_h[ii]), &(ifac_h[ii]));
  }

  interpol_emisprofile_spin_height(emis_profile_table, ifac_a, ifac_h, ind_h, dat_ind_a, tab->n_rad);

  return emis_profile_table;
}

void apply_emis_fluxboost_source_disk(emisProfile *emisProf, double a, double height, double gamma, double beta) {
  for (int ii = 0; ii < emisProf->nr; ii++) {
    emisProf->emis[ii] *= calc_fluxboost_source_disk(emisProf->re[ii], emisProf->del_emit[ii], a, height, gamma, beta);
  }
}

/**
 * @brief calculate the emissivity profile of a lamp post point source
 * Important: from the relParam input values, height and beta will be ignored (as this allows
 * the routine to be called for different values if height and beta for an extended source)
 * @param emis_profile
 * @param param
 * @param height
 * @param beta
 * @param tab
 * @param status
 */
static void calc_emis_jet_point_source(emisProfile *emis_profile, const relParam *param, double height, double beta,
                                       lpTable *tab, int *status) {
  CHECK_STATUS_VOID(*status);

  emisProfile *emis_profile_table = interpol_lptable(param->a, height, tab, status);

  rebin_emisprofile_on_radial_grid(emis_profile, emis_profile_table, status);

  // calculate the angle under which photons are emitted from the source such that they hit the outer edge of
  // the simulated accretion disk (i.e., photons with a larger emission angle are able to reach the observer)
  const double del_emit_ad_max = emis_profile_table->del_emit[tab->n_rad - 1];
  emis_profile->photon_fate_fractions =
      calc_refl_frac(emis_profile, param->rin, param->rout, del_emit_ad_max, param->beta, status);
  free(emis_profile_table->re); // is not freed by free_emisProfile
  free_emisProfile(emis_profile_table);

  apply_emis_fluxboost_source_disk(emis_profile, param->a, height, param->gamma, beta);

  emis_profile->normFactorPrimSpec = 0.0; // currently not used, calculated directly in add_primary_component
}

int modelLampPostPointsource(const relParam *param) {
  const double htopPrecLimit = 1e-3;
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
  CHECK_MALLOC_RET_STATUS(source, status, source)

  get_log_grid(source->heightArr, source->nh + 1, param->height, param->htop);

  // check and set the parameters as defined for the extended jet
  assert(param->height < param->htop);
  double beta100Rg = param->beta;
  for (int ii = 0; ii < source->nh; ii++) {
    source->heightMean[ii] = 0.5 * (source->heightArr[ii] + source->heightArr[ii + 1]);
    if (beta100Rg > 1e-6) {
      source->beta[ii] = jetSpeedConstantAccel(beta100Rg, source->heightMean[ii], param->height);
    } else {
      source->beta[ii] = 0.0;
    }
  }

  return source;
}

static void addSingleReturnFractions(lpReflFrac *reflFracAvg, lpReflFrac *singleReflFrac, double fraction) {

  reflFracAvg->refl_frac += singleReflFrac->refl_frac * fraction;
  reflFracAvg->f_ad += singleReflFrac->f_ad * fraction;
  reflFracAvg->f_inf_rest += singleReflFrac->f_inf_rest * fraction;
  reflFracAvg->f_bh += singleReflFrac->f_bh * fraction;

}

/*
 *  EXTENDED LAMP POST:  [deprecated, currently not used in any model]
 *  - if htop <= heigh=hbase we assume it's a point-like jet
 *  - the meaning of beta for the extended jet is the velocity at 100Rg, in case the
 *    profile is of interest, it will be output in the debug mode
 */
void calc_emis_jet_extended(emisProfile *emisProf,
                            const relParam *param,
                            lpTable *tab,
                            int *status) {

  extPrimSource *source = getExtendedJetGeom(param, status);
  CHECK_STATUS_VOID(*status);

  emisProfile *emisProfSingle = new_emisProfile(emisProf->re, emisProf->nr, status);

  setArrayToZero(emisProf->emis, emisProf->nr);
  emisProf->photon_fate_fractions = new_lpReflFrac(status);
  emisProf->normFactorPrimSpec = 0.0;


  for (int ii = 0; ii < source->nh; ii++) {

    calc_emis_jet_point_source(emisProfSingle,
                               param,
                               source->heightMean[ii],
                               source->beta[ii],
                               tab,
                               status);

    // assuming a constant luminosity in the frame of the jet
    double
        heightIntegrationFactor = (source->heightArr[ii + 1] - source->heightArr[ii]) / (param->htop - param->height);

    for (int jj = 0; jj < emisProf->nr; jj++) {
      emisProf->emis[jj] += emisProfSingle->emis[jj] * heightIntegrationFactor;
      emisProf->del_inc[jj] += emisProfSingle->del_inc[jj] * heightIntegrationFactor;
      emisProf->del_emit[jj] += emisProfSingle->del_emit[jj] * heightIntegrationFactor;
    }

    addSingleReturnFractions(emisProf->photon_fate_fractions,
                             emisProfSingle->photon_fate_fractions,
                             heightIntegrationFactor);

    emisProf->normFactorPrimSpec += emisProfSingle->normFactorPrimSpec * heightIntegrationFactor;

    free_lpReflFrac(&(emisProfSingle->photon_fate_fractions));
  }

  free_extendedPrimarySource(source);
  free_emisProfile(emisProfSingle);

}


static lpTable* get_lp_table(int* status){
  CHECK_STATUS_RET(*status,nullptr);

  if (cached_lp_table == nullptr) {
    read_lp_table(LPTABLE_FILENAME, &cached_lp_table, status);
    CHECK_STATUS_RET(*status, nullptr);
  }
  return cached_lp_table;
}



/**
 *  @synopsis: calculate the emissivity for broken power law defined
 *  by up to two indices
 **/
void get_emis_bkn(double *emis, const double *re, int nr,
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


/**
 *  @synopsis: calculate the emissivity for an alpha-disk (as given
 *  by Shakura & Sunyaev; 1973)
 **/
static void get_emis_alphadisk(double *emis, double *re, int nr) {

  // requires descending grid in radius
  assert(re[nr-1] < re[0]);
  assert(re[1] < re[0]);

  for (int ii = 0; ii < nr; ii++) {
    emis[ii] = 1. / pow(re[ii], 3) * (1 - 1. / sqrt(re[ii] / re[nr-1]));
  }

  norm_emis_profile(re, nr, emis);
}

static void get_emis_constant(double *emis, int n) {
  for (int ii = 0; ii < n; ii++) {
    emis[ii] = 1.0;
  }
}


/**
  * @syopsis: calculate the emissivity profile for any jet like source
  * (extended, lamp post, ...)
  **/
void get_emis_jet(emisProfile *emis_profile, const relParam *param, int *status) {

  CHECK_STATUS_VOID(*status);

  lpTable* lp_table = get_lp_table(status);

  if (modelLampPostPointsource(param)) {
    calc_emis_jet_point_source(emis_profile, param, param->height, param->beta, lp_table, status);
  } else {
    calc_emis_jet_extended(emis_profile, param, lp_table, status);
  }

}


static void add_returnrad_emis(const relParam* param, emisProfile* emis0, int* status) {

  if ( is_debug_run() ) {
    write_data_to_file("test-rrad-emis-input.dat", emis0->re, emis0->emis, emis0->nr);
  }

  emisProfile* emisReturn = get_rrad_emis_corona(emis0, param, status);

  if ( is_debug_run() ) {
    write_data_to_file("test-rrad-emis-rrad.dat",emis0->re,emisReturn->emis, emis0->nr);
  }

  assert(emisReturn->nr == emis0->nr);

  for (int ii=0; ii < emis0->nr; ii++){
    if (param->return_rad == 1) {
      emis0->emis[ii] += emisReturn->emis[ii];
    } else if (param->return_rad == -1 || param->return_rad == 2) { // only return rad (for debugging only)
      emis0->emis[ii] = emisReturn->emis[ii];
    } else {
      RELXILL_ERROR("adding returning radiation failed ", status);
      printf("    return_rad = %i is not allowed \n", param->return_rad);
    }
  }

  free_emisProfile(emisReturn);
}



/*
 *  @function: calc_emis_profile
 *  @synopsis: calculate the emissivity profile on a given radial grid,
 *  depending on the given emis-type and parameters
 */
emisProfile *calc_emis_profile(double *re, int nr, const relParam *param, int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  emisProfile *emis = new_emisProfile(re, nr, status);

  const double invalid_angle = -1.0;

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
    return nullptr;
  }


  if (abs(param->return_rad) > 1e-6) {
    add_returnrad_emis(param, emis, status);
  }


  CHECK_RELXILL_ERROR("calculating the emissivity profile failed due to previous error", status);

  return emis;
}



// constructors and destructors //

extPrimSource *new_extendedPrimarySource(int nh, int *status) {

  auto *source = new extPrimSource;
  CHECK_MALLOC_RET_STATUS(source, status, source)

  source->nh = nh;
  source->heightArr = new double[nh + 1];
  source->heightMean = new double[nh];
  source->beta = new double[nh];

  return source;
}

void free_extendedPrimarySource(extPrimSource *source) {
  if (source != nullptr) {
    free(source->heightArr);
    free(source->heightMean);
    free(source->beta);
    free(source);
  }

}

lpReflFrac *new_lpReflFrac(int *status) {

  auto *str = new lpReflFrac;
  CHECK_MALLOC_RET_STATUS(str, status, nullptr)
  str->refl_frac = 0.0;
  str->f_inf_rest = 0.0;
  str->f_ad = 0.0;
  str->f_bh = 0.0;

  return str;
}

void free_lpReflFrac(lpReflFrac **str) {
  if (*str != nullptr) {
    delete (*str);
    *str = nullptr;
  }
}

emisProfile *new_emisProfile(double *re, int nr, int *status) {

  auto *emis = (emisProfile *) malloc(sizeof(emisProfile));
  CHECK_MALLOC_RET_STATUS(emis, status, nullptr)

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
    emis->del_emit[ii] = 0.0;
    emis->del_inc[ii] = 0.0;
  }

  emis->normFactorPrimSpec = 0.0;

  emis->photon_fate_fractions = nullptr;

  return emis;
}

void free_emisProfile(emisProfile *emis_profile) {
  if (emis_profile != nullptr) {
    free(emis_profile->emis);
    free(emis_profile->del_emit);
    free(emis_profile->del_inc);

    free_lpReflFrac(&(emis_profile->photon_fate_fractions));

    free(emis_profile);
  }
}

void free_cached_lpTable() {
  free_lpTable(cached_lp_table);
}
