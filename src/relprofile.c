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
#include "relprofile.h"

#include "relutility.h"
#include "writeOutfiles.h"

cnode *cache_syspar = NULL;

/** global parameters, which can be used for several calls of the model */
relTable *ptr_rellineTable = NULL;
RelSysPar *cached_tab_sysPar = NULL;

// precision to calculate gstar from [H:1-H] instead of [0:1]
const double GFAC_H = 5e-3;

// interpolate the table in the A-MU0 plane (for one value of radius)
static void interpol_a_mu0(int ii, double ifac_a, double ifac_mu0, int ind_a,
                           int ind_mu0, RelSysPar *sysPar, relTable *reltab) {
  sysPar->gmin[ii] = interp_lin_2d_float(ifac_a, ifac_mu0,
                                         reltab->arr[ind_a][ind_mu0]->gmin[ii],
                                         reltab->arr[ind_a + 1][ind_mu0]->gmin[ii],
                                         reltab->arr[ind_a][ind_mu0 + 1]->gmin[ii],
                                         reltab->arr[ind_a + 1][ind_mu0 + 1]->gmin[ii]);

  sysPar->gmax[ii] = interp_lin_2d_float(ifac_a, ifac_mu0,
                                         reltab->arr[ind_a][ind_mu0]->gmax[ii],
                                         reltab->arr[ind_a + 1][ind_mu0]->gmax[ii],
                                         reltab->arr[ind_a][ind_mu0 + 1]->gmax[ii],
                                         reltab->arr[ind_a + 1][ind_mu0 + 1]->gmax[ii]);

  int jj;
  for (jj = 0; jj < reltab->n_g; jj++) {
    sysPar->trff[ii][jj][0] = interp_lin_2d_float(ifac_a, ifac_mu0,
                                                  reltab->arr[ind_a][ind_mu0]->trff1[ii][jj],
                                                  reltab->arr[ind_a + 1][ind_mu0]->trff1[ii][jj],
                                                  reltab->arr[ind_a][ind_mu0 + 1]->trff1[ii][jj],
                                                  reltab->arr[ind_a + 1][ind_mu0 + 1]->trff1[ii][jj]);

    sysPar->trff[ii][jj][1] = interp_lin_2d_float(ifac_a, ifac_mu0,
                                                  reltab->arr[ind_a][ind_mu0]->trff2[ii][jj],
                                                  reltab->arr[ind_a + 1][ind_mu0]->trff2[ii][jj],
                                                  reltab->arr[ind_a][ind_mu0 + 1]->trff2[ii][jj],
                                                  reltab->arr[ind_a + 1][ind_mu0 + 1]->trff2[ii][jj]);

    sysPar->cosne[ii][jj][0] = interp_lin_2d_float(ifac_a, ifac_mu0,
                                                   reltab->arr[ind_a][ind_mu0]->cosne1[ii][jj],
                                                   reltab->arr[ind_a + 1][ind_mu0]->cosne1[ii][jj],
                                                   reltab->arr[ind_a][ind_mu0 + 1]->cosne1[ii][jj],
                                                   reltab->arr[ind_a + 1][ind_mu0 + 1]->cosne1[ii][jj]);

    sysPar->cosne[ii][jj][1] = interp_lin_2d_float(ifac_a, ifac_mu0,
                                                   reltab->arr[ind_a][ind_mu0]->cosne2[ii][jj],
                                                   reltab->arr[ind_a + 1][ind_mu0]->cosne2[ii][jj],
                                                   reltab->arr[ind_a][ind_mu0 + 1]->cosne2[ii][jj],
                                                   reltab->arr[ind_a + 1][ind_mu0 + 1]->cosne2[ii][jj]);

  }
}


/* function interpolating the rel table values for rin,rout,mu0,incl   */
static RelSysPar *interpol_relTable(double a, double mu0, double rin, double rout,
                                    int *status) {

  // load tables
  if (ptr_rellineTable == NULL) {
    print_version_number();
    CHECK_STATUS_RET(*status, NULL);
    read_relline_table(RELTABLE_FILENAME, &ptr_rellineTable, status);
    CHECK_STATUS_RET(*status, NULL);
  }
  relTable *tab = ptr_rellineTable;
  assert(tab != NULL);

  double rms = kerr_rms(a);

  // make sure the desired rmin is within bounds and order correctly
  assert(rout > rin);
  assert(rin >= rms);


  /**************************************/
  /** 1 **  Interpolate in A-MU0 plane **/
  /**************************************/

  // get a structure to store the values from the interpolation in the A-MU0-plane
  if (cached_tab_sysPar == NULL) {
    cached_tab_sysPar = new_relSysPar(tab->n_r, tab->n_g, status);
    CHECK_STATUS_RET(*status, NULL);
  }

  int ind_a = binary_search_float(tab->a, tab->n_a, (float) a);
  int ind_mu0 = binary_search_float(tab->mu0, tab->n_mu0, (float) mu0);

  float ifac_a = ((float) a - tab->a[ind_a]) /
      (tab->a[ind_a + 1] - tab->a[ind_a]);
  float ifac_mu0 = ((float) mu0 - tab->mu0[ind_mu0]) /
      (tab->mu0[ind_mu0 + 1] - tab->mu0[ind_mu0]);

  // we perform tests that the input values are consistent when loading the parameters
  // this is just to double check (otherwise bad things happen, like an interpol. mu0>1!)
  assert(ifac_mu0 >= 0);
  assert(ifac_mu0 <= 1);
  assert(ifac_a >= 0);
  assert(ifac_a <= 1);

  /** get the radial grid (the radial grid only changes with A by the table definition) **/
  assert(fabsf(tab->arr[ind_a][ind_mu0]->r[tab->n_r - 1]
                   - tab->arr[ind_a][ind_mu0]->r[tab->n_r - 1]) < 1e-6);
  int ii;
  for (ii = 0; ii < tab->n_r; ii++) {
    cached_tab_sysPar->re[ii] = interp_lin_1d(ifac_a,
                                              tab->arr[ind_a][ind_mu0]->r[ii], tab->arr[ind_a + 1][ind_mu0]->r[ii]);
  }
  // we have problems for the intermost radius due to linear interpolation (-> set to RISCO)
  if ((cached_tab_sysPar->re[tab->n_r - 1] > kerr_rms(a)) &&
      ((cached_tab_sysPar->re[tab->n_r - 1] - kerr_rms(a)) / cached_tab_sysPar->re[tab->n_r - 1] < 1e-3)) {
    //		printf(" re-setting RIN from %.3f to %.3f (dr = %.2e)\n",cached_tab_sysPar->re[tab->n_r-1],kerr_rms(a),
    //		(cached_tab_sysPar->re[tab->n_r-1]-kerr_rms(a))/cached_tab_sysPar->re[tab->n_r-1]);
    cached_tab_sysPar->re[tab->n_r - 1] = kerr_rms(a);
  }

  // get the extent of the disk (indices are defined such that tab->r[ind+1] <= r < tab->r[ind]
  int ind_rmin = inv_binary_search(cached_tab_sysPar->re, tab->n_r, rin);
  int ind_rmax = inv_binary_search(cached_tab_sysPar->re, tab->n_r, rout);

  int jj;
  int kk;
  for (ii = 0; ii < tab->n_r; ii++) {
    // TODO: SHOULD WE ONLY INTERPOLATE ONLY THE VALUES WE NEED??? //
    // only interpolate values where we need them (radius is defined invers!)
    if (ii <= ind_rmin || ii >= ind_rmax + 1) {
      interpol_a_mu0(ii, ifac_a, ifac_mu0, ind_a, ind_mu0, cached_tab_sysPar, tab);
    } else {  // set everything we won't need to 0 (just to be sure)
      cached_tab_sysPar->gmin[ii] = 0.0;
      cached_tab_sysPar->gmax[ii] = 0.0;
      for (jj = 0; jj < tab->n_g; jj++) {
        for (kk = 0; kk < 2; kk++) {
          cached_tab_sysPar->trff[ii][jj][kk] = 0.0;
          cached_tab_sysPar->cosne[ii][jj][kk] = 0.0;
        }
      }
    }
  }

  /****************************/
  /** 2 **  Bin to Fine Grid **/
  /****************************/

  //  need to initialize and allocate memory
  RelSysPar *sysPar = new_relSysPar(N_FRAD, tab->n_g, status);
  CHECK_STATUS_RET(*status, NULL);
  get_fine_radial_grid(rin, rout, sysPar->re, sysPar->nr);

  /** we do not have rmax=1000.0 in the table, but just values close to it so let's do this trick**/
  double rmax = 1000.0;
  if (cached_tab_sysPar->re[ind_rmax] < rmax && cached_tab_sysPar->re[ind_rmax] * 1.01 > rmax) {
    cached_tab_sysPar->re[ind_rmax] = rmax;
  }

  // let's try to be as efficient as possible here (note that "r" DEcreases)
  assert(ind_rmin > 0); // as defined inverse, re[ind_rmin+1] is the lowest value
  assert((cached_tab_sysPar->re[ind_rmin + 1] <= rin));
  assert((cached_tab_sysPar->re[ind_rmin] >= rin));
  assert((cached_tab_sysPar->re[ind_rmax + 1] <= rout));
  assert((cached_tab_sysPar->re[ind_rmax] >= rout));
  assert(ind_rmax <= ind_rmin);
  assert(rout <= 1000.0);

  double ifac_r;
  int ind_tabr = ind_rmin;

  for (ii = sysPar->nr - 1; ii >= 0; ii--) {
    while ((sysPar->re[ii] >= cached_tab_sysPar->re[ind_tabr])) {
      ind_tabr--;
      if (ind_tabr < 0) { //TODO: construct table such that we don't need this?
        if (sysPar->re[ii] - RELTABLE_MAX_R <= 1e-6) {
          ind_tabr = 0;
          break;
        } else {
          RELXILL_ERROR("interpolation of rel_table on fine radial grid failed due to corrupted grid", status);
          printf("   --> radius %.4e ABOVE the maximal possible radius of %.4e \n",
                 sysPar->re[ii], RELTABLE_MAX_R);
          CHECK_STATUS_RET(*status, NULL);
        }
      }
    }

    ifac_r = (sysPar->re[ii] - cached_tab_sysPar->re[ind_tabr + 1])
        / (cached_tab_sysPar->re[ind_tabr] - cached_tab_sysPar->re[ind_tabr + 1]);
    // assert(ifac_r>=0.0);

    // we only allow extrapolation (i.e. ifac_r < 0) for the last bin
    if (ifac_r > 1.0 && ind_tabr > 0) {
      RELXILL_ERROR("interpolation of rel_table on fine radial grid failed due to corrupted grid", status);
      printf("   --> radius %.4e not found in [%.4e,%.4e]  \n",
             sysPar->re[ii], cached_tab_sysPar->re[ind_tabr + 1], cached_tab_sysPar->re[ind_tabr]);
      CHECK_STATUS_RET(*status, NULL);
    }

    for (jj = 0; jj < sysPar->ng; jj++) {
      for (kk = 0; kk < 2; kk++) {
        sysPar->trff[ii][jj][kk] =
            interp_lin_1d(ifac_r,
                          cached_tab_sysPar->trff[ind_tabr + 1][jj][kk],
                          cached_tab_sysPar->trff[ind_tabr][jj][kk]);

        sysPar->cosne[ii][jj][kk] =
            interp_lin_1d(ifac_r,
                          cached_tab_sysPar->cosne[ind_tabr + 1][jj][kk],
                          cached_tab_sysPar->cosne[ind_tabr][jj][kk]);

      }
    }
    sysPar->gmin[ii] =
        interp_lin_1d(ifac_r, cached_tab_sysPar->gmin[ind_tabr + 1], cached_tab_sysPar->gmin[ind_tabr]);

    sysPar->gmax[ii] =
        interp_lin_1d(ifac_r, cached_tab_sysPar->gmax[ind_tabr + 1], cached_tab_sysPar->gmax[ind_tabr]);

  }

  return sysPar;
}


/**  calculate all relativistic system parameters, including interpolation
 *   of the rel-table, and the emissivity; caching is implemented
 *   Input: relParam* param   Output: relSysPar* system_parameter_struct
 */

/* function to get the system parameters */
static RelSysPar *calculate_system_parameters(relParam *param, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  // only re-do the interpolation if rmin,rmax,a,mu0 changed
  // or if the cached parameters are NULL

  double mu0 = cos(param->incl);
  RelSysPar *sysPar = interpol_relTable(param->a, mu0, param->rin, param->rout, status);
  CHECK_STATUS_RET(*status, NULL);

  if (param->limb != 0) {
    sysPar->limb_law = param->limb;
  }


  // get emissivity profile
  sysPar->emis = calc_emis_profile(sysPar->re, sysPar->nr, param, status);

#ifdef RRAD
  if (param->return_rad == 1) {
    sysPar->emisReturn = get_rrad_emis_corona(sysPar->re, sysPar->nr, param, status);
    for (int ii=0; ii < sysPar->nr; ii++){
      sysPar->emis->emis[ii] += sysPar->emisReturn->emis[ii];
    }
  }
#endif

  if (*status != EXIT_SUCCESS) {
    RELXILL_ERROR("failed to calculate the system parameters", status);
  }

  return sysPar;
}

RelSysPar *get_system_parameters(relParam *param, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  inpar *sysinp = set_input_syspar(param, status);
  CHECK_STATUS_RET(*status, NULL);

  cache_info *ca_info = cli_check_cache(cache_syspar, sysinp, check_cache_syspar, status);
  CHECK_STATUS_RET(*status, NULL);

  RelSysPar *sysPar = NULL;
  if (ca_info->syscache == 1) {
    // system parameter values are cached, so we can take it from there
    sysPar = ca_info->store->data->relSysPar;
    if (is_debug_run()) {
      printf(" DEBUG:  SYSPAR-Cache: re-using calculated values\n");
    }
  } else {
    // NOT CACHED, so we need to calculate the system parameters
    sysPar = calculate_system_parameters(param, status);
    CHECK_STATUS_RET(*status, NULL);

    // now add (i.e., prepend) the current calculation to the cache
    set_cache_syspar(&cache_syspar, param, sysPar, status);

    if (is_debug_run() && *status == EXIT_SUCCESS) {
      printf(" DEBUG:  Adding new SYSPAR values to cache; the count is  %i \n", cli_count_elements(cache_syspar));
    }
  }

  free(ca_info);
  free(sysinp);

  CHECK_RELXILL_DEFAULT_ERROR(status);

  // make a sanity check for now
  if (*status == EXIT_SUCCESS) {
    assert(cache_syspar != NULL);
    assert(sysPar != NULL);
  }

  return sysPar;
}

/** get new structure to store the relline spectrum (possibly for several zones)
    important note: ener has n_ener+1 number of bins **/
rel_spec *new_rel_spec(int nzones, const int n_ener, int *status) {

  rel_spec *spec = (rel_spec *) malloc(sizeof(rel_spec));
  CHECK_MALLOC_RET_STATUS(spec, status, NULL)

  spec->n_zones = nzones;
  spec->n_ener = n_ener;

  spec->flux = (double **) malloc(spec->n_zones * sizeof(double *));
  CHECK_MALLOC_RET_STATUS(spec->flux, status, spec)

  int ii;

  for (ii = 0; ii < spec->n_zones; ii++) {
    spec->flux[ii] = (double *) malloc(n_ener * sizeof(double));
    CHECK_MALLOC_RET_STATUS(spec->flux[ii], status, spec)
  }

  spec->ener = (double *) malloc((spec->n_ener + 1) * sizeof(double));
  CHECK_MALLOC_RET_STATUS(spec->ener, status, spec)

  spec->rgrid = NULL; // will be allocated later
  spec->rel_cosne = NULL; // will be allocated later (only if need)

  return spec;
}

/** get new structure to store the cosne distribution spectrum (possibly for several zones) **/
RelCosne *new_rel_cosne(int nzones, int n_incl, int *status) {

  RelCosne *spec = (RelCosne *) malloc(sizeof(RelCosne));
  CHECK_MALLOC_RET_STATUS(spec, status, NULL)

  spec->n_zones = nzones;
  spec->n_cosne = n_incl;

  spec->cosne = (double *) malloc(spec->n_cosne * sizeof(double));
  CHECK_MALLOC_RET_STATUS(spec->cosne, status, spec)

  spec->dist = (double **) malloc(spec->n_zones * sizeof(double *));
  CHECK_MALLOC_RET_STATUS(spec->dist, status, spec)

  int ii;

  for (ii = 0; ii < spec->n_zones; ii++) {
    spec->dist[ii] = (double *) malloc(spec->n_cosne * sizeof(double));
    CHECK_MALLOC_RET_STATUS(spec->dist[ii], status, spec)
  }

  return spec;
}

/** initialize the rel_spec structure **/
void init_rel_spec(rel_spec **spec, relParam *param, xillTable *xill_tab, double *radialZones,
                   double **pt_ener, const int n_ener, int *status) {

  CHECK_STATUS_VOID(*status);

  /** in case of the relxill-LP model multiple zones are used **/
  int nzones = param->num_zones;

  if ((*spec) == NULL) {
    (*spec) = new_rel_spec(nzones, n_ener, status);
  } else {
    // check if the number of zones changed or number of energy bins
    if ((nzones != (*spec)->n_zones) || (n_ener != (*spec)->n_ener)) {
      // -> if yes, we need to free memory and re-allocate the appropriate space
      free_rel_spec((*spec));
      (*spec) = new_rel_spec(nzones, n_ener, status);
    }
  }

  for (int ii = 0; ii <= (*spec)->n_ener; ii++) {
    (*spec)->ener[ii] = (*pt_ener)[ii];
  }

  if (xill_tab != NULL) {
    if ((*spec)->rel_cosne == NULL) {
      (*spec)->rel_cosne = new_rel_cosne(nzones, xill_tab->n_incl, status);
    }
    for (int ii = 0; ii < (*spec)->rel_cosne->n_cosne; ii++) {
      (*spec)->rel_cosne->cosne[ii] = cos(xill_tab->incl[ii] * M_PI / 180);
    }
  }

  // if the grid changed, we called new_rel_spec
  if ((*spec)->rgrid == NULL) {
    (*spec)->rgrid = radialZones;
  } else {
    free(radialZones);
  }

  CHECK_RELXILL_DEFAULT_ERROR(status);

}

static void zero_rel_spec_flux(rel_spec *spec) {
  int ii;
  int jj;
  for (ii = 0; ii < spec->n_zones; ii++) {
    for (jj = 0; jj < spec->n_ener; jj++) {
      spec->flux[ii][jj] = 0.0;
    }
    if (spec->rel_cosne != NULL) {
      for (jj = 0; jj < spec->rel_cosne->n_cosne; jj++) {
        spec->rel_cosne->dist[ii][jj] = 0.0;
      }
    }
  }
}

/** relat. transfer function, which we will need to integrate over the energy bin then **/
static str_relb_func *new_str_relb_func(RelSysPar *sysPar, int *status) {
  str_relb_func *str = (str_relb_func *) malloc(sizeof(str_relb_func));
  CHECK_MALLOC_RET_STATUS(str, status, NULL)

  str->gstar = sysPar->gstar;
  str->ng = sysPar->ng;
  str->limb_law = 0;

  return str;
}

/** relat. function which we want to integrate **/
static double relb_func(double eg, int k, str_relb_func *str) {

  // get the redshift from the energy
  // double eg = e/line_energy;
  double egstar = (eg - str->gmin) * str->del_g;

  // find the indices in the original g-grid, but check first if they have already been calculated
  int ind;
  if (!((egstar >= str->gstar[str->save_g_ind]) && (egstar < str->gstar[str->save_g_ind + 1]))) {
    str->save_g_ind = binary_search(str->gstar, str->ng, egstar);
  }
  ind = str->save_g_ind;

  double inte = (egstar - str->gstar[ind]) / (str->gstar[ind + 1] - str->gstar[ind]);
  double inte1 = 1.0 - inte;
  double ftrf = inte * str->trff[ind][k] + inte1 * str->trff[ind + 1][k];

  double val = pow(eg, 3) / ((str->gmax - str->gmin) * sqrt(egstar - egstar * egstar)) * ftrf * str->emis;

  /** isotropic limb law by default (see Svoboda (2009)) **/
  if (str->limb_law == 0) {
    return val;
  } else {
    double fmu0 = inte * str->cosne[ind][k] + inte1 * str->cosne[ind + 1][k];
    double limb = 1.0;
    if (str->limb_law == 1) { //   !Laor(1991)
      limb = (1.0 + 2.06 * fmu0);
    } else if (str->limb_law == 2) {  //  !Haardt (1993)
      limb = log(1.0 + 1.0 / fmu0);
    }
    return val * limb;
  }
}

/** Romberg Integration Routine **/
static double romberg_integration(double a, double b, int k, str_relb_func *str) {
  const double prec = 0.02;
  double obtprec = 1.0;
  const int itermin = 0;
  int itermax = 5;
  enum {
    maxiter = 6
  };
  double t[maxiter + 1][maxiter + 1];

  int niter = 0;

  if (itermax > maxiter) {
    itermax = maxiter;
  }

  // check if this value has already been calculated
  double r;
  if (str->cached_relbf) {
    r = str->cache_val_relb_func[k];
  } else {
    r = relb_func(a, k, str);
  }

  str->cache_val_relb_func[k] = relb_func(b, k, str);
  str->cache_rad_relb_fun = str->re;
  // rb(k) = RELB_FUNC(b,k);


  int ii;

  double ta = (r + str->cache_val_relb_func[k]) / 2.0;
  double pas = b - a;
  double pasm = 1.0;
  double s;
  t[0][0] = ta * pas;
  while ((niter < itermin) || ((obtprec > prec) && (niter <= itermax))) {
    niter++;
    pas = pas / 2.0;
    pasm = pasm / 2.0;
    s = ta;
    for (ii = 1; ii <= pow(2, niter) - 1; ii++) {
      s += relb_func(a + pas * ii, k, str);
    }
    t[0][niter] = s * pas;
    r = 1.0;
    for (ii = 1; ii <= niter; ii++) {
      r *= 4.0;
      int jj = niter - ii;
      t[ii][jj] = (r * t[ii - 1][jj + 1] - t[ii - 1][jj]) / (r - 1.0);
    }
    obtprec = fabs(t[niter][0] - t[niter - 1][0]) / t[niter][0];
  }

  return t[niter][0];
}

/** function which makes an approximated integration for gstar->0/1
    this is only done within gstar=[0,H] and gstar[H,1-H]
    input:   bin_lo and bin_hi
    output:  area of the bin (= luminosity = E/dt/bin) )  **/
static double int_edge(double blo, double bhi, double h, str_relb_func *str, double line_energy) {


  // get the value of the Luminosity on the point closest to the ones to be approximated ("H")

  double hex;
  double lo;
  double hi;

  // #1: upper or lower limit -> write the corresponding gstar-value to hex
  if (blo <= 0.5) {
    hex = h;
    lo = blo;
    hi = bhi;
  } else {
    hex = 1.0 - h;
    /**  variable transformation for upper limit x -> 1 - x
         leads to the same integral
         (if the correct normalization is chosen; hex keeps track of that) **/
    lo = 1.0 - bhi;
    hi = 1.0 - blo;
  }

  // #2: get the normalization value
  int k = 0;
  double norm = 0.0;
  for (k = 0; k < 2; k++) {
    norm = norm + relb_func(gstar2ener(hex, str->gmin, str->gmax, line_energy), k, str);
  }


  // #3: thanks to variable transformation:
  //   both cases are described by the one for [0,H]
  norm = norm * sqrt(h);

  return 2 * norm * (sqrt(hi) - sqrt(lo)) * line_energy * (str->gmax - str->gmin);
}

/** function which calculates the normal integral by romberg's method
(it is acutally jsut a wrapper which sets the correct parameters
and coordinates the integration over k=1,2)
input:   bin_lo and bin_hi
output:  area of the bin (= luminosity = E/dt/bin) ) **/
static double int_romb(double lo, double hi, str_relb_func *str, double line_energy) {

  double flu = 0.0;

  /** We are doing a trick here: for the "red wing" of the line, one can show that a simple trapez integration is enough,
   * as the profile is very smooth. In order to avoid problems for narrow, double peaked lines, the red wing is defined
   * to start at 0.95*line_energy **/
  int k;
  if (lo >= line_energy * 0.95) {
    for (k = 0; k < 2; k++) {
      flu += romberg_integration(lo, hi, k, str);
    }
  } else {
    for (k = 0; k < 2; k++) {
      flu += relb_func((hi + lo) / 2.0, k, str) * (hi - lo);
    }
  }

  return flu;
}

/** integrate the flux bin (see Dauser+2010, MNRAS for details) **/
static double integ_relline_bin(str_relb_func *str, double rlo0, double rhi0) {

  double line_ener = 1.0;
  double flu = 0.0;

  double gblo = (rlo0 / line_ener - str->gmin) * str->del_g;
  if (gblo < 0.0) {
    gblo = 0.0;
  } else if (gblo > 1.0) {
    gblo = 1.0;
  }

  double gbhi = (rhi0 / line_ener - str->gmin) * str->del_g;
  if (gbhi < 0.0) {
    gbhi = 0.0;
  } else if (gbhi > 1.0) {
    gbhi = 1.0;
  }
  if (gbhi == 0) {
    return 0.0;
  }

  double rlo = rlo0;
  double rhi = rhi0;

  double hlo;
  double hhi;
  // #1: low approx. integration
  if (gblo <= GFAC_H) {
    // range of the low-app.-integration
    hlo = gblo;
    hhi = GFAC_H;
    // begin of the 'real' integration
    rlo = gstar2ener(GFAC_H, str->gmin, str->gmax, line_ener);
    // .. but also check if this integration is necessary
    if (gbhi <= GFAC_H) {
      hhi = gbhi;
      rlo = -1.0;
    }
    // approximated integration
    flu = flu + int_edge(hlo, hhi, GFAC_H, str, line_ener);
  }

  // #2: upper limit approximation to be taken into account?
  if (gbhi >= (1.0 - GFAC_H)) {
    // range of the upper-app.-integration
    hhi = gbhi;
    hlo = 1.0 - GFAC_H;

    // begin of the 'real' integration
    rhi = gstar2ener(1 - GFAC_H, str->gmin, str->gmax, line_ener);
    // .. but also check if this integration is necessary
    if (gblo >= (1.0 - GFAC_H)) {
      hlo = gblo;
      rhi = -1.0;
    }

    /** the same approximated integration as in case #1 can be
           applied, if one makes a variable transformation **/
    flu = flu + int_edge(hlo, hhi, GFAC_H, str, line_ener);
  }

  // #3: real integration (only if necessary)
  if ((rhi >= 0) && (rlo >= 0)) {

    // has the function relb_func been calculated at the lower bin boundary before?
    // (should've been upper bound before; also make sure we haven't changed the radial bin!)
    if ((fabs(rlo - str->cache_bin_ener) < CACHE_LIMIT) && (fabs(str->re - str->cache_rad_relb_fun) < CACHE_LIMIT)) {
      str->cached_relbf = 1;
    } else {
      str->cached_relbf = 0;
    }
    flu = flu + int_romb(rlo, rhi, str, line_ener);
  }

  return flu;
}

static void set_str_relbf(str_relb_func *str, double re, double gmin, double gmax, double **trff,
                          double **cosne, double emis, int limb_law) {
  str->re = re;
  str->gmin = gmin;
  str->gmax = gmax;
  str->del_g = 1. / (gmax - gmin);
  str->emis = emis;

  str->trff = trff;
  str->cosne = cosne;

  str->cache_bin_ener = -1.0;
  str->cache_rad_relb_fun = -1.0;
  str->cached_relbf = 0;

  str->limb_law = limb_law;

  str->save_g_ind = 0;
}

/** function to properly re-normalize the relline_profile **/
void renorm_relline_profile(rel_spec *spec, relParam *rel_param, const int *status) {

  CHECK_STATUS_VOID(*status);

  // normalize to 'cts/bin'
  int ii;
  int jj;
  double sum = 0.0;
  for (ii = 0; ii < spec->n_zones; ii++) {
    for (jj = 0; jj < spec->n_ener; jj++) {
      spec->flux[ii][jj] /= 0.5 * (spec->ener[jj] + spec->ener[jj + 1]);
      sum += spec->flux[ii][jj];
    }
  }

  /** only renormalize if not the relxill model or not a lamp post model **/

  if (do_renorm_model(rel_param)) {

    double relline_norm = 1;
    if (is_relxill_model(rel_param->model_type) && rel_param->emis_type == EMIS_TYPE_BKN) {
      relline_norm = norm_factor_semi_infinite_slab(rel_param->incl * 180.0 / M_PI);
    }

    if (is_debug_run()) {
      printf(" DEBUG: re-normalizing output spectrum\n");
    }
    for (ii = 0; ii < spec->n_zones; ii++) {
      for (jj = 0; jj < spec->n_ener; jj++) {
        spec->flux[ii][jj] *= relline_norm / sum;
      }
    }
  }

  if (spec->rel_cosne != NULL) {
    for (ii = 0; ii < spec->n_zones; ii++) {
      // normalize it for each zone, the overall flux will be taken care of by the normal structure
      sum = 0.0;
      for (jj = 0; jj < spec->rel_cosne->n_cosne; jj++) {
        sum += spec->rel_cosne->dist[ii][jj];
      }
      for (jj = 0; jj < spec->rel_cosne->n_cosne; jj++) {
        spec->rel_cosne->dist[ii][jj] /= sum;
      }
    }
  }
}

/** get the bin for a certain emission angle (between [0,n_incl-1] **/
int static get_cosne_bin(double mu, RelCosne *dat) {
  return ((int) (dat->n_cosne * (1 - mu) + 1)) - 1;
}

/** calculate the relline profile(s) for all given zones **/
str_relb_func *cached_str_relb_func = NULL;

static double calculate_radiallyResolvedFluxObs(str_relb_func *relb_func, rel_spec *spec, double weight) {

  double integRadialFlux = 0.0;
  for (int jj = 0; jj <= spec->n_ener; jj++) {
    integRadialFlux += integ_relline_bin(cached_str_relb_func, spec->ener[jj], spec->ener[jj + 1]);
  }

  integRadialFlux *= weight;

  return integRadialFlux;

}

static void free_str_relb_func(str_relb_func **str) {
  if (*str != NULL) {
    free(*str);
    *str = NULL;
  }
}

void relline_profile(rel_spec *spec, RelSysPar *sysPar, int *status) {

  CHECK_STATUS_VOID(*status);

  double line_ener = 1.0;

  // very important: set all fluxes to zero
  zero_rel_spec_flux(spec);

  if (cached_str_relb_func == NULL) {
    cached_str_relb_func = new_str_relb_func(sysPar, status);
  }

  // store the (energy)-integrated flux in an array for debugging
  double *radialFlux = NULL;
  if (shouldOutfilesBeWritten()) {
    radialFlux = (double *) malloc(sizeof(double) * sysPar->nr);
    CHECK_MALLOC_VOID_STATUS(radialFlux, status)
  }

  int ii;
  int jj;
  for (ii = 0; ii < sysPar->nr; ii++) {
    // gstar in [0,1] + corresponding energies (see gstar2ener for full formula)
    double egmin = sysPar->gmin[ii] * line_ener;
    double egmax = sysPar->gmax[ii] * line_ener;

    // check if the expected energy-bins are needed
    if ((egmax > spec->ener[0]) && (egmin < spec->ener[spec->n_ener])) {


      /**  make sure that integration is only done inside the
           given energy range **/
      if (egmin < spec->ener[0]) {
        egmin = spec->ener[0];
      }
      if (egmax > spec->ener[spec->n_ener]) {
        egmax = spec->ener[spec->n_ener];
      }

      /** search for the indices in the ener-array
          index is such that: ener[k]<=e<ener[k+1] **/
      int ielo = binary_search(spec->ener, spec->n_ener + 1, egmin);
      int iehi = binary_search(spec->ener, spec->n_ener + 1, egmax);

      // in which ionization bin are we?
      int izone = binary_search(spec->rgrid, spec->n_zones + 1, sysPar->re[ii]);

      // set the current parameters in a cached structure (and reset some values) [optimizes speed]
      set_str_relbf(cached_str_relb_func,
                    sysPar->re[ii], sysPar->gmin[ii], sysPar->gmax[ii],
                    sysPar->trff[ii], sysPar->cosne[ii],
                    sysPar->emis->emis[ii], sysPar->limb_law);

      /** INTEGRATION
       *   [remember: defintion of Xillver/Relxill is 1/2 * Speith Code]
       *   [remember: trapez integration returns just r*dr*PI (full integral is over dA/2)]
       *   [ -> in the end it's weigth=PI*r*dr/2 ]
       */
      double weight = trapez_integ_single(sysPar->re, ii, sysPar->nr) / 2;

      // lastly, loop over the energies
      for (jj = ielo; jj <= iehi; jj++) {
        double tmp_var = integ_relline_bin(cached_str_relb_func, spec->ener[jj], spec->ener[jj + 1]);
        spec->flux[izone][jj] += tmp_var * weight;
      }

      if (shouldOutfilesBeWritten() && spec->n_zones == 1) {
        assert(radialFlux != NULL);
        radialFlux[ii] = calculate_radiallyResolvedFluxObs(cached_str_relb_func, spec, weight);
      }

      /** only calculate the distribution if we need it here  **/
      if (spec->rel_cosne != NULL) {
        int kk;
        int imu;
        str_relb_func *da = cached_str_relb_func; // define a shortcut
        for (jj = 0; jj < sysPar->ng; jj++) {
          double g = da->gstar[jj] * (da->gmax - da->gmin) + da->gmin;
          for (kk = 0; kk < 2; kk++) {
            imu = get_cosne_bin(da->cosne[jj][kk], spec->rel_cosne);

            spec->rel_cosne->dist[izone][imu] +=
                da->re * pow(2 * M_PI * g * da->re, 2) /
                    sqrt(da->gstar[jj] - da->gstar[jj] * da->gstar[jj]) *
                    da->trff[jj][kk] * da->emis
                    * weight * sysPar->d_gstar[jj];
          }
        }
      } /** end calculating angular distribution **/


    }
  }

  /** we need to free the structure as it points to the currently cached sysPar structure
       which is freed if the cache is full and therefore causes "invalid reads" **/
  free_str_relb_func(&cached_str_relb_func);

  if (shouldOutfilesBeWritten() && spec->n_zones == 1) {
    save_relline_radial_flux_profile(sysPar->re, radialFlux, sysPar->nr);
  }

  CHECK_RELXILL_DEFAULT_ERROR(status);

}

RelSysPar *new_relSysPar(int nr, int ng, int *status) {
  RelSysPar *sysPar = (RelSysPar *) malloc(sizeof(RelSysPar));
  CHECK_MALLOC_RET_STATUS(sysPar, status, NULL)

  sysPar->ng = ng;
  sysPar->nr = nr;

  sysPar->re = (double *) malloc(nr * sizeof(double));
  CHECK_MALLOC_RET_STATUS(sysPar->re, status, sysPar)
  sysPar->gmin = (double *) malloc(nr * sizeof(double));
  CHECK_MALLOC_RET_STATUS(sysPar->gmin, status, sysPar)
  sysPar->gmax = (double *) malloc(nr * sizeof(double));
  CHECK_MALLOC_RET_STATUS(sysPar->gmax, status, sysPar)

  sysPar->emis = NULL;
  sysPar->emisReturn = NULL;

  sysPar->gstar = (double *) malloc(ng * sizeof(double));
  CHECK_MALLOC_RET_STATUS(sysPar->gstar, status, sysPar)

  // we already set the values as they are fixed
  int ii;
  int jj;
  for (ii = 0; ii < ng; ii++) {
    sysPar->gstar[ii] = GFAC_H + (1.0 - 2 * GFAC_H) / (ng - 1) * ((float) (ii));
  }

  sysPar->d_gstar = (double *) malloc(ng * sizeof(double));
  CHECK_MALLOC_RET_STATUS(sysPar->gstar, status, sysPar)
  for (ii = 0; ii < ng; ii++) {
    if ((ii == 0) || (ii == (ng - 1))) {
      sysPar->d_gstar[ii] = 0.5 * (sysPar->gstar[1] - sysPar->gstar[0]) + GFAC_H;
    } else {
      sysPar->d_gstar[ii] = sysPar->gstar[1] - sysPar->gstar[0];
    }
  }

  sysPar->trff = (double ***) malloc(nr * sizeof(double **));
  CHECK_MALLOC_RET_STATUS(sysPar->trff, status, sysPar)
  sysPar->cosne = (double ***) malloc(nr * sizeof(double **));
  CHECK_MALLOC_RET_STATUS(sysPar->cosne, status, sysPar)
  for (ii = 0; ii < nr; ii++) {
    sysPar->trff[ii] = (double **) malloc(ng * sizeof(double *));
    CHECK_MALLOC_RET_STATUS(sysPar->trff[ii], status, sysPar)
    sysPar->cosne[ii] = (double **) malloc(ng * sizeof(double *));
    CHECK_MALLOC_RET_STATUS(sysPar->cosne[ii], status, sysPar)
    for (jj = 0; jj < ng; jj++) {
      sysPar->trff[ii][jj] = (double *) malloc(2 * sizeof(double));
      CHECK_MALLOC_RET_STATUS(sysPar->trff[ii][jj], status, sysPar)
      sysPar->cosne[ii][jj] = (double *) malloc(2 * sizeof(double));
      CHECK_MALLOC_RET_STATUS(sysPar->cosne[ii][jj], status, sysPar)
    }
  }

  sysPar->limb_law = 0;

  return sysPar;
}

void free_relSysPar(RelSysPar *sysPar) {
  if (sysPar != NULL) {
    free(sysPar->re);
    free(sysPar->gmin);
    free(sysPar->gmax);
    free(sysPar->gstar);
    free(sysPar->d_gstar);

    free_emisProfile(sysPar->emis);

    if (sysPar->trff != NULL) {
      int ii;
      for (ii = 0; ii < sysPar->nr; ii++) {
        if (sysPar->trff[ii] != NULL) {
          int jj;
          for (jj = 0; jj < sysPar->ng; jj++) {
            free(sysPar->trff[ii][jj]);
          }
          free(sysPar->trff[ii]);
        }
      }
      free(sysPar->trff);
    }

    if (sysPar->cosne != NULL) {
      int ii;
      for (ii = 0; ii < sysPar->nr; ii++) {
        if (sysPar->cosne[ii] != NULL) {
          int jj;
          for (jj = 0; jj < sysPar->ng; jj++) {
            free(sysPar->cosne[ii][jj]);
          }
          free(sysPar->cosne[ii]);
        }
      }
      free(sysPar->cosne);
    }
    free(sysPar);
  }
}

void free_cached_relTable(void) {
  free_relTable(ptr_rellineTable);
}

void free_relprofile_cache(void) {
  free_relSysPar(cached_tab_sysPar);
  free_str_relb_func(&cached_str_relb_func);
}

void free_cache_syspar(void) {
  cli_delete_list(&cache_syspar);
}
