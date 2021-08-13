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
#include "IonGradient.h"

extern "C"{
#include "relprofile.h"
#include "writeOutfiles.h"
}

// calculate the log(xi) for given density and emissivity
static double cal_lxi(double dens, double emis) {
  return log10(4.0 * M_PI * emis / dens);
}

// determine the radius of maximal ionization
static double cal_lxi_max_ss73(double *re, double *emis, int nr, double rin) {

  double rad_max_lxi = pow((11. / 9.), 2)
      * rin;  // we use the same definition as Adam with r_peak = (11/9)^2 rin to be consistent (does not matter much)

      // radial AD grid is sorted descending (!)
      int kk = inv_binary_search(re, nr, rad_max_lxi);
      double interp = (rad_max_lxi - re[kk + 1]) / (re[kk] - re[kk + 1]);

      double emis_max_lxi = interp_lin_1d(interp, emis[kk + 1], emis[kk]);

      double lxi_max = cal_lxi(density_ss73_zone_a(rad_max_lxi, rin), emis_max_lxi);

      return lxi_max;
}


/** *** set log(xi) to obey the limits of the xillver table: TODO: check if we need to adjust the normalization as well  ***
 *  NOTE: with correctly set xpsec/isis limits, it is only possible to reach the lower boundary       **/
static void lxi_set_to_xillver_bounds(double *pt_lxi) {
  /**  TODO: Need to define this globally **/
  double xlxi_tab_min = 0.0;
  double xlxi_tab_max = 4.7;

  // #1: set the value of xi to the lowest value of the table
  if (*pt_lxi < xlxi_tab_min) {
    *pt_lxi = xlxi_tab_min;
  } else if (*pt_lxi > xlxi_tab_max) {
    //	#2: high ionization: we approximately assume such a highly ionized disk acts as a mirror
    *pt_lxi = xlxi_tab_max;
  }

}



ion_grad *new_ion_grad(const double *r, int n, int *status) {

  auto ion = new ion_grad;

  ion->dens = new double[n+1]{0};

  ion->r = (double *) malloc((n + 1) * sizeof(double));
  CHECK_MALLOC_RET_STATUS(ion->r, status, nullptr)
  ion->lxi = (double *) malloc((n) * sizeof(double));
  CHECK_MALLOC_RET_STATUS(ion->lxi, status, nullptr)
  ion->fx = (double *) malloc((n) * sizeof(double));
  CHECK_MALLOC_RET_STATUS(ion->fx, status, nullptr)
  ion->del_emit = (double *) malloc((n) * sizeof(double));
  CHECK_MALLOC_RET_STATUS(ion->del_emit, status, nullptr)

  ion->nbins = n;

  int ii;
  for (ii = 0; ii < n; ii++) {
    ion->r[ii] = r[ii];
    ion->lxi[ii] = 0.0;
    ion->fx[ii] = 0.0;
    ion->del_emit[ii] = M_PI / 4.; // assume default 45 deg (xillver assumption), only used if beta>0
  }
  // radius goes to n+1
  ion->r[n] = r[n];

  return ion;
}


ion_grad *calc_ion_gradient(relParam *rel_param,
                            double xlxi0,
                            double xindex,
                            int type,
                            double *rgrid,
                            int n,
                            int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  ion_grad *ion = new_ion_grad(rgrid, n, status);
  assert(ion!=nullptr);
  CHECK_STATUS_RET(*status, nullptr);

  auto rmean = new double[n];
  auto del_inc = new double[n];
  int ii;
  for (ii = 0; ii < n; ii++) {
    rmean[ii] = 0.5 * (rgrid[ii] + rgrid[ii + 1]);
  }

  if (type == ION_GRAD_TYPE_PL) {
    for (ii = 0; ii < n; ii++) {
      ion->lxi[ii] = (exp(xlxi0))
          * pow((rmean[ii] / rmean[0]), -1.0 * xindex);  // TODO: check if we need to subtract xlxi_tab_min here
          ion->lxi[ii] = log(ion->lxi[ii]);

          lxi_set_to_xillver_bounds(&(ion->lxi[ii]));

    }

  } else if (type == ION_GRAD_TYPE_ALPHA) {

    double rin = rgrid[0];

    // TODO: use a better approach to not linearly interpolate but rather average over the profile?
    auto emis_zones = new double[n];

    // we need the emissivity profile (should be cached, so no extra effort required here)
    RelSysPar *sysPar = get_system_parameters(rel_param, status);
    emisProfile *emis_profile = sysPar->emis;

    assert(emis_profile->del_inc != nullptr);
    inv_rebin_mean(emis_profile->re, emis_profile->emis, sysPar->nr, rmean, emis_zones, n, status);
    inv_rebin_mean(emis_profile->re, emis_profile->del_inc, sysPar->nr, rmean, del_inc, n, status);
    inv_rebin_mean(emis_profile->re, emis_profile->del_emit, sysPar->nr, rmean, ion->del_emit, n, status);

    // calculate the maximal ionization assuming r^-3 and SS73 alpha disk
    double lxi_max = cal_lxi_max_ss73(emis_profile->re, emis_profile->emis, emis_profile->nr, rin);

    // the maximal ionization is given as input parameter, so we need to normalize our calculation by this value
    double fac_lxi_norm = xlxi0 - lxi_max; // subtraction instead of division because of the log

    /** calculate the density for a  stress-free inner boundary condition, i.e., R0=rin in SS73)  **/
    for (ii = 0; ii < n; ii++) {
      ion->dens[ii] = density_ss73_zone_a(rmean[ii], rin);

      // now we can use the emissivity to calculate the ionization
      ion->lxi[ii] = cal_lxi(ion->dens[ii], emis_zones[ii]) + fac_lxi_norm;

      ion->lxi[ii] += log10(cos(M_PI / 4) / cos(del_inc[ii]));

      lxi_set_to_xillver_bounds(&(ion->lxi[ii]));
    }

  } else if (type
  == ION_GRAD_TYPE_CONST) {  // should not happen, as this will be approximated by 1 zone (but just in case we get here...)
    for (ii = 0; ii < n; ii++) {
      ion->lxi[ii] = xlxi0;
    }
  } else {
    printf(" *** ionization type with number %i not implemented \n", type);
    printf("     choose either %i for the PL, %i for the ALPHA-disk, or %i for constant\n",
           ION_GRAD_TYPE_PL, ION_GRAD_TYPE_ALPHA, ION_GRAD_TYPE_CONST);
    RELXILL_ERROR("unknown ionization gradient type", status);
  }

  if (is_debug_run()) {
    write_binned_data_to_file("test_ion_grad_relxill.dat", ion->r, ion->lxi, ion->nbins);
  }

  if (*status != EXIT_SUCCESS) {
    RELXILL_ERROR("calculating the ionization gradient failed due to previous error", status);
  }

  return ion;
}

void free_ion_grad(ion_grad *ion) {

  if (ion != nullptr) {
    if (ion->r != nullptr) {
      free(ion->r);
    }
    if (ion->lxi != nullptr) {
      free(ion->lxi);
    }
    if (ion->fx != nullptr) {
      free(ion->fx);
    }
    free(ion);
  }
}
