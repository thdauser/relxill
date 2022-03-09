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
#include "Relbase.h"

extern "C"{
#include "writeOutfiles.h"
}

// calculate the log(xi) for given density and emissivity
// note: delta_inc is taking the correction into account that xillver
//       is calculated for a fixed delta_inc=PI/4 (see Ingram+19)
static double cal_lxi(double dens, double emis, double delta_inc) {
  double lxi_corrfactor = cos(M_PI / 4) / cos(delta_inc);
  return log10(4.0 * M_PI * emis / dens * lxi_corrfactor);
}

/**
  * adjust "values" to be within the limits of the xillver table
  **/
void adjust_parameter_to_be_within_xillver_bounds(double* values, int nzones, double tab_min, double tab_max) {
  for (int ii=0; ii<nzones; ii++){
    if (values[ii] < tab_min) {
      values[ii] = tab_min;
    } else if (values[ii] > tab_max) {
      values[ii] = tab_max;
    }
  }
}

// determine the maximal ionization
static double cal_lxi_max_ss73(double *re, double *emis, double* delinc, int nr, double rin) {

  // we use the same definition as Adam with r_peak = (11/9)^2 rin to be consistent (does not matter much)
  double rad_max_lxi = pow((11. / 9.), 2) * rin;

  // radial AD grid is sorted descending (!)
  int kk = inv_binary_search(re, nr, rad_max_lxi);
  double interp = (rad_max_lxi - re[kk + 1]) / (re[kk] - re[kk + 1]);

  double emis_max_lxi = interp_lin_1d(interp, emis[kk + 1], emis[kk]);
  double delinc_max_lxi = interp_lin_1d(interp, delinc[kk + 1], delinc[kk]);

  double lxi_max = cal_lxi(density_ss73_zone_a(rad_max_lxi, rin), emis_max_lxi, delinc_max_lxi);

  return lxi_max;
}


void IonGradient::calc_ion_grad_alpha(emisProfile* emis_profile, double param_xlxi0, double param_density) {

  double rin = radial_grid.radius[0];

  irradiating_flux = new double[m_nzones];
  auto del_inc = new double[m_nzones];

  int status = EXIT_SUCCESS;

  assert(emis_profile->del_inc != nullptr);
  inv_rebin_mean(emis_profile->re, emis_profile->emis, emis_profile->nr, m_rmean, irradiating_flux, m_nzones, &status);
  inv_rebin_mean(emis_profile->re, emis_profile->del_inc, emis_profile->nr, m_rmean, del_inc, m_nzones, &status);
  inv_rebin_mean(emis_profile->re, emis_profile->del_emit, emis_profile->nr, m_rmean, del_emit, m_nzones, &status);

  // calculate the maximal ionization assuming r^-3 and SS73 alpha disk (see Ingram+19)
  double lxi_max = cal_lxi_max_ss73(emis_profile->re, emis_profile->emis, emis_profile->del_inc, emis_profile->nr, rin);

  // the maximal ionization is given as input parameter, so we need to normalize our calculation by this value
  double fac_lxi_norm = param_xlxi0 - lxi_max; // subtraction instead of division because of the log

  /** calculate the density for a  stress-free inner boundary condition, i.e., R0=rin in SS73)  **/
  for (int ii = 0; ii < m_nzones; ii++) {
    double density_normed = density_ss73_zone_a(m_rmean[ii], rin);
    dens[ii] = log10(density_normed) + param_density;  //addition as quantity is logarithm

    // now we can use the emissivity to calculate the ionization
    // (note we need a density which is normalized to dens(rin)=1 )
    lxi[ii] = cal_lxi(density_normed, irradiating_flux[ii], del_inc[ii]);
    lxi[ii] += fac_lxi_norm;

  }



}


void IonGradient::calc_ion_grad_pl(double xlxi0, double xindex, double inputval_dens){
  for (int ii = 0; ii < m_nzones; ii++) {
    lxi[ii] = (exp(xlxi0))
        * pow((m_rmean[ii] / m_rmean[0]), -1.0 * xindex);  // TODO: check if we need to subtract xlxi_tab_min here
        lxi[ii] = log(lxi[ii]);

    dens[ii] = inputval_dens;
  }
}


void IonGradient::calculate(emisProfile* emis_profile, xillParam *xill_param) {

  if (m_ion_grad_type == ION_GRAD_TYPE_PL) {
    calc_ion_grad_pl(xill_param->lxi, xill_param->iongrad_index, xill_param->dens);

  } else if (m_ion_grad_type == ION_GRAD_TYPE_ALPHA) {
    calc_ion_grad_alpha(emis_profile, xill_param->lxi, xill_param->dens);

  } else if (m_ion_grad_type == ION_GRAD_TYPE_CONST) {
    for (int ii = 0; ii < m_nzones; ii++) {
      lxi[ii] = xill_param->lxi;
      dens[ii] = xill_param->dens;
    }

  } else {
    printf(" *** ionization type with number %i not implemented \n", m_ion_grad_type);
    printf("     choose either %i for the PL, %i for the ALPHA-disk, or %i for constant\n",
           ION_GRAD_TYPE_PL, ION_GRAD_TYPE_ALPHA, ION_GRAD_TYPE_CONST);
    throw std::exception();
  }

  // need to make sure lxi and logn are within the xillver table values //  TODO: Need to define this automatically
  adjust_parameter_to_be_within_xillver_bounds(lxi, m_nzones, 0.0, 4.7);
//  adjust_parameter_to_be_within_xillver_bounds(dens, m_nzones, 15, 19);

  if  (shouldOutfilesBeWritten()) {
    IonGradient::write_to_file("__relxillOutput_iongrad.dat");
  }

}

void IonGradient::write_to_file(const char* fout) {

  FILE *fp = fopen(fout, "w+");

  fprintf(fp, "# rlo [R_g]\t rhi [R_g] \t log(xi) \t log(N) \t delta_emit \n");
    for (int ii = 0; ii < m_nzones; ii++) {
      fprintf(fp, " %e \t %e \t %e \t %e \t %e \n", radial_grid.radius[ii], radial_grid.radius[ii + 1],
              lxi[ii], dens[ii], del_emit[ii]);
    }
    fclose_errormsg(fp, fout);
}


const double* RadialGrid::calculate_radial_grid(double rmin, double rmax, const int nzones, double h) {

  auto rgrid =  new double[nzones + 1];

  if (nzones == 1) {
    rgrid[0] = rmin;
    rgrid[1] = rmax;
  } else {

    double r_transition = rmin;
    int indr = 0;

    if (h > rmin) { // if h > rmin we choose a log grid for r<h
      r_transition = h;

      get_log_grid(rgrid, nzones + 1, rmin, rmax);
      indr = binary_search(rgrid, nzones + 1, r_transition);

      r_transition = rgrid[indr];
    }

    if (indr < nzones) {
      double rlo = r_transition;
      double rhi = rmax; // radius[nzones];
      // add 1/r for larger radii
      int ii;
      for (ii = indr; ii < nzones + 1; ii++) {
        rgrid[ii] = 1.0 * (ii - indr) / (nzones - indr) * (1.0 / rhi - 1.0 / rlo) + 1.0 / rlo;
        rgrid[ii] = fabs(1.0 / rgrid[ii]);
      }

    }

  }

  return rgrid;
}
