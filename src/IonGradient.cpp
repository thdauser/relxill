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
static double cal_lxi_max_ss73(double *re, double *emis, int nr, double rin) {

  // we use the same definition as Adam with r_peak = (11/9)^2 rin to be consistent (does not matter much)
  double rad_max_lxi = pow((11. / 9.), 2) * rin;

  // radial AD grid is sorted descending (!)
  int kk = inv_binary_search(re, nr, rad_max_lxi);
  double interp = (rad_max_lxi - re[kk + 1]) / (re[kk] - re[kk + 1]);

  double emis_max_lxi = interp_lin_1d(interp, emis[kk + 1], emis[kk]);
  double lxi_max = cal_lxi(density_ss73_zone_a(rad_max_lxi, rin), emis_max_lxi);

  return lxi_max;
}


void IonGradient::calc_ion_grad_alpha(relParam *rel_param, double param_xlxi0, double param_density) {

  double rin = m_radius[0];

  auto emis_zones = new double[m_nzones];
  auto del_inc = new double[m_nzones];

  int status = EXIT_SUCCESS;

  // we need the emissivity profile (should be cached, so no extra effort required here)
  RelSysPar *sysPar = get_system_parameters(rel_param, &status);
  emisProfile *emis_profile = sysPar->emis;

  assert(emis_profile->del_inc != nullptr);
  inv_rebin_mean(emis_profile->re, emis_profile->emis, sysPar->nr, m_rmean, emis_zones, m_nzones, &status);
  inv_rebin_mean(emis_profile->re, emis_profile->del_inc, sysPar->nr, m_rmean, del_inc, m_nzones, &status);
  inv_rebin_mean(emis_profile->re, emis_profile->del_emit, sysPar->nr, m_rmean, del_emit, m_nzones, &status);

  // calculate the maximal ionization assuming r^-3 and SS73 alpha disk
  double lxi_max = cal_lxi_max_ss73(emis_profile->re, emis_profile->emis, emis_profile->nr, rin);

  // the maximal ionization is given as input parameter, so we need to normalize our calculation by this value
  double fac_lxi_norm = param_xlxi0 - lxi_max; // subtraction instead of division because of the log

  /** calculate the density for a  stress-free inner boundary condition, i.e., R0=rin in SS73)  **/
  for (int ii = 0; ii < m_nzones; ii++) {
    dens[ii] = density_ss73_zone_a(m_rmean[ii], rin);  //addition as quantity is logarithm
    dens[ii] /= dens[0];

    // now we can use the emissivity to calculate the ionization
    lxi[ii] = cal_lxi(dens[ii], emis_zones[ii]) + fac_lxi_norm;
    lxi[ii] += log10(cos(M_PI / 4) / cos(del_inc[ii]));

    dens[ii] = log10(dens[ii]) + param_density;  //addition as quantity is logarithm


  }



}


void IonGradient::calc_ion_grad_pl(double xlxi0, double xindex){
  for (int ii = 0; ii < m_nzones; ii++) {
    lxi[ii] = (exp(xlxi0))
        * pow((m_rmean[ii] / m_rmean[0]), -1.0 * xindex);  // TODO: check if we need to subtract xlxi_tab_min here
        lxi[ii] = log(lxi[ii]);
  }
}


void IonGradient::calculate(relParam *rel_param, xillParam *xill_param) {

  if (m_ion_grad_type == ION_GRAD_TYPE_PL) {
    calc_ion_grad_pl(xill_param->lxi, xill_param->ion_grad_index);

  } else if (m_ion_grad_type == ION_GRAD_TYPE_ALPHA) {
    calc_ion_grad_alpha(rel_param, xill_param->lxi, xill_param->dens);

  } else if (m_ion_grad_type == ION_GRAD_TYPE_CONST) {
    for (int ii = 0; ii < m_nzones; ii++) {
      lxi[ii] = xill_param->lxi;
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
    IonGradient::write_to_file("test_ion_grad_relxill.dat");
  }

}

void IonGradient::write_to_file(const char* fout) {

  FILE *fp = fopen(fout, "w+");

  fprintf(fp, "# rlo \t rhi \t logxi \t density \t delta_emit \n");
    for (int ii = 0; ii < m_nzones; ii++) {
      fprintf(fp, " %e \t %e \t %e \t %e \t %e \n", m_radius[ii], m_radius[ii + 1],
              lxi[ii], dens[ii], del_emit[ii]);
    }
    fclose_errormsg(fp, fout);
}

