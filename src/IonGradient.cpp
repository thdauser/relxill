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
#include "IonGradient.h"

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
void adjust_parameter_to_be_within_xillver_bounds(double *values, int nzones, double tab_min, double tab_max) {
  for (int ii = 0; ii < nzones; ii++) {
    if (values[ii] < tab_min) {
      values[ii] = tab_min;
    } else if (values[ii] > tab_max) {
      values[ii] = tab_max;
    }
  }
}

// we use the same definition as Adam with r_peak = (11/9)^2 rin to be consistent (does not matter much)
static double get_radius_lxi_max_ss73(double rin) {
  return pow((11. / 9.), 2) * rin;
}

/**
 * @brief calculate the emissivit profile (emissivity, del_inc, del_emit)for the give radius
 * from interpolation of the emissivity profile
 * @param rad            [input]
 * @param emis_profile   [input]
 */
static emisProfile *get_emis_at_rad(double rad, const emisProfile &emis_profile) {

  // we only need one bin of the emissivity profile
  const int nbins = 1;
  const int ibin = 0;

  int status = EXIT_SUCCESS;
  auto emis_struct = new_emisProfile(&rad, nbins, &status);

  // radial AD grid is sorted descending (!)
  const int kk = inv_binary_search(emis_profile.re, emis_profile.nr, rad);
  const double interp = (rad - emis_profile.re[kk + 1]) / (emis_profile.re[kk] - emis_profile.re[kk + 1]);

  emis_struct->emis[ibin] = interp_lin_1d(interp, emis_profile.emis[kk + 1], emis_profile.emis[kk]);
  emis_struct->del_inc[ibin] = interp_lin_1d(interp, emis_profile.del_inc[kk + 1], emis_profile.del_inc[kk]);
  emis_struct->del_emit[ibin] = interp_lin_1d(interp, emis_profile.del_emit[kk + 1], emis_profile.del_emit[kk]);

  return emis_struct;
}

/**
 * @brief calculates the emissivity profile at the radius of maximal ionization [only 1 value]
 * @param rel_param
 * @param emis_profile
 * @return
 */
static emisProfile *get_emis_profile_lxi_max(const relParam &rel_param, const emisProfile &emis_profile) {

  const double radius_lxi_max = get_radius_lxi_max_ss73(rel_param.rin);
  auto emis_profile_rad = get_emis_at_rad(radius_lxi_max, emis_profile);

  return emis_profile_rad;
}

// determine the maximal ionization
static double cal_lxi_max_ss73(const emisProfile &emis_profile, double rin) {

  const double rad_lxi = get_radius_lxi_max_ss73(rin);
  auto emis_profile_rad = get_emis_at_rad(rad_lxi, emis_profile);  // 1dim arrays

  const double lxi_max = cal_lxi(density_ss73_zone_a(rad_lxi, rin),
                                 emis_profile_rad->emis[0],
                                 emis_profile_rad->del_inc[0]);

  free_emisProfile(emis_profile_rad);

  return lxi_max;
}

/**
 * @brief: calculate the ionization gradient for the alpha type models
 * @param emis_profile
 * @param param_xlxi0
 * @param param_density
 */
void IonGradient::calc_ion_grad_alpha(const emisProfile &emis_profile, double param_xlxi0, double param_density) {

  const double rin = radial_grid.radius[0];

  irradiating_flux = new double[m_nzones];
  auto del_inc = new double[m_nzones];

  int status = EXIT_SUCCESS;

  assert(emis_profile.del_inc != nullptr);
  inv_rebin_mean(emis_profile.re, emis_profile.emis, emis_profile.nr, m_rmean, irradiating_flux, m_nzones, &status);
  inv_rebin_mean(emis_profile.re, emis_profile.del_inc, emis_profile.nr, m_rmean, del_inc, m_nzones, &status);

  if (status != EXIT_SUCCESS) {
    throw std::exception();
  }

  // calculate the maximal ionization assuming r^-3 and SS73 alpha disk (see Ingram+19)
  const double lxi_max = cal_lxi_max_ss73(emis_profile, rin);

  // the maximal ionization is given as input parameter, so we need to normalize our calculation by this value
  const double fac_lxi_norm = param_xlxi0 - lxi_max; // subtraction instead of division because of the log

  /** calculate the density for a  stress-free inner boundary condition, i.e., R0=rin in SS73)  **/
  for (int ii = 0; ii < m_nzones; ii++) {
    const double density_normed = density_ss73_zone_a(m_rmean[ii], rin);
    dens[ii] = log10(density_normed) + param_density;  //addition as quantity is logarithm

    // now we can use the emissivity to calculate the ionization
    // (note we need a density which is normalized to dens(rin)=1 )
    lxi[ii] = cal_lxi(density_normed, irradiating_flux[ii], del_inc[ii]);
    lxi[ii] += fac_lxi_norm;

  }

  delete[] del_inc;

}


void IonGradient::calc_ion_grad_pl(double xlxi0, double xindex, double inputval_dens) {
  for (int ii = 0; ii < m_nzones; ii++) {
    lxi[ii] = (exp(xlxi0))
        * pow((m_rmean[ii] / m_rmean[0]), -1.0 * xindex);
    lxi[ii] = log(lxi[ii]);

    dens[ii] = inputval_dens;
  }
}

void IonGradient::calc_energy_shift_from_source_to_disk(const relParam *rel_param) const {
  assert (del_emit != nullptr);

  for (int ii = 0; ii < m_nzones; ii++) {
    if (rel_param->emis_type == EMIS_TYPE_LP) {
      m_energy_shift_source_disk[ii] = energy_shift_source_disk(rel_param, m_rmean[ii], del_emit[ii]);
    } else {
      m_energy_shift_source_disk[ii] = 1.0;
    }
  }
}

xillTableParam **IonGradient::get_xill_param_zone(const xillTableParam *primary_source_xilltab_params) const {

  auto xill_param_zone = new xillTableParam *[m_nzones];

  for (int ii = 0; ii < m_nzones; ii++) {
    xill_param_zone[ii] = new xillTableParam;
    (*xill_param_zone[ii]) = (*primary_source_xilltab_params);  // shallow copy (enough for the parameter structure)

    // set xillver parameters for the given zone
    xill_param_zone[ii]->ect = primary_source_xilltab_params->ect * m_energy_shift_source_disk[ii];
    xill_param_zone[ii]->lxi = lxi[ii];
    xill_param_zone[ii]->dens = dens[ii];
  }

  return xill_param_zone;
}

double calc_flux_rin_cgs(double l_source, double height, double Rin, double mass) {

  const double CONST_rgIncm = 1.4822e5;  // G / c^2 * (M/Msolar)

  const double dRay_rg = sqrt(pow(height, 2) + pow(Rin, 2));
  const double dRay_cm = dRay_rg * CONST_rgIncm * mass;

  return l_source / (4 * M_PI * pow(dRay_cm, 2));
}

/**
 * @brief: for a given primary source configuration (mass, distance, observed flux), we calculate the flux at the
 * given radius on the disk, using no GR for the flux conversion from source to disk (i.e., no energy shifts and
 * and straight photons paths)
 * @param photon_fate_fractions
 * @param primary_source
 * @param rad
 * @return
 */
double calc_flux_disk_newton_cgs(const lpReflFrac *photon_fate_fractions,
                                 const PrimarySourceParameters &primary_source,
                                 double rad) {

  const auto &rel_param = (*primary_source.rel_param());

  // calculate the source luminosity for boost to get a higher/lower source luminosity
  // note: at add_primary_component this still needs to be added to the normalization of the reflection spectrum
  double l_source_cgs = primary_source.luminosity_source_cgs((*photon_fate_fractions));
  const double boost = primary_source.get_boost_parameter(photon_fate_fractions);
  l_source_cgs *= boost;

  return calc_flux_rin_cgs(l_source_cgs, rel_param.height, rad, primary_source.mass_msolar());
}

/**
 * @brief calculate the flux on the disk at the radius given by the emissivity profile
 * (require the emissivity profile to be only of length=1!)
 * @param photon_fate_fractions
 * @param emis_profile_lxi_max  (is only allowed to contain one value)
 * @param rel_param
 * @return
 */
double calc_flux_disk_cgs(const lpReflFrac *photon_fate_fractions,
                          const emisProfile &emis_profile_lxi_max,
                          const PrimarySourceParameters &primary_source) {
  assert(emis_profile_lxi_max.nr == 1);

  const double radius_lxi_max = emis_profile_lxi_max.re[0];
  const auto &rel_param = (*primary_source.rel_param());

  // calculate the flux boost wrt to Newton from the given emissivity profile, as it is defined such that for large
  // radii (i.e., neglecting GR effects it coincides with the Newtonian version normalized as \int \emis dA = 1)
  double const boost_newton_to_gr =
      photon_fate_fractions->refl_frac * emis_profile_lxi_max.emis[0]
          / calc_lp_emissivity_newton(rel_param.height, radius_lxi_max);

  const double flux_rad_lxi_max_cgs =
      calc_flux_disk_newton_cgs(photon_fate_fractions, primary_source, radius_lxi_max)
          * boost_newton_to_gr;

  return flux_rad_lxi_max_cgs;
}

/**
 * @brief: calculates the value of maximal value of log(xi) for an alpha disk density gradient for a given
 * density at the Rin (note that the maximal value of log(xi) is not at Rin, but (11/9)^2*Rin)
 * @param emis_profile
 * @param primary_source
 * @param log_density_rin
 * @return
 */
double IonGradient::calculate_lxi_max_from_distance(const emisProfile &emis_profile,
                                                    const PrimarySourceParameters &primary_source,
                                                    double log_density_rin) {

  const auto &rel_param = (*primary_source.rel_param());

  auto emis_profile_lxi_max = get_emis_profile_lxi_max(rel_param, emis_profile);
  const double radius_lxi_max = emis_profile_lxi_max->re[0];

  const double flux_rad_lxi_max_cgs =
      calc_flux_disk_cgs(emis_profile.photon_fate_fractions, (*emis_profile_lxi_max), primary_source);

  // the routine "density_ss73_zone_a" returns a density of 1 at rin
  const double density_lxi_max = density_ss73_zone_a(radius_lxi_max, rel_param.rin) * pow(10, log_density_rin);

  const double lxi_max = cal_lxi(density_lxi_max, flux_rad_lxi_max_cgs, emis_profile_lxi_max->del_inc[0]);

  free_emisProfile(emis_profile_lxi_max);

  return lxi_max;
}

/**
 * @brief calculates the physical emissivity (meaning in units of erg/cm2/s) for the given source configuration
 * needs the irradiating_flux to be calculated
 * @param primary_source
 * @param emis_profile
 */
void IonGradient::calculate_physical_emissivity(const PrimarySourceParameters &primary_source,
                                                const emisProfile &emis_profile) {

  if (irradiating_flux == nullptr) {
    return;
  }

  const auto &rel_param = (*primary_source.rel_param());

  auto emis_profile_lxi_max = get_emis_profile_lxi_max(rel_param, emis_profile);
  const double radius_lxi_max = emis_profile_lxi_max->re[0];

  const double flux_rad_lxi_max_cgs =
      calc_flux_disk_cgs(emis_profile.photon_fate_fractions, (*emis_profile_lxi_max), primary_source);

  const double physical_emis_normalization_factor = flux_rad_lxi_max_cgs / emis_profile_lxi_max->emis[0];

  m_physical_irrad_flux = new double[m_nzones];
  for (int ii = 0; ii < m_nzones; ii++) {
    m_physical_irrad_flux[ii] = irradiating_flux[ii] * physical_emis_normalization_factor;
  }

}


void IonGradient::calculate_gradient(const emisProfile &emis_profile,
                                     const PrimarySourceParameters &primary_source_params) {

  set_del_emit_for_each_zone(emis_profile);

  auto rel_param = primary_source_params.rel_param();
  auto xill_param = primary_source_params.xilltab_param();

  calc_energy_shift_from_source_to_disk(rel_param);

  if (m_ion_grad_type == ION_GRAD_TYPE_PL) {
    calc_ion_grad_pl(xill_param->lxi, m_ion_grad_index, xill_param->dens);

  } else if (m_ion_grad_type == ION_GRAD_TYPE_ALPHA) {

    const double lxi_rin = (is_alpha_model(rel_param->model_type))
                           ? IonGradient::calculate_lxi_max_from_distance(emis_profile,
                                                                          primary_source_params,
                                                                          xill_param->dens)
                           : xill_param->lxi;

    calc_ion_grad_alpha(emis_profile, lxi_rin, xill_param->dens);

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

  calculate_physical_emissivity(primary_source_params, emis_profile);

  if (shouldOutfilesBeWritten()) {
    if (is_alpha_model(rel_param->model_type)) {
      calculate_physical_emissivity(primary_source_params, emis_profile);
    }
    IonGradient::write_to_file("__relxillOutput_iongrad.dat");
  }

  // need to make sure lxi and logn are within the xillver table values //  TODO: Need to define this automatically
  adjust_parameter_to_be_within_xillver_bounds(lxi, m_nzones, 0.0, 4.7);

}

/**
 * @brief: write the ionization gradient to file; for the alpha type gradient, als the irradiating flux
 * is written
 * @param fout
 */
void IonGradient::write_to_file(const char *fout) const {

  FILE *fp = fopen(fout, "w+");

  if (irradiating_flux != nullptr) {
    fprintf(fp, "# rlo [R_g]\t rhi [R_g] \t log(xi) \t log(N) \t delta_emit \n F_x\n");

    // if we calculated the physical flux, use this for the output
    const double *irrad_flux = (m_physical_irrad_flux != nullptr) ? m_physical_irrad_flux : irradiating_flux;

    for (int ii = 0; ii < m_nzones; ii++) {
      fprintf(fp, " %e \t %e \t %e \t %e \t %e \t %e \n", radial_grid.radius[ii], radial_grid.radius[ii + 1],
              lxi[ii], dens[ii], del_emit[ii], irrad_flux[ii]);
    }
  } else {
    fprintf(fp, "# rlo [R_g]\t rhi [R_g] \t log(xi) \t log(N) \t delta_emit \n ");
    for (int ii = 0; ii < m_nzones; ii++) {
      fprintf(fp, " %e \t %e \t %e \t %e \t %e  \n", radial_grid.radius[ii], radial_grid.radius[ii + 1],
              lxi[ii], dens[ii], del_emit[ii]);
    }
  }
  fclose_errormsg(fp, fout);
}

void IonGradient::set_del_emit_for_each_zone(const emisProfile &emis_profile) {
  int status = EXIT_SUCCESS;
  if (del_emit == nullptr) {
    del_emit = new double[m_nzones];
  }
  inv_rebin_mean(emis_profile.re, emis_profile.del_emit, emis_profile.nr, m_rmean, del_emit, m_nzones, &status);
  if (status != EXIT_SUCCESS) {
    throw std::exception();
  }
}

/** @brief calculate the Ecut/kTe value as seen by the disk zone (if nzones>1)
 *        - for nzones=1 the value at the primary source is returned
 *        - the cutoff is calculated for the (linear) middle of the radial zone
 * @param rel_param
 * @param ecut_primary [ecut/kTe value at the primary source]
 * @param izone [zone index]
 * @return
 */
double IonGradient::get_ecut_disk_zone(const relParam *rel_param, double ecut_primary, int izone) const {

  if (radial_grid.num_zones == 1) {
    return ecut_primary; // TODO: not obvious what "ecut0" is (it is the input value, from the model fitting)
  } else {
    double rzone = 0.5 * (radial_grid.radius[izone] + radial_grid.radius[izone + 1]);

    if (del_emit == nullptr) { // del_emit relevant if beta!=0 (doppler boosting)
      printf(" *** error in IonGradient: ionization gradient not calculated \n");
      printf("       can not calculate Ecut/kTe on the disk \n\n");
      throw std::exception();
    }

    return ecut_primary * gi_potential_lp(rzone, rel_param->a, rel_param->height, rel_param->beta, del_emit[izone]);
  }
}

const double *RadialGrid::calculate_radial_grid(double rmin, double rmax, const int nzones, double h) {

  auto rgrid = new double[nzones + 1];

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
