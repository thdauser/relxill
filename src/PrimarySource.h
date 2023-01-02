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

#ifndef RELXILL_SRC_PRIMARYSOURCE_H_
#define RELXILL_SRC_PRIMARYSOURCE_H_

#include <string>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <utility>

#include "ModelParams.h"
#include "Relphysics.h"
#include "Relbase.h"

extern "C" {
#include "common.h"
#include "xilltable.h"
}

#include <vector>



////////////////////////////////////////////
class PrimarySourceParameters {
  ////////////////////////////////////////////

 public:

  explicit PrimarySourceParameters(const ModelParams &inp_param) :
      m_inp_param{inp_param},
      m_rel_param{get_rel_params(m_inp_param)},
      m_refl_frac{inp_param.get_otherwise_default(XPar::refl_frac, 0)},
      m_normalization_factor{inp_param.get_otherwise_default(XPar::norm_factor, 1)},
      m_distance{inp_param.get_otherwise_default(XPar::distance, 0)},
      m_interpret_reflfrac_as_boost{
          static_cast<int>(lround(inp_param.get_otherwise_default(XPar::switch_switch_reflfrac_boost, 0)))} {

    m_energy_shift_source_observer = (m_rel_param != nullptr && m_rel_param->emis_type == EMIS_TYPE_LP) ?
                                     energy_shift_source_obs(m_rel_param) : 1.0;

    m_xilltab_param = m_get_xillver_params_primary_source(m_rel_param, m_energy_shift_source_observer);

  }

  ~PrimarySourceParameters() {
    delete m_rel_param;
    delete m_xilltab_param;
  }

  /**
   * @brief: calculate the flux boost from the source spectrum to the observered spectrum
   * @param f_inf
   * @return
   */
  double luminosity_source_cgs(double f_inf) {
    const double lum_observed = 4 * M_PI * m_distance * m_distance * m_normalization_factor;

    // add energy shift and flux boost
    const double flux_boost_source_observer =
        pow(m_energy_shift_source_observer, m_rel_param->gamma) * f_inf / 0.5 * doppler_factor_source_obs(m_rel_param);

    return lum_observed / flux_boost_source_observer;
  }

  // return full xilltab and rel_param structures (TODO: refactoring methods such that this is not needed, at least for rel_param?)
  [[nodiscard]] const xillTableParam *xilltab_param() const {
    return m_xilltab_param;
  }
  [[nodiscard]] const relParam *rel_param() const {
    return m_rel_param;
  }

  int model_type() const {
    return m_xilltab_param->model_type;
  }
  int prim_type() const {
    return m_xilltab_param->prim_type;
  }
  int emis_type() const { // need to make sure it works for xillver models as well
    return (m_rel_param != nullptr) ? m_rel_param->emis_type : EMIS_TYPE_NONE;
  }

  xillTableParam *get_xillver_params_primary_source() {
    return m_get_xillver_params_primary_source(m_rel_param, m_energy_shift_source_observer);
  }

  double normalization_factor() const {
    return m_normalization_factor;
  }
  double distance() const {
    return m_distance;
  }
  double refl_frac() const {
    return m_refl_frac;
  }
  int interpret_reflfrac_as_boost() const {
    return m_interpret_reflfrac_as_boost;
  }
  double energy_shift_source_observer() const {
    return m_energy_shift_source_observer;
  }

  // functions for parameters we need
  [[nodiscard]] double gam() const {
    return m_xilltab_param->gam;
  }
  [[nodiscard]] double beta() const {
    return m_rel_param->beta;
  }
  [[nodiscard]] double ect() const {
    return m_xilltab_param->ect;
  }

 private:

  double m_normalization_factor;
  double m_distance;
  double m_refl_frac;
  int m_interpret_reflfrac_as_boost;
  double m_energy_shift_source_observer;

  ModelParams m_inp_param;
  relParam *m_rel_param = nullptr;
  xillTableParam *m_xilltab_param = nullptr;

  xillTableParam *m_get_xillver_params_primary_source(relParam *_rel_param, double _energy_shift_source_observer) {
    int status = EXIT_SUCCESS;
    auto xill_param = get_xill_params(m_inp_param);
    auto xilltab_param = get_xilltab_param(xill_param, &status);
    delete xill_param;

    // special case, for the LP model with the Ecut spectrum: cutoff energy is given in the observer frame
    // -> convert it such that ecut is also given in the source frame
    if (_rel_param != nullptr && _rel_param->emis_type == EMIS_TYPE_LP && xilltab_param->prim_type == PRIM_SPEC_ECUT) {
      xilltab_param->ect /= _energy_shift_source_observer;
    }

    return xilltab_param;
  }
};

// TODO:
//  - properly include the reflection fraction and a nice error message
////////////////////////////////////////////
class PrimarySource {
  ////////////////////////////////////////////

 public:
  PrimarySource(const ModelParams &_model_params, RelSysPar *sys_par, const XspecSpectrum &_xspec_spec) :
      source_parameters{_model_params},
      m_struct_refl_frac((sys_par == nullptr) ? nullptr : sys_par->emis->photon_fate_fractions),
      m_prime_spec_observer{calculate_observed_primary_spectrum(_xspec_spec, source_parameters)} {

    //  m_prime_spec_xillver_norm_factor = 1. / calculate_source_prime_spec_xillver_norm();
    // source to obs normalization
    //
  }

  double calculate_source_prime_spec_xillver_norm() const {
    int status = EXIT_SUCCESS;
    // need to create a specific energy grid for the primary component normalization to
    // fulfill the XILLVER NORM condition (Dauser+2016)
    EnerGrid *egrid = get_coarse_xillver_energrid(&status);
    auto prim_spec_source = new double[egrid->nbins];
    calc_primary_spectrum(prim_spec_source, egrid->ener, egrid->nbins,
                          source_parameters.xilltab_param(), &status);

    // calculate the normalization at the primary source
    double const
        norm_xillver_prim_spec = calcNormWrtXillverTableSpec(prim_spec_source, egrid->ener, egrid->nbins, &status);
    delete[] prim_spec_source;

    return norm_xillver_prim_spec;
  }

  //
  static Spectrum calculate_observed_primary_spectrum(const XspecSpectrum &_xspec_spec,
                                                      const PrimarySourceParameters &_parameters) {
    int status = EXIT_SUCCESS;
    auto spec = Spectrum(_xspec_spec.energy, _xspec_spec.num_flux_bins());
    calc_primary_spectrum(spec.flux, spec.energy(), spec.num_flux_bins,
                          _parameters.xilltab_param(), &status,
                          _parameters.energy_shift_source_observer());
    return spec;
  }

  void print_reflection_strength(const XspecSpectrum &spectrum, const double *pl_flux) const;
  void add_primary_spectrum(const XspecSpectrum &spectrum);

  PrimarySourceParameters source_parameters;

 private:

  lpReflFrac *m_struct_refl_frac = nullptr;
  Spectrum m_prime_spec_observer;  // not normalized
  // double m_prime_spec_xillver_norm_factor = 1.0;

  // calculate the flux of the primary spectrum in ergs/cm2/sec
  double get_normalized_primary_spectrum_flux_in_ergs() const {
    int status = EXIT_SUCCESS;
    EnerGrid *egrid = get_coarse_xillver_energrid(&status);
    auto prim_spec_source = new double[egrid->nbins];
    double *pl_flux = calc_normalized_primary_spectrum(egrid->ener,
                                                       egrid->nbins,
                                                       source_parameters.rel_param(),
                                                       source_parameters.xilltab_param(),
                                                       &status);

    double ener_flux = 0.0;
    for (int ii = 0; ii < egrid->nbins; ii++) {
      if (egrid->ener[ii] >= EMIN_XILLVER_NORMALIZATION && egrid->ener[ii + 1] < EMAX_XILLVER_NORMALIZATION) {
        ener_flux += pl_flux[ii] * 0.5 * (egrid->ener[ii] + egrid->ener[ii + 1]);
      }
    }

    delete[] prim_spec_source;

    return ener_flux * CONVERT_KEV2ERG;
  }

};

#endif //RELXILL_SRC_PRIMARYSOURCE_H_
