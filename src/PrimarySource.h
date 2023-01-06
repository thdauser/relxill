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
      m_interpret_reflfrac_as_boost{
          static_cast<int>(lround(inp_param.get_otherwise_default(XPar::switch_switch_reflfrac_boost, 0)))} {
    m_energy_shift_source_observer = (m_rel_param != nullptr && m_rel_param->emis_type == EMIS_TYPE_LP) ?
                                     energy_shift_source_obs(m_rel_param) : 1.0;

    auto xill_param = get_xill_params(m_inp_param);
    m_refl_frac = xill_param->refl_frac;
    m_norm_flux_cgs = xill_param->norm_flux_cgs;
    m_distance_kpc = xill_param->distance;
    m_mass_msolar = xill_param->mass_msolar;
    delete xill_param;

    m_xilltab_param = m_get_xilltab_params_primary_source(m_rel_param, m_energy_shift_source_observer);
  }

  ~PrimarySourceParameters() {
    delete m_rel_param;
    delete m_xilltab_param;
  }

  // delete the copy constructor and copy assignment operator [deleted as it crashes gitlab]
  // PrimarySourceParameters(const PrimarySourceParameters &) = delete;
  // PrimarySourceParameters &operator=(const PrimarySourceParameters &) = delete;

  /**
   * @brief calculates the boost of the flux (i.e. also the spectrum normalization) from source to observer
   * due to (1) the energy shift from source to observer and (2) a potential velocity of the primary source
   * @param lp_refl_frac
   * @return
   */
  double flux_boost_source_to_observer(const lpReflFrac &lp_refl_frac) const {
    double prime_spec_factors_source_to_observer =
        lp_refl_frac.f_inf_rest / 0.5 * pow(m_energy_shift_source_observer, PrimarySourceParameters::gam());

    // flux boost of primary radiation taking into account here (therefore we need f_inf_rest above)
    if (PrimarySourceParameters::beta() > 1e-4) {
      prime_spec_factors_source_to_observer *= pow(doppler_factor_source_obs(m_rel_param), 2);
    }
    return prime_spec_factors_source_to_observer;
  }

  /**
   * @brief: calculate the flux boost from the source spectrum to the observed spectrum
   * @param f_inf
   * @return luminosity source [ergs/s]
   */
  double luminosity_source_cgs(const lpReflFrac &lp_refl_frac, double boost = 1.0) const {

    const double CONST_cm2kpc = 3.2407792700054E-22;

    const double distance_cm = m_distance_kpc / CONST_cm2kpc;
    const double lum_observed = 4 * M_PI * distance_cm * distance_cm * m_norm_flux_cgs;

    // add energy shift and flux boost
    return lum_observed / flux_boost_source_to_observer(lp_refl_frac);
  }

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
    return m_get_xilltab_params_primary_source(m_rel_param, m_energy_shift_source_observer);
  }

  double norm_flux_cgs() const {
    return m_norm_flux_cgs;
  }
  double distance() const {
    return m_distance_kpc;
  }

  double mass_msolar() const {
    return m_mass_msolar;
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

  /**
   * @brief: get the absolute value of the boost parameter for the current configuration
   **/
  double get_boost_parameter(lpReflFrac *lp_refl_frac) const {

    if (lp_refl_frac == nullptr) {
      printf(
          " *** error: can not calculate the boost parameter if the reflection fraction information is not available\n");
      throw std::exception();
    }

    if (m_interpret_reflfrac_as_boost) {
      return fabs(m_refl_frac);
    } else {
      return fabs(m_refl_frac) / lp_refl_frac->refl_frac;
    }

  }

 private:

  double m_norm_flux_cgs;
  double m_distance_kpc;
  double m_refl_frac;
  double m_mass_msolar;
  int m_interpret_reflfrac_as_boost;
  double m_energy_shift_source_observer;

  ModelParams m_inp_param;
  relParam *m_rel_param;
  xillTableParam *m_xilltab_param = nullptr;

  xillTableParam *m_get_xilltab_params_primary_source(relParam *_rel_param, double _energy_shift_source_observer) {
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


////////////////////////////////////////////
class PrimarySource {
  ////////////////////////////////////////////

 public:
  PrimarySource(const ModelParams &_model_params, RelSysPar *sys_par, const XspecSpectrum &_xspec_spec) :
      source_parameters{_model_params},
      m_lp_refl_frac((sys_par == nullptr) ? nullptr : sys_par->emis->photon_fate_fractions),
      m_prime_spec_observer{get_observed_primary_spectrum(_xspec_spec, source_parameters)} {
  }

  /**
   * @brief calculates the observed primary spectrum, normalized according to Dauser+2016
   * in the source frame
   * @details The normalization is always given in the source frame (for the parameters of
   * Ecut/kTe in this frame).
   * @param _xspec_spec
   * @param _parameters
   * @return
   */
  Spectrum get_observed_primary_spectrum(const XspecSpectrum &_xspec_spec,
                                         const PrimarySourceParameters &_parameters) {
    int status = EXIT_SUCCESS;
    auto spec = Spectrum(_xspec_spec.energy, _xspec_spec.num_flux_bins());
    calc_primary_spectrum(spec.flux, spec.energy(), spec.num_flux_bins,
                          _parameters.xilltab_param(), &status,
                          _parameters.energy_shift_source_observer());

    // take the xillver normalization factor into account
    spec.multiply_flux_by(calc_normalization_factor_source());

    return spec;
  }

  double get_boost_parameter() const {
    return source_parameters.get_boost_parameter(m_lp_refl_frac);
  }

  void print_reflection_strength(const XspecSpectrum &refl_spec, const Spectrum &primary_spec) const;
  void add_primary_spectrum(const XspecSpectrum &reflection_spectrum);

  PrimarySourceParameters source_parameters;

 private:

  lpReflFrac *m_lp_refl_frac = nullptr;
  Spectrum m_prime_spec_observer;  // (xillver) normalized spectrum


  /**
   * @brief calculate the normalization factor of the observed primary source spectrum
   * according to the xillver normalization given in Dauser+2016.
   *
   * It is defined such that F_E(E) = calc_primary_spectrum*norm_factor follows the
   * definition of the normalization
   * @details note that for the LP the normalization is applied for the spectrum in the frame
   * of the source (and NOT the observer); the energy shift is automatically taken into account
   * @return norm_factor [double]
   */
  double calc_normalization_factor_source() const {
    int status = EXIT_SUCCESS;

    // need to use a specific energy grid for the primary component normalization to
    // fulfill the XILLVER NORM condition (Dauser+2016)
    EnerGrid *egrid = get_coarse_xillver_energrid(&status);
    auto prim_spec_source = new double[egrid->nbins];
    calc_primary_spectrum(prim_spec_source, egrid->ener, egrid->nbins, source_parameters.xilltab_param(), &status);

    // calculate the normalization at the primary source (i.e., not shifted by the energy shift)
    double const norm_factor_prim_spec =
        1. / calcNormWrtXillverTableSpec(prim_spec_source, egrid->ener, egrid->nbins, &status);
    delete[] prim_spec_source;
    return norm_factor_prim_spec;
  }

  // calculate the flux of the observed primary spectrum in ergs/cm2/sec
  double get_normalized_primary_spectrum_flux_in_ergs() const {

    const double *ener = m_prime_spec_observer.energy();

    double ener_flux = 0.0;
    for (int ii = 0; ii < m_prime_spec_observer.num_flux_bins; ii++) {
      if (ener[ii] >= EMIN_XILLVER_NORMALIZATION && ener[ii + 1] < EMAX_XILLVER_NORMALIZATION) {
        ener_flux += m_prime_spec_observer.flux[ii] * 0.5 * (ener[ii] + ener[ii + 1]);
      }
    }

    return ener_flux * CONVERT_KEV2ERG;
  }

};

#endif //RELXILL_SRC_PRIMARYSOURCE_H_
