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

#include "ModelDefinition.h"
#include "Relphysics.h"
#include "Xillspec.h"
#include "Relcache.h"
//#include "Rellp_Extended.h" // to calculate energy shift from source to observer (cannot load??)

extern "C" {
#include "common.h"
#include "xilltable.h"
}

#include <vector>



////////////////////////////////////////////
class PrimarySourceParameters {
  ////////////////////////////////////////////

 public:

  explicit PrimarySourceParameters(const ModelDefinition &inp_param) :
      m_inp_param{inp_param},
      m_rel_param{get_rel_params(m_inp_param)},
      m_interpret_reflfrac_as_boost{
          static_cast<int>(lround(inp_param.get_otherwise_default(XPar::switch_switch_reflfrac_boost, 0)))},
      m_energy_shift_source_observer{(m_rel_param != nullptr && m_rel_param->emis_type == EMIS_TYPE_LP) ?
                                     energy_shift_source_obs(m_rel_param) : 1.0},
      m_xilltab_param{m_get_xilltab_params_primary_source(m_rel_param, m_energy_shift_source_observer)} {
    //m_energy_shift_source_observer = (m_rel_param != nullptr && m_rel_param->emis_type == EMIS_TYPE_LP) ?
    //                                 energy_shift_source_obs(m_rel_param) : 1.0;

    auto xill_param = get_xill_params(m_inp_param);
    m_refl_frac = xill_param->refl_frac;
    m_distance_kpc = xill_param->distance;
    m_mass_msolar = xill_param->mass_msolar;
    m_luminosity_primary_source = xill_param->luminosity_primary_source;
    delete xill_param;

    // m_xilltab_param = m_get_xilltab_params_primary_source(m_rel_param, m_energy_shift_source_observer);
  }

  ~PrimarySourceParameters() {
    delete m_rel_param;
    free(m_xilltab_param);
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
    double prime_spec_factors_source_to_observer;
    if (m_rel_param->model_type == MOD_TYPE_RELXILLLPEXT) {
      // general lensing_factor and energy shift to observer (now calculated in Rellp_Extended and has form lensing_factor * g_{so}^Gamma)
      prime_spec_factors_source_to_observer = lp_refl_frac.lensing_and_boost_factor;
  // For some reason I cannot import Rellp_Ext.h in here. Therefore cannot use this function and need to put boost into lensing:
      //prime_spec_factors_source_to_observer *= get_energyshift_gamma_ext_source_obs(m_rel_param);
      // todo: ideally, m_energy_shift_source_observer^Gamma should be calculated from extended table for ext models
      // then lp_refl_frac.lensing_factor will be indeed lensing_factor

      assert(prime_spec_factors_source_to_observer > 0.0);
    } else {
      // as it is in LP RELXILL
      prime_spec_factors_source_to_observer =
          lp_refl_frac.f_inf_rest / 0.5 * pow(m_energy_shift_source_observer, PrimarySourceParameters::gam());
      /// todo after Gullo and Adam discussion: for nthcomp the situation is different, boost is not g^gamma but just g, check that it is correct and update of Cp models required
    }

    // flux boost of primary radiation taking into account here (therefore we need f_inf_rest above)
    if (PrimarySourceParameters::beta() > 1e-4) {
      // todo need to approve/modify for extended source
      prime_spec_factors_source_to_observer *= pow(doppler_factor_source_obs(m_rel_param), 2);
    }
    return prime_spec_factors_source_to_observer;
  }

  /**
   * @brief: return the intrinsic source luminosity
   * @return luminosity source [ergs/s]
   */
  double luminosity_source_cgs() const {
    return m_luminosity_primary_source;
  }

  /**
 * @brief: calculate the flux boost from the source spectrum to the observed spectrum
 * @param f_inf
 * @return luminosity source [ergs/s]
 */
  double flux_observed_cgs(const lpReflFrac &lp_refl_frac) const {

    const double CONST_cm2kpc = 3.2407792700054E-22;

    const double distance_cm = m_distance_kpc / CONST_cm2kpc;
    const double flux_observed_cgs = m_luminosity_primary_source / (4 * M_PI * distance_cm * distance_cm)
        * flux_boost_source_to_observer(lp_refl_frac);

    // add energy shift and flux boost
    return flux_observed_cgs;
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
  double get_boost_parameter(const lpReflFrac *lp_refl_frac) const {

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

  ModelDefinition m_inp_param;
  relParam *m_rel_param = nullptr;
  int m_interpret_reflfrac_as_boost;

  double m_energy_shift_source_observer;

  double m_luminosity_primary_source;
  double m_distance_kpc;
  double m_refl_frac;
  double m_mass_msolar;

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
  PrimarySource(const ModelDefinition &_model_params, RelSysPar *sys_par) :
      source_parameters{_model_params},
      m_lp_refl_frac{sys_par->emis->photon_fate_fractions}
      {  }

  /**
   * @brief calculates the observed primary spectrum, normalized according to Dauser+2016
   * in the source frame
   * @details The normalization is always given in the source frame (for the parameters of
   * Ecut/kTe in this frame).
   * @param spec
   * @param _parameters
   * @return
   */
  Spectrum get_observed_primary_spectrum(const Spectrum &relxill_spec) const {
    int status = EXIT_SUCCESS;
    auto prim_spec = Spectrum(relxill_spec.energy(), relxill_spec.num_flux_bins);
    calc_primary_spectrum(prim_spec.flux, prim_spec.energy(), prim_spec.num_flux_bins,
                          source_parameters.xilltab_param(), &status,
                          source_parameters.energy_shift_source_observer());

    // take the xillver normalization factor into account
    prim_spec.multiply_flux_by(calc_normalization_factor_source());

    return prim_spec;
  }

  double get_boost_parameter() const {
    return source_parameters.get_boost_parameter(m_lp_refl_frac);
  }

  void print_reflection_strength(const Spectrum &refl_spec, const Spectrum &primary_spec) const;

  void add_primary_spectrum(const RelxillSpec &relxill_spec);

  PrimarySourceParameters source_parameters;

 private:

  lpReflFrac *m_lp_refl_frac = nullptr;
  double m_emin_norm = get_env_otherwise_default("RELXILL_ALPHA_CGSNORM_EMIN", EMIN_XILLVER_NORMALIZATION);
  double m_emax_norm = get_env_otherwise_default("RELXILL_ALPHA_CGSNORM_EMAX", EMAX_XILLVER_NORMALIZATION);

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
    EnerGrid *egrid = get_coarse_xillver_energrid();
    auto prim_spec_source = new double[egrid->nbins];
    calc_primary_spectrum(prim_spec_source, egrid->ener, egrid->nbins, source_parameters.xilltab_param(), &status);

    // calculate the normalization at the primary source (i.e., not shifted by the energy shift)
    double const norm_factor_prim_spec =
        1. / calcXillverNormFromPrimarySpectrum(prim_spec_source, egrid->ener, egrid->nbins, &status);
    delete[] prim_spec_source;
    return norm_factor_prim_spec;
  }

  /**
   * @brief calculate the flux of the observed primary spectrum in ergs/cm2/sec
   *
   * @details the energy band is given by the Xillver Normalization Definition (from 0.1-1000 keV),
   * or can be set by using the env variables RELXILL_ALPHA_CGSNORM_EMIN and
   * RELXILL_ALPHA_CGSNORM_EMIN
   */
  double get_normalized_primary_spectrum_flux_in_ergs() const {

    EnerGrid *egrid = get_coarse_xillver_energrid();
    const auto spec = Spectrum(egrid->ener, egrid->nbins); // egrid has num_flux_bins+1 energy bins
    //  const auto xspec_spec = spec.get_xspec_spectrum();

    const auto prime_spec = PrimarySource::get_observed_primary_spectrum(spec);
    const auto ener = prime_spec.energy();

    //    const double emin_norm = get_env_otherwise_default("RELXILL_ALPHA_CGSNORM_EMIN", EMIN_XILLVER_NORMALIZATION);
    //    const double emax_norm = get_env_otherwise_default("RELXILL_ALPHA_CGSNORM_EMIN", EMAX_XILLVER_NORMALIZATION);

    // need this for the whole energy range
    double ener_flux = 0.0;
    for (size_t ii = 0; ii < prime_spec.num_flux_bins; ii++) {
      if (ener[ii] >= m_emin_norm && ener[ii + 1] < m_emax_norm) {
        ener_flux += prime_spec.flux[ii] * 0.5 * (ener[ii] + ener[ii + 1]);
      }
    }

    return ener_flux * CONVERT_KEV2ERG;
  }

};

#endif //RELXILL_SRC_PRIMARYSOURCE_H_
