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



class PrimarySourceParameters {

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
class PrimarySource {

 public:
  PrimarySource(const ModelParams &_model_params, RelSysPar *sys_par, const XspecSpectrum &_xspec_spec) :
      parameters{_model_params},
      m_struct_refl_frac((sys_par == nullptr) ? nullptr : sys_par->emis->photon_fate_fractions),
      m_prime_spec_observer{calculate_observed_primary_spectrum(_xspec_spec, parameters)} {

    m_prime_spec_xillver_norm_factor = 1. / calculate_source_prime_spec_xillver_norm();
    // source to obs normalization
    //
  }

  /*
   * @brief: print the reflection strength on the screen (requires struct_refl_frac to be set)
   */
  void print_reflection_strength(const XspecSpectrum &spectrum, const double *pl_flux) const {

    if (m_struct_refl_frac != nullptr) {
      return;
    } // will do nothing if the refl_frac structure is not set

    const relParam *rel_param = parameters.rel_param();

    int const imin = binary_search(spectrum.energy, spectrum.num_flux_bins() + 1, RSTRENGTH_EMIN);
    int const imax = binary_search(spectrum.energy, spectrum.num_flux_bins() + 1, RSTRENGTH_EMAX);

    double sum_pl = 0.0;
    double sum = 0.0;
    for (int ii = imin; ii <= imax; ii++) {
      sum_pl += pl_flux[ii];
      sum += spectrum.flux[ii];
    }

    printf("For a = %.3f, Rin = %.3f, and h = %.2f rg", rel_param->a, rel_param->rin, rel_param->height);
    if (is_iongrad_model(rel_param->ion_grad_type) || rel_param->beta > 1e-6) {
      printf(" and beta=%.3f v/c", rel_param->beta);
    }
    printf(" (using boost=1): \n - reflection fraction  %.3f \n - reflection strength is: %.3f \n",
           m_struct_refl_frac->refl_frac,
           sum / sum_pl);
    printf(" - photons falling into the black hole or plunging region: %.2f%%\n", m_struct_refl_frac->f_bh * 100);
    printf(" - energy shift from the primary source to the observer is %.3f\n", energy_shift_source_obs(rel_param));
  }

  double calculate_source_prime_spec_xillver_norm() {
    int status = EXIT_SUCCESS;
    // need to create a specific energy grid for the primary component normalization to
    // fulfill the XILLVER NORM condition (Dauser+2016)
    EnerGrid *egrid = get_coarse_xillver_energrid(&status);
    auto prim_spec_source = new double[egrid->nbins];
    calc_primary_spectrum(prim_spec_source, egrid->ener, egrid->nbins,
                          parameters.xilltab_param(), &status);

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

  // @brief: adds primary spectrum to the input spectrum
  void add_primary_spectrum(const XspecSpectrum &spectrum) {

    int status = EXIT_SUCCESS;
    double *pl_flux = calc_normalized_primary_spectrum(spectrum.energy, spectrum.num_flux_bins(),
                                                       parameters.rel_param(), parameters.xilltab_param(), &status);
    CHECK_STATUS_VOID(status);

    // For the non-relativistic model and if not the LP geometry, we simply multiply by the reflection fraction
    if (is_xill_model(parameters.model_type()) || parameters.emis_type() != EMIS_TYPE_LP) {
      spectrum.multiply_flux_by(fabs(parameters.refl_frac()));

    } else { // we are in the LP geometry

      assert(m_struct_refl_frac != nullptr);

      double refl_frac = parameters.refl_frac();  // todo: disentangle refl_frac and boost
      if (parameters.interpret_reflfrac_as_boost()) {
        // if set, it is given as boost, wrt predicted refl_frac
        refl_frac *= m_struct_refl_frac->refl_frac;
      }

      double prime_spec_factors_source_to_observer =
          m_struct_refl_frac->f_inf_rest / 0.5 * pow(parameters.energy_shift_source_observer(), parameters.gam());
      // flux boost of primary radiation taking into account here (therfore we need f_inf_rest above)
      if (parameters.beta() > 1e-4) {
        prime_spec_factors_source_to_observer *= pow(doppler_factor_source_obs(parameters.rel_param()), 2);
      }

      // if the user sets the refl_frac parameter manually, we need to calculate the ratio
      // to end up with the correct normalization
      double const norm_fac_refl = (fabs(refl_frac)) / m_struct_refl_frac->refl_frac;

      if (is_alpha_model(parameters.model_type())) {
        //  for the alpha model, which includes the distance, the normalization is defined over "normalization_factor"
        double const
            prime_norm_factor = parameters.normalization_factor() / get_normalized_primary_spectrum_flux_in_ergs();

        for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) {
          pl_flux[ii] *= prime_norm_factor;
        }
        // multiply the reflection spectrum also with the normalization factor and the inverse of primary source
        // factors to take those into account in the reflection spectrum
        spectrum.multiply_flux_by(prime_norm_factor * norm_fac_refl / prime_spec_factors_source_to_observer);
      } else {
        spectrum.multiply_flux_by(norm_fac_refl);
        for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) {
          pl_flux[ii] *= prime_spec_factors_source_to_observer;
        }
      }

      /** 5 ** if desired, we ouput the reflection fraction and strength (as defined in Dauser+2016) **/
      if (shouldAuxInfoGetPrinted()) {
        print_reflection_strength(spectrum, pl_flux);
      }

    }

    // Finally, add the power law component if refl_frac >= 0
    if (parameters.refl_frac() >= 0) {
      for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) {
        spectrum.flux[ii] += pl_flux[ii];
      }
    }

    delete[] pl_flux;

  }

  PrimarySourceParameters parameters;

 private:

  lpReflFrac *m_struct_refl_frac = nullptr;
  Spectrum m_prime_spec_observer;
  double m_prime_spec_xillver_norm_factor = 1.0;

  // calculate the flux of the primary spectrum in ergs/cm2/sec
  double get_normalized_primary_spectrum_flux_in_ergs() const {
    int status = EXIT_SUCCESS;
    EnerGrid *egrid = get_coarse_xillver_energrid(&status);
    auto prim_spec_source = new double[egrid->nbins];
    double *pl_flux = calc_normalized_primary_spectrum(egrid->ener, egrid->nbins,
                                                       parameters.rel_param(), parameters.xilltab_param(), &status);

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
