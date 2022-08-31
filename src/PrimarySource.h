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

#include "common.h"
#include "ModelParams.h"
#include "Relphysics.h"

#include <vector>

typedef struct {
} primeSourceParam;

class PrimarySourceParameters {

 public:

  explicit PrimarySourceParameters(const ModelParams &inp_param) {

    rel_param = get_rel_params(inp_param);  // can be nullptr if xillver model
    xill_param = get_xill_params(inp_param);

    prim_type = xill_param->prim_type;
    emis_type = rel_param->emis_type;

    ect = xill_param->ect;  // can be Ecut or kTe depending on the xillver model
    kTbb = xill_param->kTbb;

    // energy shiftset to "1" for non-LP and xillver models
    energy_shift_source_observer = (emis_type == EMIS_TYPE_LP) ?
                                   energy_shift_source_obs(rel_param) :
                                   1.0;

    // special case, for the LP model and the Ecut model, the cutoff energy is given in the observer frame
    // -> convert it such that ecut is also given in the source frame
    if (emis_type == EMIS_TYPE_LP && prim_type == PRIM_SPEC_ECUT) {
      ect /= energy_shift_source_observer;
    }

    // those values should never be used, unless it is set by the model
    refl_frac = inp_param.get_otherwise_default(XPar::refl_frac, 0);
    mass = inp_param.get_otherwise_default(XPar::mass, 0);
    distance = inp_param.get_otherwise_default(XPar::distance, 0);
    lbol = inp_param.get_otherwise_default(XPar::lbol, 0);

    interpret_reflfrac_as_boost =
        static_cast<int>(lround(inp_param.get_otherwise_default(XPar::switch_switch_reflfrac_boost, 0)));

  }

  ~PrimarySourceParameters() {
    delete rel_param;
    delete xill_param;
  }

  double lbol;
  double mass;
  double distance;
  int prim_type;
  int emis_type;
  double refl_frac;
  int interpret_reflfrac_as_boost;
  double ect;
  double kTbb;
  double energy_shift_source_observer;

  relParam *rel_param = nullptr;
  xillParam *xill_param = nullptr;

};

// TODO:
//  - properly include the reflection fraction and a nice error message
class PrimarySource {

 public:
  explicit PrimarySource(const ModelParams &_model_params) :
      m_param{PrimarySourceParameters(_model_params)} {
  }

  PrimarySource(const ModelParams &_model_params, lpReflFrac *struct_refl_frac) :
      m_param{_model_params}, m_struct_refl_frac{struct_refl_frac} {
  }

  void print_reflection_strength(const XspecSpectrum &spectrum,
                                 const double *pl_flux) const {

    assert(m_struct_refl_frac != nullptr);

    relParam *rel_param = m_param.rel_param;

    // todo: all this to be set by a ENV
    int imin = binary_search(spectrum.energy, spectrum.num_flux_bins() + 1, RSTRENGTH_EMIN);
    int imax = binary_search(spectrum.energy, spectrum.num_flux_bins() + 1, RSTRENGTH_EMAX);

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

  // @brief: adds primary spectrum to the input spectrum
  void add_primary_spectrum(const XspecSpectrum &spectrum) {

    int status = EXIT_SUCCESS;

    xillTableParam *xill_table_param =
        get_xilltab_param(m_param.xill_param, &status);  // make this separate for the primary source

    double *pl_flux = calc_normalized_primary_spectrum(spectrum.energy, spectrum.num_flux_bins(),
                                                       m_param.rel_param, xill_table_param, &status);
    free(xill_table_param);
    CHECK_STATUS_VOID(status);

    // For the non-relativistic model and if not the LP geometry, we simply multiply by the reflection fraction
    if (is_xill_model(m_param.rel_param->model_type) || m_param.emis_type != EMIS_TYPE_LP) {
      spectrum.multiply_flux_by(fabs(m_param.refl_frac));

    } else { // we are in the LP geometry

      assert(m_struct_refl_frac != nullptr);

      if (m_param.interpret_reflfrac_as_boost) {
        // if set, it is given as boost, wrt predicted refl_frac
        m_param.refl_frac *= m_struct_refl_frac->refl_frac;
      }

      double prim_fac =
          m_struct_refl_frac->f_inf_rest / 0.5 * pow(m_param.energy_shift_source_observer, m_param.xill_param->gam);
      // flux boost of primary radiation taking into account here (therfore we need f_inf_rest above)
      if (m_param.rel_param->beta > 1e-4) {
        prim_fac *= pow(doppler_factor_source_obs(m_param.rel_param), 2);
      }

      // if the user sets the refl_frac parameter manually, we need to calculate the ratio
      // to end up with the correct normalization
      double norm_fac_refl = (fabs(m_param.xill_param->refl_frac)) / m_struct_refl_frac->refl_frac;

      spectrum.multiply_flux_by(norm_fac_refl);
      for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) {
        pl_flux[ii] *= prim_fac;
      }

      /** 5 ** if desired, we ouput the reflection fraction and strength (as defined in Dauser+2016) **/
      if (shouldAuxInfoGetPrinted()) {
        print_reflection_strength(spectrum, pl_flux);
      }

    }

    // Finally, add the power law component if refl_frac >= 0
    if (m_param.refl_frac >= 0) {
      for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) {
        spectrum.flux[ii] += pl_flux[ii];
      }
    }

    delete[] pl_flux;

  }

 private:
  PrimarySourceParameters m_param;
  lpReflFrac *m_struct_refl_frac = nullptr;
};

#endif //RELXILL_SRC_PRIMARYSOURCE_H_
