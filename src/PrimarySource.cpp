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
#include "PrimarySource.h"

/*
  * @brief: print the reflection strength on the screen (requires struct_refl_frac to be set)
  */
void PrimarySource::print_reflection_strength(const XspecSpectrum &spectrum, const double *pl_flux) const {

  if (m_struct_refl_frac != nullptr) {
    return;
  } // will do nothing if the refl_frac structure is not set

  const relParam *rel_param = source_parameters.rel_param();

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

// @brief: adds primary spectrum to the input spectrum
void PrimarySource::add_primary_spectrum(const XspecSpectrum &spectrum) {

  int status = EXIT_SUCCESS;
  double *pl_flux = calc_normalized_primary_spectrum(spectrum.energy,
                                                     spectrum.num_flux_bins(),
                                                     source_parameters.rel_param(),
                                                     source_parameters.xilltab_param(),
                                                     &status);
  CHECK_STATUS_VOID(status);

  // For the non-relativistic model and if not the LP geometry, we simply multiply by the reflection fraction
  if (is_xill_model(source_parameters.model_type()) || source_parameters.emis_type() != EMIS_TYPE_LP) {
    spectrum.multiply_flux_by(fabs(source_parameters.refl_frac()));

  } else { // we are in the LP geometry

    assert(m_struct_refl_frac != nullptr);

    double refl_frac = source_parameters.refl_frac();  // todo: disentangle refl_frac and boost
    if (source_parameters.interpret_reflfrac_as_boost()) {
      // if set, it is given as boost, wrt predicted refl_frac
      refl_frac *= m_struct_refl_frac->refl_frac;
    }

    double prime_spec_factors_source_to_observer =
        m_struct_refl_frac->f_inf_rest / 0.5
            * pow(source_parameters.energy_shift_source_observer(), source_parameters.gam());
    // flux boost of primary radiation taking into account here (therefore we need f_inf_rest above)
    if (source_parameters.beta() > 1e-4) {
      prime_spec_factors_source_to_observer *= pow(doppler_factor_source_obs(source_parameters.rel_param()), 2);
    }

    // if the user sets the refl_frac parameter manually, we need to calculate the ratio
    // to end up with the correct normalization
    double const norm_fac_refl = (fabs(refl_frac)) / m_struct_refl_frac->refl_frac;

    if (is_alpha_model(source_parameters.model_type())) {
      //  for the alpha model, which includes the distance, the normalization is defined over "normalization_factor"
      double const
          prime_norm_factor = source_parameters.normalization_factor() / get_normalized_primary_spectrum_flux_in_ergs();

      for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) {
        pl_flux[ii] *= prime_norm_factor;
      }
      // multiply the reflection spectrum also with the normalization factor and the inverse of primary source
      // factors to take those into account in the reflection spectrum
      spectrum.multiply_flux_by(prime_norm_factor * norm_fac_refl / prime_spec_factors_source_to_observer);

    } else { // not the ALPHA model
      spectrum.multiply_flux_by(norm_fac_refl);
      for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) {
        pl_flux[ii] *= prime_spec_factors_source_to_observer;
      }
    }

    /** 5 ** if desired, we output the reflection fraction and strength (as defined in Dauser+2016) **/
    if (shouldAuxInfoGetPrinted()) {
      print_reflection_strength(spectrum, pl_flux);
    }

  }

  // Finally, add the power law component if refl_frac >= 0
  if (source_parameters.refl_frac() >= 0) {
    for (int ii = 0; ii < spectrum.num_flux_bins(); ii++) {
      spectrum.flux[ii] += pl_flux[ii];
    }
  }

  delete[] pl_flux;

}
