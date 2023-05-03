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
#include "Relbase.h"

/*
  * @brief: print the reflection strength on the screen (requires struct_refl_frac to be set)
  * note, the energy grids of both spectra need to be the same
  */
void PrimarySource::print_reflection_strength(const XspecSpectrum &refl_spec, const Spectrum &primary_spec) const {

  if (m_lp_refl_frac != nullptr) {
    return;
  } // will do nothing if the refl_frac structure is not set

  // energy grids need to be the same
  assert(refl_spec.num_flux_bins() == static_cast<int>(primary_spec.num_flux_bins));

  const relParam *rel_param = source_parameters.rel_param();

  int const imin = binary_search(refl_spec.energy, refl_spec.num_flux_bins() + 1, RSTRENGTH_EMIN);
  int const imax = binary_search(refl_spec.energy, refl_spec.num_flux_bins() + 1, RSTRENGTH_EMAX);

  double sum_pl = 0.0;
  double sum = 0.0;
  for (int ii = imin; ii <= imax; ii++) {
    sum_pl += primary_spec.flux[ii];
    sum += refl_spec.flux[ii];
  }

  printf("For a = %.3f, Rin = %.3f, and h = %.2f rg", rel_param->a, rel_param->rin, rel_param->height);
  if (is_iongrad_model(rel_param->ion_grad_type) || rel_param->beta > 1e-6) {
    printf(" and beta=%.3f v/c", rel_param->beta);
  }
  printf(" (using boost=1): \n - reflection fraction  %.3f \n - reflection strength is: %.3f \n",
         m_lp_refl_frac->refl_frac,
         sum / sum_pl);
  printf(" - photons falling into the black hole or plunging region: %.2f%%\n", m_lp_refl_frac->f_bh * 100);
  printf(" - energy shift from the primary source to the observer is %.3f\n", energy_shift_source_obs(rel_param));
}

// @brief: adds primary reflection_spectrum to the input reflection_spectrum
void PrimarySource::add_primary_spectrum(const XspecSpectrum &reflection_spectrum) {

  // get the primary spectrum on the energy grid of the reflection spectrum
  auto primary_spectrum = get_observed_primary_spectrum(reflection_spectrum);

  // For the non-relativistic model and if not the LP geometry, we simply multiply by the reflection fraction
  if (is_xill_model(source_parameters.model_type()) || source_parameters.emis_type() != EMIS_TYPE_LP) {
    reflection_spectrum.multiply_flux_by(fabs(source_parameters.refl_frac()));

  } else { // we are in the LP geometry

    assert(m_lp_refl_frac != nullptr);

    double refl_frac_input = source_parameters.refl_frac();
    // if set, it is given as boost, wrt predicted refl_frac_input
    if (source_parameters.interpret_reflfrac_as_boost())
      refl_frac_input *= m_lp_refl_frac->refl_frac;


    const double prime_spec_factors_source_to_observer =
        source_parameters.flux_boost_source_to_observer((*m_lp_refl_frac));

    // if the user sets the refl_frac_input parameter manually, we need to calculate the ratio
    // to end up with the correct normalization
    double const norm_fac_refl = (fabs(refl_frac_input)) / m_lp_refl_frac->refl_frac;

    if (is_alpha_model(source_parameters.model_type())) {
      //  for the alpha model, which includes the distance, the normalization is defined over "norm_flux_cgs"
      double const
          prime_norm_factor = source_parameters.norm_flux_cgs() / get_normalized_primary_spectrum_flux_in_ergs();

      primary_spectrum.multiply_flux_by(prime_norm_factor);

      // multiply the reflection reflection_spectrum also with the normalization factor and the inverse of primary source
      // factors to take those into account in the reflection reflection_spectrum
      reflection_spectrum.multiply_flux_by(prime_norm_factor * norm_fac_refl / prime_spec_factors_source_to_observer);

    } else { // not the ALPHA model
      reflection_spectrum.multiply_flux_by(norm_fac_refl);
      primary_spectrum.multiply_flux_by(prime_spec_factors_source_to_observer);
    }

    // if desired, we output the reflection fraction and strength (as defined in Dauser+2016)
    if (shouldAuxInfoGetPrinted()) {
      print_reflection_strength(reflection_spectrum, primary_spectrum);
    }
  }

  // Finally, add the power law component if refl_frac >= 0
  if (source_parameters.refl_frac() >= 0) {
    for (int ii = 0; ii < reflection_spectrum.num_flux_bins(); ii++) {
      reflection_spectrum.flux[ii] += primary_spectrum.flux[ii];
    }
  }

}
