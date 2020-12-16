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

    Copyright 2020 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "cppmodels.h"
#include "cppspectrum.h"

#include <stdexcept>
#include <iostream>

void LocalModel::line_model(const XspecSpectrum &spectrum) {

  puts(" I am LINE ");

  // get the relevant parameters
  relParam *rel_param = getRelParamStruct(m_param, m_name, m_info);

  // shift the spectrum such that we can calculate the line for 1 keV
  spectrum.shift_energy_grid_1keV(rel_param->lineE, rel_param->z);

  int status = EXIT_SUCCESS;

  relline_base(spectrum.energy(), spectrum.flux(), spectrum.num_flux_bins(), rel_param, &status);

  if (status != EXIT_SUCCESS) {
    throw std::exception();
  }

  delete rel_param;
}

void LocalModel::relxill_model(const XspecSpectrum &spectrum) {

  // get the relevant parameters
  // relParam *rel_param = getRelParamStruct(m_param, m_name, m_info);
  // relParam *xill_param = getXillParamStruct(m_param, m_name, m_info);


  //int status = EXIT_SUCCESS;
  //relxill_kernel(spectrum.energy(), spectrum.flux(), spectrum.num_flux_bins(), xill_param, rel_param, &status);
  //  double *ener_shifted = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);
  //  rebin_spectrum(ener_shifted, photar, n_ener0, ener, flux, n_ener0);
}

void LocalModel::conv_model(const XspecSpectrum &spectrum) {
  // shift the spectrum such that we can calculate the line for 1 keV
  //  double *ener1keV = shift_energ_spec_1keV(ener, n_ener, param_struct->lineE, param_struct->z, status);
  //  CHECK_STATUS_VOID(*status);
  //
  //  relconv_kernel(ener1keV, photar, n_ener, rel_param, status);
}

void LocalModel::xillver_model(const XspecSpectrum &spectrum) {
  //  xillver_base(ener0, n_ener0, photar, xill_param, status);
}

/**
 * Wrapper function, which can be called from any C-type local model function as
 * required by Xspec
 * @param model_name: unique name of the model
 * @param parameter: input parameter array (size can be determined from the model definition for the model_name)
 * @param xspec_energy[num_flux_bins]: input energy grid
 * @param xspec_flux[num_flux_bins]: - flux array (already allocated), used to return the calculated values
 *                    - for convolution models this is also the input flux
 * @param num_flux_bins
 */
void xspec_C_wrapper_eval_model(ModelName model_name,
                                const double *parameter,
                                double *xspec_flux,
                                int num_flux_bins,
                                const double *xspec_energy) {

  try {
    auto const model_definition = ModelDatabase::instance().get(model_name);
    ModelParams params{model_definition.input_parameters(), parameter};

    LocalModel model{params, model_name};

    XspecSpectrum spectrum{xspec_energy, xspec_flux, static_cast<size_t>(num_flux_bins)};
    model.eval_model(spectrum);

  } catch (std::exception &e) {
    std::cout << " *** relxill-error: model evaluation failed " << std::endl;
    throw e;
  }

}
