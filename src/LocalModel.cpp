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

#include "LocalModel.h"
#include "XspecSpectrum.h"

#include <stdexcept>
#include <iostream>



/*
 * @brief calculate line model
 */
void LocalModel::line_model(const XspecSpectrum &spectrum) {

  relParam *rel_param = LocalModel::get_rel_params();

  // relline_base calculates the line for 1keV -> shift the energy grid accordingly
  spectrum.shift_energy_grid_1keV(rel_param->lineE);

  int status = EXIT_SUCCESS;

  relline_base(spectrum.energy, spectrum.flux, spectrum.num_flux_bins(), rel_param, &status);

  if (status != EXIT_SUCCESS) {
    throw std::exception();
  }

  delete rel_param;
}

/*
 * @brief calculate relxill model
 */
void LocalModel::relxill_model(const XspecSpectrum &spectrum) {


  int status = EXIT_SUCCESS;
  if (m_model_params.irradiation() == T_Irrad::Const){
    relParam *rel_param = LocalModel::get_rel_params();
    xillParam *xill_param = LocalModel::get_xill_params();
    relxill_bb_kernel(spectrum.energy, spectrum.flux, spectrum.num_flux_bins(), xill_param, rel_param, &status);
    delete rel_param;
    delete xill_param;
  } else {
    relxill_kernel(spectrum, m_model_params, &status);
  }

  if (status != EXIT_SUCCESS) {
    throw std::exception();
  }

}

/**
 * @brief calculate convolution of a given spectrum
 */
void LocalModel::conv_model(const XspecSpectrum &spectrum) {

  if (calcSum(spectrum.flux, spectrum.num_flux_bins()) <= 0.0) {
    throw ModelEvalFailed("input flux for convolution model needs to be >0");
  }

  relParam *rel_param = LocalModel::get_rel_params();

  int status = EXIT_SUCCESS;
  relconv_kernel(spectrum.energy, spectrum.flux, spectrum.num_flux_bins(), rel_param, &status);

  if (status != EXIT_SUCCESS) {
    throw ModelEvalFailed("executing convolution failed");
  }

  delete rel_param;
}


/**
 * @brief calculate xillver model
 */
void LocalModel::xillver_model(const XspecSpectrum &spectrum) {

  xillParam *xill_param = LocalModel::get_xill_params();

  int status = EXIT_SUCCESS;
  xillver_base(spectrum.energy, spectrum.num_flux_bins(), spectrum.flux, xill_param, &status);

  if (status != EXIT_SUCCESS) {
    throw std::exception();
  }

  delete xill_param;
}


/**
 * Wrapper function, which can be called from any C-type local model function as
 * required by Xspec
 * @param model_name: unique name of the model
 * @param parameter_values: input parameter_values array (size can be determined from the model definition for the model_name)
 * @param xspec_energy[num_flux_bins]: input energy grid
 * @param xspec_flux[num_flux_bins]: - flux array (already allocated), used to return the calculated values
 *                    - for convolution models this is also the input flux
 * @param num_flux_bins
 */
void xspec_C_wrapper_eval_model(ModelName model_name,
                                const double *parameter_values,
                                double *xspec_flux,
                                int num_flux_bins,
                                const double *xspec_energy) {

  try {
    LocalModel local_model{parameter_values, model_name};

    XspecSpectrum spectrum{xspec_energy, xspec_flux, static_cast<size_t>(num_flux_bins)};
    local_model.eval_model(spectrum);

  } catch (ModelNotFound &e) {
    std::cout << e.what();
  }
    // TODO: what should we do if the evaluation fails? return zeros?

}
