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

void LocalModel::line_model(CppSpectrum &spectrum) {

  puts(" I am LINE ");

  // get the relevant parameters
  relParam *rel_param = getRelParamStruct(m_param, m_name, m_info);

  // shift the spectrum such that we can calculate the line for 1 keV
  spectrum.shift_energy_grid_1keV(rel_param->lineE, rel_param->z);

  int status = EXIT_SUCCESS;

  auto flux_double = spectrum.flux_double();

  relline_base(spectrum.energy_double(), flux_double, spectrum.nener_bins(), rel_param, &status);

  if (status != EXIT_SUCCESS) {
    throw std::exception();
  }

  // need to (internally) copy the values from the double-array to the Array structures
  spectrum.copy_doubleArrays2array();

  delete rel_param;
}

void relxill_model(const CppSpectrum &spectrum) {
  //  relxill_kernel(ener, flux, n_ener0, xill_param, rel_param, status);
  //  double *ener_shifted = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);
  //  rebin_spectrum(ener_shifted, photar, n_ener0, ener, flux, n_ener0);
}

void conv_model(const CppSpectrum &spectrum) {
  // shift the spectrum such that we can calculate the line for 1 keV
  //  double *ener1keV = shift_energ_spec_1keV(ener, n_ener, param_struct->lineE, param_struct->z, status);
  //  CHECK_STATUS_VOID(*status);
  //
  //  relconv_kernel(ener1keV, photar, n_ener, rel_param, status);
}

void xillver_model(const CppSpectrum &spectrum) {
  //  xillver_base(ener0, n_ener0, photar, xill_param, status);
}

void xspec_wrapper_eval_model(ModelName model_name, const Array &energy, Array &flux, const Array &parameter) {

  try {
    auto const model_definition = ModelDatabase::instance().get(model_name);

    CppSpectrum spectrum{energy, flux};

    ModelParams params{model_definition.input_parameters(), parameter};

    LocalModel model{params, model_name};
    model.eval_model(spectrum);

    flux = spectrum.flux();

  } catch (std::out_of_range &e) {
    std::cout << " *** relxill-error: required model not found in database " << std::endl;
    throw e;
  } catch (std::bad_alloc &e) {
    std::cout << " *** relxill-error: allocation error " << std::endl;
    throw e;
  }

}

extern "C" {

void lmodcpprelline(const Array &energy, const Array &parameter,
                                     int spectrum, Array &flux, Array &fluxError,
                                     const string &init) {
  xspec_wrapper_eval_model(ModelName::relline, energy, flux, parameter);
}

void lmodcpprelxill(const Array &energy, const Array &parameter,
                                     int spectrum, Array &flux, Array &fluxError,
                                     const string &init) {
  xspec_wrapper_eval_model(ModelName::relxill, energy, flux, parameter);
}

}

