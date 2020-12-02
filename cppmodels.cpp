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
  }

}

extern "C" {

[[maybe_unused]] void lmodcpprelline(const Array &energy, const Array &parameter,
                                     int spectrum, Array &flux, Array &fluxError,
                                     const string &init) {
  xspec_wrapper_eval_model(ModelName::relline, energy, flux, parameter);
}

[[maybe_unused]] void lmodcpprelxill(const Array &energy, const Array &parameter,
                                     int spectrum, Array &flux, Array &fluxError,
                                     const string &init) {
  xspec_wrapper_eval_model(ModelName::relxill, energy, flux, parameter);
}

}

