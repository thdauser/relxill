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

#include <stdexcept>
#include <iostream>

// namespace  relxill{

void eval_model_xspec(ModelName model, const Array &energy, const Array &flux, const Array &parameter) {

  try {
    auto const model_info = LocalModels::instance().getModelInfo(model);

    //lmod{model, parameter};

    // Spectrum spec{energy, flux};
    //  Parameters pars{parameter, std::move(par_names)};


    switch (model_info.type()) {
      case T_Model::LineModel:puts("I am LINE ");
        break;
      case T_Model::RelxillModel:puts(" I am RELXILL ");
        break;
      case T_Model::ConvModel:puts(" I am CONV ");
        break;
      case T_Model::XillModel:puts(" I am XILL ");
        break;
    }

  } catch (std::out_of_range &e) {
    std::cout << " *** relxill-error: required model not found in database " << std::endl;
    throw e;
  }

}

extern "C" {

[[maybe_unused]] void lmodcpprelline(const Array &energy, const Array &parameter,
                                     int spectrum, Array &flux, Array &fluxError,
                                     const string &init) {
  eval_model_xspec(ModelName::relxilllp, energy, flux, parameter);
}

[[maybe_unused]] void lmodcpprelxill(const Array &energy, const Array &parameter,
                                     int spectrum, Array &flux, Array &fluxError,
                                     const string &init) {
  eval_model_xspec(ModelName::relxill, energy, flux, parameter);
}

}

// }   // namespace relxill

