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

#ifndef RELXILL__CPPMODELS_H_
#define RELXILL__CPPMODELS_H_

#include <valarray>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// #include <xsTypes.h>

extern "C" {
#include "common.h"
#include "relmodels.h"
}

#include "cppModelDatabase.h"
#include "cppparameters.h"

// namespace relxill {



void line_model(ModelParams param, CppSpectrum spectrum);
void relxill_model(const CppSpectrum &spectrum);
void conv_model(const CppSpectrum &spectrum);
void xillver_model(const CppSpectrum &spectrum);

class LocalModel {

 public:
  LocalModel(const ModelParams &par, ModelName model_name)
      : m_param{par}, m_name{model_name}, m_info{ModelDatabase::instance().get(model_name).model_info()} {
  };

  explicit LocalModel(ModelName model_name) :
      LocalModel(ModelParams(), model_name) {
  };

  void line_model(CppSpectrum &spectrum);

  /**
   * Evaluate the LocalModel (as defined in the class) and overwrite the "flux" array
   * of Spectrum with the output values
   * @param spectrum
   */
  void eval_model(CppSpectrum &spectrum) {

    try {
      switch (m_info.type()) {
        case T_Model::Line:line_model(spectrum);
          break;

        case T_Model::Relxill:puts(" I am RELXILL ");
          relxill_model(spectrum);
          break;

        case T_Model::Conv:puts(" I am CONV ");
          conv_model(spectrum);
          break;

        case T_Model::Xill:puts(" I am XILL ");
          xillver_model(spectrum);
          break;
      }

    } catch (std::exception &e) {
      puts(" *** relxill-error : evaluating model failed ");
      throw e;
    }

  }

 private:
  ModelParams m_param;
  ModelInfo m_info;
  ModelName m_name;
};


void xspec_wrapper_eval_model(ModelName model, const Array &energy, const Array &flux, const Array &parameter);

extern "C" {
[[maybe_unused]] void lmodcpprelline(const Array &energy, const Array &parameter,
                                     int spectrum, Array &flux, Array &fluxError,
                                     const string &init);

[[maybe_unused]] void lmodcpprelxill(const Array &energy, const Array &parameter,
                                     int spectrum, Array &flux, Array &fluxError,
                                     const string &init);
}

// } // namespace relxill


#endif //RELXILL__CPPMODELS_H_
