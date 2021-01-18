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

    Copyright 2021 Thomas Dauser, Remeis Observatory & ECAP
*/

#ifndef RELXILL__CPPMODELS_H_
#define RELXILL__CPPMODELS_H_

#include <valarray>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C" {
#include "common.h"
#include "relmodels.h"
}

#include "cppModelDatabase.h"
#include "cppparameters.h"

/**
 * exception if the model evaluation failed
 */
class ModelEvalFailed : public std::exception {
 public:
  ModelEvalFailed() = default;

  explicit ModelEvalFailed(const std::string &_msg) {
    m_msg += _msg;
  }

  [[nodiscard]] const char *what() const noexcept override {
    return m_msg.c_str();
  }

 private:
  std::string m_msg{"*** relxill-error: "};
};



class LocalModel {

 public:
  LocalModel(const ModelParams &par, ModelName model_name)
      : m_name{model_name},
        m_param{par},
        m_info{ModelDatabase::instance().get_model_definition(model_name).model_info()}
  {  };

  LocalModel(double* inp_param, ModelName model_name)
      : LocalModel( ModelParams(ModelDatabase::instance().get_model_definition(model_name).parameter_names(), inp_param),
                    model_name )
  {  };

  explicit LocalModel(ModelName model_name) :
      LocalModel(ModelParams(), model_name)
  {  };

  void line_model(const XspecSpectrum &spectrum);
  void relxill_model(const XspecSpectrum &spectrum);
  void conv_model(const XspecSpectrum &spectrum);
  void xillver_model(const XspecSpectrum &spectrum);

  /**
   * Evaluate the LocalModel and overwrite the "flux" array
   * of Spectrum with the output values
   * @param spectrum
   */
  void eval_model(XspecSpectrum &spectrum) {

    switch (m_info.type()) {
      case T_Model::Line: line_model(spectrum);
        break;

      case T_Model::Relxill: relxill_model(spectrum);
        break;

      case T_Model::Conv: conv_model(spectrum);
        break;

      case T_Model::Xill: xillver_model(spectrum);
        break;
    }

  }

 private:
  ModelName m_name;
  ModelParams m_param;
  ModelInfo m_info;
};

void xspec_C_wrapper_eval_model(ModelName model_name,
                                const double *parameter_values,
                                double *xspec_flux,
                                int num_flux_bins,
                                const double *xspec_energy);




#endif //RELXILL__CPPMODELS_H_
