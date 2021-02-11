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
#include "relxill.h"
}

#include "ModelDatabase.h"
#include "ModelParams.h"

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
        m_model_params{par},
        m_info{ModelDatabase::instance().model_info(model_name)}
  {  };

  LocalModel(const double* inp_param, ModelName model_name)
      : LocalModel( ModelParams(ModelDatabase::instance().param_names(model_name), inp_param),
                    model_name )
  {  };

  explicit LocalModel(ModelName model_name) :
      LocalModel(ModelParams(), model_name)
  {  };

  /**
   * @brief set the value of a single parameter
   * @param XPar param
   * @param double value
   */
  void set_par(const XPar param, double value){
    m_model_params.set(param,value);
  }

  /**
   * Evaluate the LocalModel (in the Rest Frame of the Source)
   * (applies the redshift to the energy grid)
   * @param spectrum
   * @output spectrum.flux
   */
  void eval_model(XspecSpectrum &spectrum) {

    spectrum.shift_energy_grid_redshift(m_model_params[XPar::z]);

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
  ModelParams m_model_params;
  ModelInfo m_info;

  void line_model(const XspecSpectrum &spectrum);
  void relxill_model(const XspecSpectrum &spectrum);
  void conv_model(const XspecSpectrum &spectrum);
  void xillver_model(const XspecSpectrum &spectrum);

};

void xspec_C_wrapper_eval_model(ModelName model_name,
                                const double *parameter_values,
                                double *xspec_flux,
                                int num_flux_bins,
                                const double *xspec_energy);




#endif //RELXILL__CPPMODELS_H_
