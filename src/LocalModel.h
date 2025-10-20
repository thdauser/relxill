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

#ifndef RELXILL__CPPMODELS_H_
#define RELXILL__CPPMODELS_H_

#include <utility>
#include <valarray>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C" {
#include "common.h"
}

#include "Relreturn_BlackBody.h"
#include "Relxill.h"
#include "ModelDatabase.h"
#include "ModelDefinition.h"

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

/*
explicit ParamList(ModelName model_name) :
ParamList(model_name, &ModelDatabase::instance().get_default_values_array(model_name)[0]) {
};

explicit ParamList(ModelName model_name, const double* parvalues) {
auto parnames = ModelDatabase::instance().param_names(model_name);
for (size_t ii = 0; ii < parnames.size() ; ++ii){
param.insert(std::make_pair(parnames[ii],parvalues[ii]));
}
}; */

/**
   * class LocalModel
   */
class LocalModel {

 public:
  LocalModel(const ParamList& par, ModelName model_name):
    m_model_params{ModelDefinition(par, model_name, ModelDatabase::instance().model_info(model_name))},
    m_rel_params{nullptr},
    m_xill_params{nullptr}
  {  };

    LocalModel(const double* inp_param, ModelName model_name)
        : LocalModel(ModelDatabase::instance().param_list(model_name), model_name )
    {
      set_input_params(inp_param);
    };

    explicit LocalModel(ModelName model_name) :
        LocalModel(ModelDatabase::instance().param_list(model_name), model_name)
    {  };

    // Destructor to clean up allocated parameters
  ~LocalModel()
  {
    if (m_rel_params != nullptr)
    {
      delete m_rel_params;
      m_rel_params = nullptr;
    }
    if (m_xill_params != nullptr)
    {
      delete m_xill_params;
      m_xill_params = nullptr;
    }
  }

  // Delete copy constructor and assignment operator to prevent double-free
  LocalModel(const LocalModel&) = delete;
  LocalModel& operator=(const LocalModel&) = delete;

  // Allow move operations
  LocalModel(LocalModel&& other) noexcept
    : m_model_params(std::move(other.m_model_params)),
      m_rel_params(other.m_rel_params),
      m_xill_params(other.m_xill_params)
  {
    other.m_rel_params = nullptr;
    other.m_xill_params = nullptr;
  }

  // Move assignment operator - swap idiom
  LocalModel& operator=(LocalModel&& other) noexcept
  {
    if (this != &other)
    {
      // Clean up existing resources
      delete m_rel_params;
      delete m_xill_params;

      m_model_params = std::move(other.m_model_params);
      m_rel_params = other.m_rel_params;
      m_xill_params = other.m_xill_params;

      // Reset other
      other.m_rel_params = nullptr;
      other.m_xill_params = nullptr;
    }
    return *this;
  }

  /** set the value of a single parameter
   *  this routine ensures that also xill_params and rel_params is updated
   * @param XPar param
   * @param double value
   */
  void set_par(const XPar param, double value)
  {
    m_model_params.set_par(param, value);

      // make sure the resulting rel and xill params are always updated and reflect the changes in m_model_params
    delete m_rel_params;
    m_rel_params = ::get_rel_params(m_model_params);

    delete m_xill_params;
    m_xill_params = ::get_xill_params(m_model_params);
  }

  std::string get_model_string()
  {
    return ModelDatabase::instance().model_string(m_model_params.get_model_name());
    }

    /**
     * Evaluate the LocalModel (in the Rest Frame of the Source)
     * (applies the redshift to the energy grid)
     * @param spectrum
     * @output spectrum.flux
     */
    void eval_model(XspecSpectrum &spectrum) {

      spectrum.shift_energy_grid_redshift(m_model_params.get_otherwise_default(XPar::z,0));

      try {
        switch (m_model_params.model_type()) {
          case T_Model::Line: line_model(spectrum);
            break;
          case T_Model::Relxill: relxill_model(spectrum);
            break;
          case T_Model::Conv: conv_model(spectrum);
            break;
          case T_Model::Xill: xillver_model(spectrum);
            break;
        }
      } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        throw ModelEvalFailed("model evaluation failed");
      }

    }

  relParam *rel_params() const {
    if (m_rel_params == nullptr)
    {
      m_rel_params = ::get_rel_params(m_model_params);
    }
    return m_rel_params;
  }

  xillParam* xill_params() const
  {
    if (m_xill_params == nullptr)
    {
      m_xill_params = ::get_xill_params(m_model_params);
    }
    return m_xill_params;
  }

  const ModelDefinition& get_model_params() const
  {
    return m_model_params;
  }

 private:
  ModelDefinition m_model_params;
  mutable relParam* m_rel_params; // mutable because automatically synchronized
  mutable xillParam* m_xill_params; // mutable because automatically synchronized

  void line_model(const XspecSpectrum& spectrum);
  void add_jedsad_parameters_to_model_params();
  bool handle_zero_jedsad_parameters(const XspecSpectrum& spectrum);
  static void check_relxill_error_default(int status);
  void relxill_model(const XspecSpectrum &spectrum);
  void conv_model(const XspecSpectrum &spectrum) const;
  void xillver_model(const XspecSpectrum &spectrum) const;

  void set_input_params(const double *inp_par_values) {
    auto parnames = m_model_params.get_parnames();
    for (size_t ii = 0; ii < m_model_params.num_params() ; ++ii){
      m_model_params.set_par(parnames[ii], inp_par_values[ii]);
    }
  }

};

void xspec_C_wrapper_eval_model(ModelName model_name,
                                const double *parameter_values,
                                double *xspec_flux,
                                int num_flux_bins,
                                const double *xspec_energy);




#endif //RELXILL__CPPMODELS_H_
