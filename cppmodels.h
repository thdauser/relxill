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

//extern "C" {
//#include "relbase.h"
//}

#include "cppparameters.h"

// namespace relxill {


class ModelType {

 public:
  ModelType(T_Model type, T_Irrad irrad, T_PrimSpec prim)
      : m_type{type}, m_irradiation{irrad}, m_primeSpec{prim} {
  };

  [[nodiscard]] T_Model type() const {
    return m_type;
  }

  [[nodiscard]] bool is_model_type(T_Model const model_type) const {
    return typeid(model_type) == typeid(m_type);
  }

 private:
  const T_Model m_type;
  const T_Irrad m_irradiation;
  const T_PrimSpec m_primeSpec;
};

// TODO: make the inheritance the other way around (no need for parameters)
class ModelInfo : ModelType {
 public:
  ModelInfo(ModelParamVector _par, T_Model type, T_Irrad irrad, T_PrimSpec prim)
      : m_param{std::move(_par)}, ModelType(type, irrad, prim) {
  }

  const ModelParamVector &parameters() const {
    return m_param;
  }

 private:
  ModelParamVector m_param{};
};

class ModelDatabase {

 public:
  static ModelDatabase &instance() {
    static auto *instance = new ModelDatabase();
    return *instance;
  }

  ModelInfo getModelInfo(ModelName name) {
    return m_models.at(name);
  }

  std::unordered_map<ModelName, ModelInfo> all() const {
    return m_models;
  }

 private:
  ModelDatabase() = default;  //hidden constructor and destructor
  ~ModelDatabase() = default;
  const std::unordered_map<ModelName, ModelInfo> m_models =
      {
          {ModelName::relline,
           ModelInfo({
                         XPar::linee,
                         XPar::index1,
                         XPar::index2,
                         XPar::rbr,
                         XPar::a,
                         XPar::incl,
                         XPar::rin,
                         XPar::rout,
                         XPar::z,
                         XPar::limb
                     },
                     T_Model::Line, T_Irrad::BknPowerlaw, T_PrimSpec::None)
          },

          {ModelName::relxill,
           ModelInfo({
                         XPar::index1,
                         XPar::index2,
                         XPar::rbr,
                         XPar::a,
                         XPar::incl,
                         XPar::rin,
                         XPar::rout,
                         XPar::z,
                         XPar::gamma,
                         XPar::logxi,
                         XPar::afe,
                         XPar::ecut,
                         XPar::refl_frac
                     },
                     T_Model::Relxill, T_Irrad::BknPowerlaw, T_PrimSpec::CutoffPl)},

          {ModelName::relconv,
           ModelInfo({
                         XPar::index1,
                         XPar::index2,
                         XPar::rbr,
                         XPar::a,
                         XPar::incl,
                         XPar::rin,
                         XPar::rout,
                         XPar::z,
                         XPar::gamma
                     },
                     T_Model::Conv, T_Irrad::BknPowerlaw, T_PrimSpec::CutoffPl)},

          {ModelName::xillver,
           ModelInfo({
                         XPar::index1,
                         XPar::index2,
                         XPar::rbr,
                         XPar::a,
                         XPar::incl,
                         XPar::rin,
                         XPar::rout,
                         XPar::z,
                         XPar::gamma,
                         XPar::logxi,
                         XPar::afe,
                         XPar::ecut,
                         XPar::refl_frac
                     },
                     T_Model::Xill, T_Irrad::BknPowerlaw, T_PrimSpec::CutoffPl)},
      };

};

class LocalModel {

  LocalModel(const ParamMap &par, const ModelInfo &info)
      : m_param{par}, m_info{info} {
  };

 private:
  ParamMap m_param;
  const ModelInfo m_info;
};

void eval_model_xspec(ModelName model, const Array &energy, const Array &flux, const Array &param_values);

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
