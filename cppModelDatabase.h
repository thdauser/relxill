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

#ifndef RELXILL__CPPMODELDATABASE_H_
#define RELXILL__CPPMODELDATABASE_H_

#include <typeinfo>
#include <utility>
#include "cppparameters.h"
#include "cppTypes.h"

/**
 * contains information to define one local model fully:
 *  - the input parameter array
 *  - the type of the model (class "ModelType")
 *
 */
class ModelDefinition {

 public:
  ModelDefinition(ModelParamVector _par, ModelInfo _type)
      : m_param{std::move(_par)}, m_model{_type} {
  }

  [[nodiscard]] const ModelParamVector &input_parameters() const {
    return m_param;
  }

  [[nodiscard]] ModelInfo model_info() const {
    return m_model;
  }

 private:
  ModelParamVector m_param{};
  const ModelInfo m_model;
};

/**
 * not a real class, but rather a global database (container) for all possible
 * local models with, using instances of the class "ModelDefinition":
 *  - input parameters as given in the order in the lmodel.dat file
 *  - specify the type of the model, including it's relevant physical components
 */
class ModelDatabase {

 public:
  static ModelDatabase &instance() {
    static auto *instance = new ModelDatabase();
    return *instance;
  }

  ModelDefinition get(ModelName name) {
    return m_models.at(name);
  }

  std::unordered_map<ModelName, ModelDefinition> all() const {
    return m_models;
  }

 private:
  ModelDatabase() = default;  //hidden constructor and destructor
  ~ModelDatabase() = default;

  const std::unordered_map<ModelName, ModelDefinition> m_models =
      {
          {ModelName::relline,
           ModelDefinition({
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
                           ModelInfo(T_Model::Line, T_Irrad::BknPowerlaw))
          },

          {ModelName::relxill,
           ModelDefinition({
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
                           ModelInfo(T_Model::Relxill, T_Irrad::BknPowerlaw, T_PrimSpec::CutoffPl)
           )},

          {ModelName::relconv,
           ModelDefinition({
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
                           ModelInfo(T_Model::Conv, T_Irrad::BknPowerlaw))
          },

          {ModelName::xillver,
           ModelDefinition({
                               XPar::gamma,
                               XPar::afe,
                               XPar::ecut,
                               XPar::logxi,
                               XPar::z,
                               XPar::incl,
                               XPar::refl_frac
                           },
                           ModelInfo(T_Model::Xill, T_PrimSpec::CutoffPl))
          },
      };

};

#endif //RELXILL__CPPMODELDATABASE_H_
