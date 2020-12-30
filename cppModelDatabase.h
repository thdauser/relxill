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

#include "xspec_wrapper_lmodels.h"


/**
 * @class ModelNotFound:
 * exception if anything is wrong with the input parameters
 */
class ModelNotFound : public std::exception {
 public:
  ModelNotFound() = default;

  explicit ModelNotFound(const std::string &_msg) {
    m_msg += _msg;
  }

  [[nodiscard]] const char *what() const noexcept override {
    return m_msg.c_str();
  }

 private:
  std::string m_msg{"*** model not found: "};
};

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
    try {
      return ModelDefinition(lmodel_database.params(name), lmodel_info.at(name));
    } catch (std::out_of_range &e) {
      throw ModelNotFound(lmodel_database.name_string(name));
    }
  }

  const XspecModelDatabase &database() const {
    return lmodel_database;
  }

 private:
  ModelDatabase() = default;  //hidden constructor and destructor to avoid initialization
  ~ModelDatabase() = default;

  XspecModelDatabase lmodel_database{}; // automatically created Class

  // TODO: make this an own class?
  // TODO: how to treat the high density table stuff
  const std::unordered_map<ModelName, ModelInfo> lmodel_info = {
      {ModelName::relline, ModelInfo(T_Model::Line, T_Irrad::BknPowerlaw)},
      {ModelName::rellinelp, ModelInfo(T_Model::Line, T_Irrad::LampPost)},

      {ModelName::relconvlp, ModelInfo(T_Model::Conv, T_Irrad::LampPost)},
      {ModelName::relconv, ModelInfo(T_Model::Conv, T_Irrad::BknPowerlaw)},

      {ModelName::relxill, ModelInfo(T_Model::Relxill, T_Irrad::BknPowerlaw, T_PrimSpec::CutoffPl)},
      {ModelName::relxillCp, ModelInfo(T_Model::Relxill, T_Irrad::BknPowerlaw, T_PrimSpec::Nthcomp)},
      {ModelName::relxillD, ModelInfo(T_Model::Relxill, T_Irrad::BknPowerlaw, T_PrimSpec::CutoffPl)},

      {ModelName::relxilllp, ModelInfo(T_Model::Relxill, T_Irrad::LampPost, T_PrimSpec::CutoffPl)},
      {ModelName::relxilllpCp, ModelInfo(T_Model::Relxill, T_Irrad::LampPost, T_PrimSpec::Nthcomp)},
      {ModelName::relxilllpD, ModelInfo(T_Model::Relxill, T_Irrad::LampPost, T_PrimSpec::CutoffPl)},
      {ModelName::relxilllpion, ModelInfo(T_Model::Relxill, T_Irrad::LampPost, T_PrimSpec::CutoffPl)},
      {ModelName::relxilllpionCp, ModelInfo(T_Model::Relxill, T_Irrad::LampPost, T_PrimSpec::Nthcomp)},

      {ModelName::xillver, ModelInfo(T_Model::Xill, T_PrimSpec::CutoffPl)},
      {ModelName::xillverD, ModelInfo(T_Model::Xill, T_PrimSpec::CutoffPl)},
      {ModelName::xillverCp, ModelInfo(T_Model::Xill, T_PrimSpec::Nthcomp)},

  };
};

#endif //RELXILL__CPPMODELDATABASE_H_
