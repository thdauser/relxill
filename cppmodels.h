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

  //  ModelType(T_Model type)
  //      : ModelType(type, T_Irrad::None, T_PrimSpec::None)
  //  {};

  ModelType(T_Model type, T_Irrad irrad)
      : ModelType(type, irrad, T_PrimSpec::None) {
  };

  ModelType(T_Model type, T_PrimSpec prim)
      : ModelType(type, T_Irrad::None, prim) {
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
class ModelDefinition {
 public:
  ModelDefinition(ModelParamVector _par, ModelType _type)
      : m_param{std::move(_par)}, m_model{_type} {
  }

  const ModelParamVector &input_parameters() const {
    return m_param;
  }

  const ModelType &model() {
    return m_model;
  }

 private:
  ModelParamVector m_param{};
  const ModelType m_model;
};

/*  == Class ModelDatabase ==
 *
 * not a real class, but rather a global database (container) for all possible
 * local models with:
 *  - input parameters as given in the order in the lmodel.dat file
 *  - specify the type of the model, including it's relevant physical components
 *
 */
class ModelDatabase {

 public:
  static ModelDatabase &instance() {
    static auto *instance = new ModelDatabase();
    return *instance;
  }

  ModelDefinition getInfo(ModelName name) {
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
                           ModelType(T_Model::Line, T_Irrad::BknPowerlaw))
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
                           ModelType(T_Model::Relxill, T_Irrad::BknPowerlaw, T_PrimSpec::CutoffPl)
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
                           ModelType(T_Model::Conv, T_Irrad::BknPowerlaw))
          },

          {ModelName::xillver,
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
                           ModelType(T_Model::Xill, T_PrimSpec::CutoffPl))
          },
      };

};

class LocalModel {

  LocalModel(const ParamMap &par, const ModelDefinition &info)
      : m_param{par}, m_info{info} {
  };

 private:
  ParamMap m_param;
  const ModelDefinition m_info;
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
