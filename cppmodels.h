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

#include <string>
#include <map>
#include <utility>
#include <vector>
#include <xsTypes.h>

extern "C" {
#include "relbase.h"
}

#include "src/cppparameters.h"

namespace relxill {

enum class ModelName {
  relline,
  relxill,
  relxilllp
};

enum class T_Model {
  LineModel,
  ConvModel,
  XillModel,
  RelxillModel
};

enum class T_Irrad {
  BknPowerlaw,
  LampPost,
  BlackBody,
  None
};

enum class T_PrimSpec {
  CutoffPl,
  Nthcomp,
  Blackbody,
  None
};

class ModelInfo {
 public:
  ModelInfo(T_Model type, T_Irrad irrad, T_PrimSpec prim)
      : m_type{type}, m_irradiation{irrad}, m_primeSpec{prim} {
  }

  void all_models() {

  }

 private:
  T_Model m_type{};
  T_Irrad m_irradiation{};
  T_PrimSpec m_primeSpec{};
};

typedef std::vector<std::string> StringVector;

// TODO:: I think this should be done in the "lmod" call (for everything else we can use something better)
class ModelParameters {

 public:
  static ModelParameters &instance() {
    static auto *instance = new ModelParameters();
    return *instance;
  }

  StringVector getModelInfo(ModelName name) {
    return m_param.at(name);
  }

  std::unordered_map<ModelName, StringVector> all() const {
    return m_param;
  }

 private:
  ModelParameters() = default;  //hidden constructor and destructor
  ~ModelParameters() = default;
  const std::unordered_map<ModelName, StringVector> m_param =
      {
          {ModelName::relline, {"a", "Rin", "Rout"}},
          {ModelName::relxill, {"a", "Rin", "Rout", "lxi", "refl_frac"}}
      };

};

class LocalModels {

 public:
  static LocalModels &instance() {
    static auto *instance = new LocalModels();
    return *instance;
  }

  ModelInfo getModelInfo(ModelName name) {
    return m_models.at(name);
  }

  std::unordered_map<ModelName, ModelInfo> all() const {
    return m_models;
  }

 private:
  LocalModels() = default;  //hidden constructor and destructor
  ~LocalModels() = default;
  const std::unordered_map<ModelName, ModelInfo> m_models =
      {
          {ModelName::relline, ModelInfo(T_Model::LineModel, T_Irrad::BknPowerlaw, T_PrimSpec::None)},
          {ModelName::relxill, ModelInfo(T_Model::RelxillModel, T_Irrad::BknPowerlaw, T_PrimSpec::CutoffPl)}
      };

};

//class LineModel: ModelInfo{
// private:
//  T_Irrad m_irradiation{};
//};
//class XillModel: ModelInfo{
// private:
//  T_PrimSpec m_primeSpec{};
//};
//class RelxillModel: LineModel, XillModel{ };
//class ConvModel: ModelInfo{};
//class XillModel: ModelInfo{};


class LocalModel {

 public:
  LocalModel(ModelName name, std::vector<std::string> par_names) {
    model = ModelInfo(name);
  }

 private:
  ModelInfo model;
  Parameters par{};

};

template<typename T_Model>
//, typename T_Irradiation, typename T_PrimeSpec>
void eval_model_xspec(T_Model model, const Array &energy, const Array &flux,
                      const Array &parameter, std::vector<std::string> par_names) {

  relxill::LocalModels::instance().getModelInfo(model);

  //lmod{model, parameter};

  Spectrum spec{energy, flux};
  Parameters pars{parameter, std::move(par_names)};



  //  eval_model(lmod, spec);


  if (typeid(model) == typeid(LineModel)) {
    puts("I am LINE ");
  } else if (typeid(model) == typeid(RelxillModel)) {
    puts(" I am RELXILL ");
  }

  puts(" not doing anything ");
}

}
#endif //RELXILL__CPPMODELS_H_
