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

// #include "cppparameters.h"

// namespace relxill {

enum class ModelName {
  relline,
  relxill,
  relxilllp
};

enum class XPar {
  linee,
  index1,
  index2,
  rbr,
  a,
  rin,
  rout,
  incl,
  z,
  limb,
  gamma,
  logxi,
  afe,
  ecut,
  refl_frac
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

typedef std::vector<std::string> StringVector;
typedef std::vector<XPar> ModelParamVector;
// typedef RealArray Array;

typedef std::valarray<double> Array;
typedef std::string string;

class ModelInfo {
 public:
  ModelInfo(ModelParamVector par, T_Model type, T_Irrad irrad, T_PrimSpec prim)
      : m_param{std::move(par)}, m_type{type}, m_irradiation{irrad}, m_primeSpec{prim} {
  }

  void all_models() {
  }

  [[nodiscard]] T_Model type() const {
    return m_type;
  }

  [[nodiscard]] bool is_model_type(T_Model const model_type) const {
    return typeid(model_type) == typeid(m_type);
  }

 private:
  ModelParamVector m_param{};
  T_Model m_type{};
  T_Irrad m_irradiation{};
  T_PrimSpec m_primeSpec{};
};


//class XspecInputParameters {
//
// public:
//  static XspecInputParameters &instance() {
//    static auto *instance = new XspecInputParameters();
//    return *instance;
//  }
//
//  std::vector<XPar> get(ModelName name) const {
//    return m_param.at(name);
//  }
//
//  void add(ModelName name, std::vector<XPar> params){
//    m_param.try_emplace(name, params);
//  }
//
// private:
//  XspecInputParameters() = default;  //hidden constructor and destructor
//  ~XspecInputParameters() = default;
//  std::unordered_map<ModelName, ModelParamVector> m_param;
//};


//class XspecParamStrings {
//
// public:
//  static XspecParamStrings &instance() {
//    static auto *instance = new XspecParamStrings();
//    return *instance;
//  }
//
//  std::string get(InputParam id) const {
//    return m_strings.at(id);
//  }
//
//
// private:
//  XspecParamStrings() = default;  //hidden constructor and destructor
//  ~XspecParamStrings() = default;
//  std::unordered_map<InputParam, std::string> m_strings ={
//      {InputParam::index1, "Index1"}
//
//  };
//};



class LocalModels {

 public:
  static LocalModels &instance() {
    static auto *instance = new LocalModels();
    return *instance;
  }

  ModelInfo getModelInfo(ModelName name) {
    //if (std::any_of(m_models.begin(), m_models.end(), ))
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
                     T_Model::LineModel, T_Irrad::BknPowerlaw, T_PrimSpec::None)
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
                     T_Model::RelxillModel, T_Irrad::BknPowerlaw, T_PrimSpec::CutoffPl)}
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


//class LocalModel {
//
// public:
//  LocalModel(ModelName name, std::vector<std::string> par_names) {
//    model = ModelInfo(name);
//  }
//
// private:
//  ModelInfo model;
//  Parameters par{};
//
//};

void eval_model_xspec(ModelName model, const Array &energy, const Array &flux, const Array &parameter);

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
