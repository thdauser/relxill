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

#ifndef RELXILL_SRC_MODELPARAMS_H_
#define RELXILL_SRC_MODELPARAMS_H_

#include <string>
#include <cassert>
#include <iostream>
#include <unordered_map>

#include "XspecSpectrum.h" //only to get the typedef of Array
#include "ModelInfo.h"
#include "common.h"

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
  logn,
  afe,
  a_co,
  ecut,
  kte,
  refl_frac,
  h,
  htop,
  d_offaxis,
  beta,
  return_rad,
  frac_pl_bb,
  ktbb,
  xi_index,
  switch_fixreflfrac,
  switch_ion_grad_type
};

typedef std::vector<XPar> ModelParamVector;

/**
 * @brief exception if anything is wrong with the input parameters
 */
class ParamInputException : public std::exception {
 public:
  ParamInputException() = default;

  explicit ParamInputException(const std::string &_msg) {
    m_msg += _msg;
  }

 private:
  [[nodiscard]] const char *what() const noexcept override {
    return m_msg.c_str();
  }

 private:
  std::string m_msg{"*** input parameter error: "};
};

/**
 * @brief class storing the name and the parameter list
 * - used by the python script from the lmodel.dat file to create the xspec_wrapper
 */
class LmodelParamList {

 public:
  LmodelParamList(std::string _name, ModelParamVector _params)
      : m_name{std::move(_name)}, m_params{std::move(_params)} {
  };

  [[nodiscard]] std::string name() const {
    return m_name;
  }

  [[nodiscard]] ModelParamVector params() const {
    return m_params;
  }

 private:
  std::string m_name;
  ModelParamVector m_params;

};

/**
 * @brief class to store all parameters of the model
 * @remark contains default values for all parameters,
 * which will get overwritten for the given input parameters.
 */
class ModelParams {

 public:
  ModelParams() = default;
  ~ModelParams() = default;

  /**
   * @param pars: list of parameters  (XPar)
   * @param values: input values corresponding to the parameter list
   */
  ModelParams(ModelParamVector pars, const double *values) {
    for (size_t ii = 0; ii < pars.size(); ii++) {
      m_param.at(pars[ii]) = values[ii];
    }
  }

  void set(XPar name, double value){
    try {
      m_param.at(name) = value;
    } catch (std::out_of_range &e){
      throw ParamInputException("parameter not found");
    }
  }

  auto &operator[](const XPar &name) const {
    return m_param.at(name);
  }

 private:
  std::unordered_map<XPar, double> m_param{
      {XPar::linee, 1.0},
      {XPar::index1, 3.0},
      {XPar::index2, 3.0},
      {XPar::a, 0.998},
      {XPar::rin, -1.0},
      {XPar::rbr, 30},
      {XPar::rout, 400},
      {XPar::incl, 30},
      {XPar::logxi, 3.0},
      {XPar::logn, 15.0},
      {XPar::afe, 1.0},
      {XPar::a_co, 1.0},
      {XPar::refl_frac, 1.0},
      {XPar::limb, 0.0},
      {XPar::z, 0.0},
      {XPar::gamma, 2.0},
      {XPar::ecut, 300},
      {XPar::kte, 100},
      {XPar::h, 3.0},
      {XPar::htop, 0.0},
      {XPar::beta, 0.0},
      {XPar::return_rad, 0.0},
      {XPar::frac_pl_bb, -1.0},
      {XPar::ktbb, -1.0},
      {XPar::xi_index, 0.0},
      {XPar::switch_fixreflfrac, -1},
      {XPar::switch_ion_grad_type, ION_GRAD_TYPE_CONST}
  };

};

// *** Function Definitions ***
int convertModelType(ModelName name);
int convertIrradType(T_Irrad name);
int convertPrimSpecType(T_PrimSpec name);

relParam *getRelParamStruct(const ModelParams &params, ModelName model_name, ModelInfo model_info);
xillParam *getXillParamStruct(const ModelParams &params, ModelName model_name, ModelInfo model_info);
const double *get_xspec_default_parameter_array(ModelName model_name);

#endif //RELXILL_SRC_MODELPARAMS_H_
