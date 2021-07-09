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
  boost,
  h,
  htop,
  d_offaxis,
  beta,
  frac_pl_bb,
  ktbb,
  xi_index,
  switch_fixreflfrac,
  switch_return_rad,
  switch_ion_grad_type
};

typedef std::unordered_map<XPar, double> ParamList;

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
  LmodelParamList(std::string _name, ParamList _params)
      : m_name{std::move(_name)}, m_params{std::move(_params)} {
  };

  [[nodiscard]] std::string name() const {
    return m_name;
  }

  [[nodiscard]] ParamList params() const {
    return m_params;
  }

 private:
  std::string m_name;
  ParamList m_params;

};


// *** Function Definitions ***
int convertModelType(ModelName name);
int convertIrradType(T_Irrad name);
int convertPrimSpecType(T_PrimSpec name);

//relParam *getRelParamStruct(const ModelParams &params, ModelName model_name, ModelInfo model_info);
//xillParam *getXillParamStruct(const ModelParams &params, ModelName model_name, ModelInfo model_info);
const double *get_xspec_default_parameter_array(ModelName model_name);

#endif //RELXILL_SRC_MODELPARAMS_H_
