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

#ifndef RELXILL_SRC_CPPPARAMETERS_H_
#define RELXILL_SRC_CPPPARAMETERS_H_

#include <string>
#include <cassert>
#include <iostream>
#include <unordered_map>

#include "cppspectrum.h" //only to get the typedef of Array
#include "cppTypes.h"
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
  dens,
  afe,
  ecut,
  refl_frac,
  fixReflFrac,
  height,
  htop,
  beta,
  return_rad,
  frac_pl_bb,
  kTbb,
  ion_grad_type,
  ion_grad_index
};

typedef std::vector<XPar> ModelParamVector;

/**
 * exception if anything is wrong with the input parameters
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


class ModelParams {

 public:
  ModelParams() = default;
  ~ModelParams() = default;

  ModelParams(ModelParamVector pars, const double *values) {
    for (size_t ii = 0; ii < pars.size(); ii++) {
      m_param.at(pars[ii]) = values[ii];
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
      {XPar::dens, 15.0},
      {XPar::afe, 1.0},
      {XPar::refl_frac, 1.0},
      {XPar::fixReflFrac, -1},
      {XPar::limb, 0.0},
      {XPar::z, 0.0},
      {XPar::gamma, 2.0},
      {XPar::ecut, 300},
      {XPar::height, 6.0},
      {XPar::htop, 6.0},
      {XPar::beta, 0.0},
      {XPar::return_rad, 0.0},
      {XPar::frac_pl_bb, -1.0},
      {XPar::kTbb, -1.0},
      {XPar::ion_grad_type, ION_GRAD_TYPE_CONST},
      {XPar::ion_grad_index, 0.0}
  };

};

// *** Function Definitions ***
int convertModelType(ModelName name);
int convertIrradType(T_Irrad name);
int convertPrimSpecType(T_PrimSpec name);

relParam *getRelParamStruct(const ModelParams &params, ModelName model_name, ModelInfo model_info);
xillParam *getXillParamStruct(const ModelParams &params, ModelName model_name, ModelInfo model_info);

#endif //RELXILL_SRC_CPPPARAMETERS_H_
