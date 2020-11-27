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
#include <xsTypes.h>
#include <cassert>
#include <iostream>

enum class ModelName {
  relline,
  relxill,
  relconv,
  xillver
};

enum class T_Model {
  Line,
  Conv,
  Xill,
  Relxill
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

typedef std::vector<std::string> StringVector;
typedef std::vector<XPar> ModelParamVector;
// typedef RealArray Array;

typedef std::valarray<double> Array;
typedef std::string string;

typedef RealArray Array;

class ModelParams {

 public:
  ModelParams() = default;
  ~ModelParams() = default;

  ModelParams(ModelParamVector pars, Array values) {
    std::cout << pars.size() << "--- " << values.size() << std::endl;
    assert(pars.size() == values.size());
    for (int ii = 0; ii < pars.size(); ii++) {
      m_param.at(pars[ii]) = values[ii];
    }
  }

  auto &operator[](const XPar &name) {
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
      {XPar::afe, 1.0},
      {XPar::refl_frac, 1.0},
      {XPar::limb, 0.0},
      {XPar::z, 0.0},
      {XPar::gamma, 2.0},
      {XPar::ecut, 300}
  };

};



#endif //RELXILL_SRC_CPPPARAMETERS_H_
