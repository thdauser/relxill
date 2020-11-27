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

class ParamMap {

 public:
  ParamMap() = default;
  ~ParamMap() = default;

  ParamMap(ModelParamVector pars, Array values) {
    assert(pars.size() == values.size());
    for (int ii = 0; ii < pars.size(); ii++) {
      m_param[pars[ii]] = values[ii];
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
      {XPar::a, 0.998}
  };

};

//
//class RelParameter{
//
//  void linee(double _v){  m_a = _v;  }
//  void a(double _v){  m_a = _v;  }
//  void a(double _v){  m_a = _v;  }
//
//  void set(XPar pname, double val){
//    switch(pname){
//
//      case XPar::linee:     m_linee=val;     break;
//      case XPar::index1:    m_index1=val;    break;
//      case XPar::index2:    m_index2=val;    break;
//      case XPar::rbr:       m_rbr=val;       break;
//      case XPar::a:         m_a=val;         break;
//      case XPar::rin:       m_rin=val;       break;
//      case XPar::rout:      m_rout=val;      break;
//      case XPar::incl:      m_incl=val;      break;
//      case XPar::z:         m_z=val;         break;
//      case XPar::gamma:     m_gamma=val;     break;
//      case XPar::logxi:     m_logxi=val;     break;
//      case XPar::afe:       m_afe=val;       break;
//      case XPar::ecut:      m_ecut=val;      break;
//      case XPar::refl_frac: m_refl_frac=val; break;
//
//      case XPar::limb:      m_limb=static_cast<int>(val);     break;
//
//    }
//
//  }
//
// private:
//  double m_linee{1.0};
//  double m_a{0.988};
//  double m_index1{3.0};
//  double m_index2{3.0};
//  double m_rin{-1.0};
//  double m_rbr{10.0};
//  double m_rout{400.0};
//  double m_incl{40.0};
//  double m_z{0.0};
//  double m_gamma{2.0};
//  double m_logxi{3.0};
//  double m_afe{1.0};
//  double m_ecut{100.0};
//  double m_refl_frac{1.0};
//
//  int    m_limb{0};
//
//};

//class Parameters : RelParameter  {
// public:
//  Parameters() = default;
//  Parameters(Array values, std::vector<ModelName> pnames) {
//  }
//
//
//
//
// private:
//  std::map<std::string, double> param{
//      {"a", 0.998}
//
//  };
//
//};
//


//
//class Spectrum {
// public:
//  Spectrum(const Array energy, const Array flu) :
//      energyArray{energy}, fluxArray{flu} {
//  }
//
//  Array energy() {
//    return energyArray;
//  }
//
//  Array flux() {
//    return fluxArray;
//  }
//
// private:
//  const Array energyArray;
//  const Array fluxArray;
//};

#endif //RELXILL_SRC_CPPPARAMETERS_H_
