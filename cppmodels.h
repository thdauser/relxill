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
#include <vector>
#include <xsTypes.h>

extern "C" {
#include "relbase.h"
}

typedef RealArray Array;

namespace relxill {

//class ModelType{};
//class LineModel: ModelType{};
//class ConvModel: ModelType{};
//class XillModel: ModelType{};
//class RelxillModel: ModelType{};

enum class RelxillModel {
  Relxill,
  RelxillLp
} RelxillModel;

enum class LineModel {
  Relline,
  RellineLp
} LineModel;

enum IrradiationType {
  BknPowerlaw,
  LampPost,
  BlackBody
};

enum PrimarySpecType {
  CutoffPl,
  Nthcomp,
  Blackbody
};

class Param {
 public:
  Param(double inp_val, std::string inp_name) :
      val{inp_val}, name{inp_name} {
  }

 private:
  double val;
  std::string name{};
};

class RelParameters {

 private:
  Param a{0.998, "a"};
  Param rin{-1, "Rin"};
};

class Parameters : RelParameters {
 public:
  Parameters(Array values, std::vector<std::string> pnames) {
    // ....
  }

 private:
  std::map<std::string, double> param{};

};

class LocalModel {
 private:
  std::string name{};
  // Parameters par{};

 public:
  void eval_model();

};

class Spectrum {
 public:
  Spectrum(const Array energy, const Array flu) :
      energyArray{energy}, fluxArray{flu} {
  }

  Array energy() {
    return energyArray;
  }

  Array flux() {
    return fluxArray;
  }

 private:
  const Array energyArray;
  const Array fluxArray;
};

template<typename T_Model>
//, typename T_Irradiation, typename T_PrimeSpec>
void eval_model_xspec(T_Model model, std::string name, const Array &energy, const Array &flux,
                      const Array &parameter, std::vector<std::string> par_names) {

  Spectrum spec{energy, flux};
  Parameters pars{parameter, par_names};

  if (typeid(model) == typeid(LineModel)) {
    puts("I am LINE ");
  } else if (typeid(model) == typeid(RelxillModel)) {
    puts(" I am RELXILL ");
  }

  puts(" not doing anything ");
}

}

#endif //RELXILL__CPPMODELS_H_
