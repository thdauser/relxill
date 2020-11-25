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

typedef RealArray Array;

class Param {
 public:
  Param(double inp_val, std::string inp_name) :
      val{inp_val}, name{std::move(inp_name)} {
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
  Parameters() = default;
  Parameters(Array values, std::vector<std::string> pnames) {
    // ....
  }

 private:
  std::map<std::string, double> param{};

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

#endif //RELXILL_SRC_CPPPARAMETERS_H_
