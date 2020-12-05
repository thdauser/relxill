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

#include "cppmodels.h"
#include <valarray>
#include <string>

// namespace relxill{

typedef std::valarray<double> Array;

class LmodTest {

 public:
  const Array energy{0.1, 1.0, 10.0};
  Array flux{0.0, 0.0, 0.0};
  Array empty{};
};

int main() {


  auto spec = LmodTest();
  int ispec = 0;
  const string empty_string;
  Array def_param;
  def_param.resize(10);

  lmodcpprelline(spec.energy, def_param, ispec, spec.flux, spec.empty, empty_string);

  return EXIT_SUCCESS;
}
// } // namespace relxill
