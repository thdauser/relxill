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

relxill::IrradiationType T_Irrad;
relxill::PrimarySpecType T_PrimeSpec;

extern "C" {
void lmodcpprelxill(const Array &energy, const Array &parameter,
                    int spectrum, Array &flux, Array &fluxError,
                    const string &init) {

  const std::vector<std::string> names = {"a", "Rin"};

  relxill::eval_model_xspec(relxill::RelxillModel::Relxill, energy, flux, parameter, names);
  // Model code:  Should resize flux RealArray to energy.size()-1.
  // Do the same for fluxError array if calculating errors, otherwise
  // leave it at size 0.
}

}