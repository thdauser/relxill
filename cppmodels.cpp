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



//class Global {
// public:
//  Global(){
//    m_models[ModelName::relline] = ModelInfo{T_Model::LineModel, T_Irrad::BknPowerlaw, T_PrimSpec::CutoffPl};
//
//  }
//
// private:
//  const std::map<ModelName, ModelInfo> m_models;
//}
//        {ModelName::relxill, {T_Model::RelxillModel, T_Irrad::BknPowerlaw, T_PrimSpec::CutoffPl}},



std::map<ModelName, ModelInfo> models;

extern "C" {
void lmodcpprelxill(const Array &energy, const Array &parameter,
                    int spectrum, Array &flux, Array &fluxError,
                    const string &init) {

  eval_model_xspec(ModelName::relxill, energy, flux, parameter);

}

}
