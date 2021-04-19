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

#ifndef RELXILL__CPPTYPES_H_
#define RELXILL__CPPTYPES_H_

/**
 * stores all possible models
 */
enum class ModelName {
  relline,
  relline_lp,
  relconv,
  relconv_lp,
  relxill,
  relxillNS,
  relxillCO,
  relxillCp,
  relxillD,
  relxillDCp,
  relxilllp,
  relxilllpCp,
  relxilllpD,
  relxilllpDCp,
  relxilllpion,
  relxilllpionCp,
  relxilllpRet,
  rellinelpRet,
  xillver,
  xillverD,
  xillverCp,
  xillverDCp,
  xillverNS,
  xillverCO
};

enum class T_Irrad {
  BknPowerlaw,
  LampPost,
  BlackBody,
  Const,
  None
};

enum class T_PrimSpec {
  CutoffPl,
  Nthcomp,
  Blackbody,
  None
};

enum class T_Model {
  Line,
  Conv,
  Xill,
  Relxill
};

/**
 * @function holds information about the type of the model
 * @param T_Model, T_Irrad, T_PrimSpec
 *  - T_Model: line model, convolution, xillver,...
 *  - T_Irrad: type of irradiation (e.g., lamp post)
 *  - T_PrimSpec: spectral shape type of primary source
 */
class ModelInfo {

 public:
  ModelInfo(T_Model _type, T_Irrad _irrad, T_PrimSpec _prim)
      : m_type{_type}, m_irradiation{_irrad}, m_primeSpec{_prim} {
  };

  ModelInfo(T_Model type, T_Irrad irrad)
      : ModelInfo(type, irrad, T_PrimSpec::None) {
  };

  ModelInfo(T_Model type, T_PrimSpec prim)
      : ModelInfo(type, T_Irrad::None, prim) {
  };

  [[nodiscard]] T_Model type() const {
    return m_type;
  }

  [[nodiscard]] T_Irrad irradiation() const {
    return m_irradiation;
  }

  [[nodiscard]] T_PrimSpec primeSpec() const {
    return m_primeSpec;
  }

  [[nodiscard]] bool is_model_type(T_Model const model_type) const {
    return typeid(model_type) == typeid(m_type);
  }

 private:
  const T_Model m_type;
  const T_Irrad m_irradiation;
  const T_PrimSpec m_primeSpec;
};


/**
 * types to map to the previously defined C-types
 */


// std::map<Model


#endif //RELXILL__CPPTYPES_H_
