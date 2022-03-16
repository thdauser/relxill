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

    Copyright 2022 Thomas Dauser, Remeis Observatory & ECAP
*/

#ifndef RELXILL_SRC_MODELPARAMS_H_
#define RELXILL_SRC_MODELPARAMS_H_

#include <string>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <utility>

#include "XspecSpectrum.h" //only to get the typedef of Array
#include "ModelInfo.h"
#include "common.h"

#include <vector>

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
  iongrad_index,
  switch_switch_reflfrac_boost,
  switch_switch_returnrad,
  switch_iongrad_type,
  shifttmaxrrad  // only for testing relxillBB
};


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


// *** Function Definitions ***
int convertModelType(ModelName name);
int convertIrradType(T_Irrad name);
int convertPrimSpecType(T_PrimSpec name);

//relParam *getRelParamStruct(const ModelParams &params, ModelName model_name, ModelInfo model_info);
//xillParam *getXillParamStruct(const ModelParams &params, ModelName model_name, ModelInfo model_info);
const double *get_xspec_default_parameter_array(ModelName model_name);



/**
 * @brief class storing the name and the parameter list
 * - used by the python script from the lmodel.dat file to create the xspec_wrapper
 */

/**
 * class to store all input parameters of the model (explicit or hidden)
 */
class ParamList {

 public:

  explicit ParamList( std::vector<XPar> parnames, std::vector<double> parvalues) :
      m_parnames{std::move(parnames)}
  {

    for (size_t ii = 0; ii < m_parnames.size() ; ++ii){
      m_param.insert(std::make_pair(m_parnames[ii],parvalues[ii]));
    }
  }

  void set(XPar name, double value){
    try {
      m_param.at(name) = value;
    } catch (std::out_of_range &e){
      throw ParamInputException("parameter not found");
    }
  }

  auto &operator[](const XPar &name) const {
    try {
      return m_param.at(name);
    } catch (std::exception &e){
      throw ParamInputException("parameter not found");
    }
  }

  auto num_params(){
    return m_parnames.size();
  }

  auto get_parnames(){
    return m_parnames;
  }

  /**
   * get the parameter value for "name", otherwise return the default
   * value "def_value"
   * @param name
   * @param def_value
   * @return
   */
  double get_otherwise_default(const XPar &name, double def_value) const {
    if (m_param.find(name) != m_param.end()){
      return m_param.at(name);
    } else {
      return def_value;
    }
  }

 private:
  std::vector<XPar> m_parnames = {};  // need this as this needs to be in order given by Xspec TODO: could be fixed in xspec_wrapper
  std::unordered_map<XPar, double> m_param = {};
};

/* class InputParamList: public ParamList{

}; */

// class Model


#endif //RELXILL_SRC_MODELPARAMS_H_
