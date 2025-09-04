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

    Copyright 2024 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef JEDSAD_H_
#define JEDSAD_H_

#include <utility>
#include <array>
#include <numeric>
#include "ModelDefinition.h"

extern "C" {
#include "relutility.h"
}

enum JSPar{
  mass,
  rj,
  mdot,
  risco,
  omega,
  mu,
  ms,
  b,
  p
};

const int num_params_jedsad_table = 9;

class JedsadTableInformation{

 public:
  /**
 * use this instance to always and securely access the database, without
 * needing to initialize it anywhere else
 */
  static JedsadTableInformation &instance() {
    static JedsadTableInformation inst;
    return inst;
  }

  // this is the order and name of the parameters in the table
  const std::vector<JSPar> param_names_table = {
      JSPar::mass,
      JSPar::rj,
      JSPar::mdot,
      JSPar::risco,
      JSPar::omega,
      JSPar::mu,
      JSPar::ms,
      JSPar::b,
      JSPar::p
  };

  string get_filename(){
    return m_filename;
  }

  const int n_data = 2;

 private:
  const string m_filename = "Refl_XrB_DauserTest_lb3.fits";  // TODO: make it readable by env variable

};


class JedsadParams{

 public:
  explicit JedsadParams(const ModelDefinition& _model_params){
    jedsad_params[JSPar::mass] = _model_params[XPar::mass];
    jedsad_params[JSPar::rj] = _model_params[XPar::rj];
    jedsad_params[JSPar::mdot] = _model_params[XPar::mdot];
    jedsad_params[JSPar::risco] = _model_params[XPar::rin]; // this jedsad parameter is simply Rg=GM/c^2
    jedsad_params[JSPar::omega] = _model_params[XPar::omega];
    jedsad_params[JSPar::mu] = _model_params[XPar::mu];
    jedsad_params[JSPar::ms] = _model_params[XPar::ms];
    jedsad_params[JSPar::b] = _model_params[XPar::b];
    jedsad_params[JSPar::p] = _model_params[XPar::p];
  }

  double* get_param_array(){
    auto array = new double[jedsad_params.size()];
    for (int ii=0; ii<jedsad_params.size(); ii++){
      // use the correct order as defined in the Table (in JedsadTableInformation)
      array[ii] = jedsad_params[JedsadTableInformation::instance().param_names_table[ii]];
    }
    return array;
  }

  std::unordered_map<JSPar, double> jedsad_params;
};

class PrimespecParams { ;

 public:
  PrimespecParams(double _gamma, double _ecut) : gamma{_gamma}, ecut{_ecut} {};

  double gamma;
  double ecut;
};


class JedsadTable{

 public:

  size_t num_params() {
    return param_names_table.size();
  };

  // Singleton Pattern
  static JedsadTable &instance(){
    static JedsadTable tab;
    return tab;
  }
  void read_table();

  std::vector<double> interpolate(const double *param_array);

  int get_data_index(std::vector<int> pind);


  int *num_param_vals;
  std::vector<double*> param_vals;
  std::vector<double *> raw_data;

  void update_filename(){
    // check env variable (TBD)
    // if different, update it
  }

  ~JedsadTable() {
    for (auto entry1: param_vals) {
      delete[] entry1;
    }

    for (auto entry2: raw_data)
      delete[] entry2;
    delete[] num_param_vals;
  }

  [[nodiscard]] int num_data_vals() const{
    return m_num_data_vals;
  }

 private:
  string m_fullfilename;
  const std::vector<JSPar> param_names_table = JedsadTableInformation::instance().param_names_table;
  void read_jstable_params(fitsfile *fptr);
  int m_num_data_vals = 0;

  JedsadTable(){
    m_fullfilename = JedsadTableInformation::instance().get_filename();  // the standard filename
    const int _num_params = num_params();
    assert(_num_params == num_params_jedsad_table);
    num_param_vals = new int[_num_params];
    read_table();
  };
  JedsadTable( const JedsadTable& );

} ;




PrimespecParams convert_jedsad_to_primaryspec_params(const ModelDefinition& model_definition);


#endif /* JEDSAD_H_ */