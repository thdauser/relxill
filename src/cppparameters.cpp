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

#include "cppparameters.h"
#include "cppModelDatabase.h"
extern "C" {
#include "relmodels.h"
#include "relutility.h"
}

/**
 * maps the new C++ model type definition to the C-integers
 *  - TODO: use ModelName types throughout the code, then this function is obsolete
 * @param name
 * @return (int) MODEL_TYPE
 */
int convertModelType(ModelName name) {
  switch (name) {
    case ModelName::relline: return MOD_TYPE_RELLINE;
    case ModelName::relconv: return MOD_TYPE_RELCONV;
    case ModelName::relline_lp: return MOD_TYPE_RELLINELP;
    case ModelName::relconv_lp: return MOD_TYPE_RELCONVLP;
    case ModelName::relxill  : return MOD_TYPE_RELXILL;
    case ModelName::relxillCp: return MOD_TYPE_RELXILL;
    case ModelName::relxillD: return MOD_TYPE_RELXILLDENS;
    case ModelName::relxilllp: return MOD_TYPE_RELXILLLP;
    case ModelName::relxilllpCp: return MOD_TYPE_RELXILLLP;
    case ModelName::relxilllpD: return MOD_TYPE_RELXILLLPDENS;
    case ModelName::relxilllpion  : return MOD_TYPE_RELXILLLPION;
    case ModelName::relxilllpionCp: return MOD_TYPE_RELXILLLPION;
    case ModelName::xillver: return MOD_TYPE_XILLVER;
    case ModelName::xillverCp: return MOD_TYPE_XILLVER_NTHCOMP;
    case ModelName::xillverD: return MOD_TYPE_XILLVERDENS;
    case ModelName::xillverNS: return MOD_TYPE_XILLVERNS;
    case ModelName::xillverCO: return MOD_TYPE_XILLVERCO;
    case ModelName::relxillNS: return MOD_TYPE_RELXILLNS;
    case ModelName::relxillCO: return MOD_TYPE_RELXILLCO;
    case ModelName::relxillDCp: return MOD_TYPE_RELXILLDENS_NTHCOMP;
    case ModelName::relxilllpDCp: return MOD_TYPE_RELXILLLPDENS_NTHCOMP;
    case ModelName::xillverDCp: return MOD_TYPE_XILLVERDENS_NTHCOMP;
  }
  puts(" *** relxill-error: unknown ModelName, converting model name to integer failed ");
  exit(EXIT_FAILURE);
}

/**
 * maps the new C++ irrad type definition to the C-integers
 * @param name
 * @return (int) MODEL_TYPE
 */
int convertIrradType(T_Irrad name) {
  switch (name) {
    case T_Irrad::BknPowerlaw: return EMIS_TYPE_BKN;
    case T_Irrad::LampPost: return EMIS_TYPE_LP;
    case T_Irrad::BlackBody: return EMIS_TYPE_ALPHA;
    case T_Irrad::Const: return EMIS_TYPE_CONST;
    case T_Irrad::None:puts(" *** relxill-error: not possible to construct a model with Irradiation-Type <None> ");
      break;
  }
  exit(EXIT_FAILURE);
}

/**
 * maps the new C++ primary spectrum type definition to the C-integers
 * @param name
 * @return (int) MODEL_TYPE
 */
int convertPrimSpecType(T_PrimSpec name) {
  switch (name) {
    case T_PrimSpec::CutoffPl: return PRIM_SPEC_ECUT;
    case T_PrimSpec::Nthcomp: return PRIM_SPEC_NTHCOMP;
    case T_PrimSpec::Blackbody: return PRIM_SPEC_BB;
    case T_PrimSpec::None:puts(" *** relxill-error: not possible construct a model with Primary-Spectrum-Type <None> ");
      break;
  }
  exit(EXIT_FAILURE);
}

/**
 * get a new RELATIVISITC PARAMETER STRUCTURE and initialize it with DEFAULT VALUES
 */
relParam *getRelParamStruct(const ModelParams &params, ModelName model_name, ModelInfo model_info) {
  auto *param = new relParam;

  param->model_type = convertModelType(model_name);
  param->emis_type = convertIrradType(model_info.irradiation());

  param->a = params[XPar::a];
  param->incl = params[XPar::incl] * M_PI / 180;  // conversion to rad is heritage from the old code
  param->emis1 = params[XPar::index1];
  param->emis2 = params[XPar::index2];
  param->rbr = params[XPar::rbr];
  param->rin = params[XPar::rin];
  param->rout = params[XPar::rout];
  param->lineE = params[XPar::linee];
  param->z = params[XPar::z];
  param->height = params[XPar::h];
  param->gamma = params[XPar::gamma];
  param->beta = params[XPar::beta];
  param->htop = params[XPar::htop];
  param->limb = static_cast<int>(lround(params[XPar::limb]));
  param->return_rad = static_cast<int>(lround(params[XPar::return_rad]));

  // this is set by the environment variable "RELLINE_PHYSICAL_NORM"
  param->do_renorm_relline = do_renorm_model(param);

  int status = EXIT_SUCCESS;
  check_parameter_bounds(param, &status);
  if (status != EXIT_SUCCESS) {
    puts(" *** relxill-error: problem interpreting the input parameter values");
    throw ParamInputException();
  }

  // set depending on model/emis type and ENV "RELXILL_NUM_RZONES"
  param->num_zones = get_num_zones(param->model_type, param->emis_type, ION_GRAD_TYPE_CONST);

  return param;
}


/**
 * get a new XILLVER PARAMETER STRUCTURE and initialize it with DEFAULT VALUES
 */
xillParam *getXillParamStruct(const ModelParams &params, ModelName model_name, ModelInfo model_info) {
  auto *param = new xillParam;

  param->model_type = convertModelType(model_name);
  param->prim_type = convertPrimSpecType(model_info.primeSpec());

  param->gam = params[XPar::gamma];
  param->afe = params[XPar::afe];
  param->lxi = params[XPar::logxi];
  param->ect = (model_info.primeSpec() == T_PrimSpec::Nthcomp) ? params[XPar::kte]
                                                               : params[XPar::ecut];  // TODO: make kTe own parameter
  param->dens = params[XPar::logn];
  param->incl = params[XPar::incl];
  param->z = params[XPar::z];
  param->refl_frac = params[XPar::refl_frac];
  param->fixReflFrac = static_cast<int>(lround(params[XPar::switch_fixreflfrac]));
  param->frac_pl_bb = params[XPar::frac_pl_bb];
  param->kTbb = params[XPar::ktbb];
  param->ion_grad_type = static_cast<int>(lround(params[XPar::switch_ion_grad_type]));
  param->ion_grad_index = params[XPar::xi_index];

  // special definition of the xillver-co table
  if(is_co_model(param->model_type)){
    param->dens = 17;
    param->afe = params[XPar::a_co];
    param->lxi = 0.0;
  }

  return param;
}


/**
 * Return Default Parameters Array for a given Model
 * It will be return as double array, the same as given as input from Xspec
 * @param ModelName model_name
 * @return double param_array
 */
const double *get_xspec_default_parameter_array(ModelName model_name) {

  auto const model_parameters = ModelDatabase::instance().get_model_definition(model_name).parameter_names();

  auto default_param_values = ModelParams();
  auto output_param_array = new double[model_parameters.size()];

  for (size_t ii = 0; ii < model_parameters.size(); ii++) {
    output_param_array[ii] = default_param_values[model_parameters[ii]];
  }

  return output_param_array;
}
