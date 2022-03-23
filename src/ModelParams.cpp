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

#include "ModelParams.h"
#include "ModelDatabase.h"
#include "Relmodels.h"
int get_iongrad_type(const ModelParams &params);
extern "C" {
#include "relutility.h"
}

/**
 * @brief maps the new C++ model type definition to the C-integers
 *  // TODO: use ModelName types throughout the code, then this function is obsolete
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
    case ModelName::relxillAlpha: return MOD_TYPE_RELXILLALPHA;
    case ModelName::relxilllpAlpha: return MOD_TYPE_RELXILLLPALPHA;
    case ModelName::relxillBB: return MOD_TYPE_RELXILLBBRET;
  }
  puts(" *** relxill-error: unknown ModelName, converting model name to integer failed ");
  exit(EXIT_FAILURE);
}

/**
 * @brief maps the new C++ irrad type definition to the C-integers
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
 * @brief maps the new C++ primary spectrum type definition to the C-integers
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
 * @brief Return Default Parameters Array for a given Model
 * It will be return as double array, the same as given as input from Xspec
 * @param ModelName model_name
 * @return double param_array
 */
const double *get_xspec_default_parameter_array(ModelName model_name) {

  auto const default_values = ModelDatabase::instance().get_default_values_array(model_name);

  auto output_param_array = new double[default_values.size()];

  // TODO: Need to change this to the actual input parameters

  for (size_t ii = 0; ii < default_values.size(); ii++){
    output_param_array[ii] = default_values[ii];
  }

  return output_param_array;
}



static int get_returnrad_switch(const ModelParams& model_params){

  // by default activated for Lamp Post, otherwise noe
  int default_switch = ( model_params.irradiation()==T_Irrad::LampPost) ? 1 : 0;

  // if env is set, we use the value given there, it takes precedence over default values
  char *env = getenv("RELXILL_RETURNRAD_SWITCH");
  if (env != nullptr) {
    int value_switch = atof(env);
    if (value_switch == 1 ) {
      default_switch = value_switch;
    } else {
      default_switch = 0; // set to zero for any other value
    }
  }

  // model parameter values take precedence over env variable
  // TODO: warning if conflict with env variable (??)
  return static_cast<int>(lround(model_params.get_otherwise_default(XPar::switch_switch_returnrad,default_switch)));
}



/**
 * @brief get a new RELATIVISITC PARAMETER STRUCTURE and initialize it with DEFAULT VALUES
 */
relParam* get_rel_params(const ModelParams& inp_param) {
  auto *param = new relParam;

  param->model_type = convertModelType(inp_param.get_model_name());
  param->emis_type = convertIrradType(inp_param.irradiation());

  // these parameters have to be given for any relativistic parameter structure
  try {
    param->a = inp_param.get_par(XPar::a);
    param->incl = inp_param.get_par(XPar::incl) * M_PI / 180;  // conversion to rad is heritage from the old code
    param->rin = inp_param.get_par(XPar::rin);
    param->rout = inp_param.get_par(XPar::rout);
  } catch (ParamInputException &e) {
    throw ParamInputException("get_rel_params: model evaluation failed due to missing relativistic parameters");
  }

  // those values should never be used, unless it is set by the model
  param->emis1 = inp_param.get_otherwise_default(XPar::index1, 0);
  param->emis2 = inp_param.get_otherwise_default(XPar::index2, 0);
  param->rbr = inp_param.get_otherwise_default(XPar::rbr, 0);
  param->lineE = inp_param.get_otherwise_default(XPar::linee, 0);
  param->gamma = inp_param.get_otherwise_default(XPar::gamma, 0);
  param->height = inp_param.get_otherwise_default(XPar::h, 0);
  param->htop = inp_param.get_otherwise_default(XPar::htop, 0);

  // important default values
  param->z = inp_param.get_otherwise_default(XPar::z, 0);
  param->beta = inp_param.get_otherwise_default(XPar::beta, 0);
  param->limb = static_cast<int>(lround(inp_param.get_otherwise_default(XPar::limb, 0)));
  param->return_rad = get_returnrad_switch(inp_param );

  param->rrad_corr_factors = nullptr; //

  // this is set by the environment variable "RELLINE_PHYSICAL_NORM"
  param->do_renorm_relline = do_renorm_model(param);

  int status = EXIT_SUCCESS;
  check_parameter_bounds(param, &status);
  if (status != EXIT_SUCCESS) {
    puts(" *** relxill-error: problem interpreting the input parameter values");
    throw ParamInputException();
  }

  param->ion_grad_type = get_iongrad_type(inp_param);

  // set depending on model/emis type and ENV "RELXILL_NUM_RZONES"
  param->num_zones = get_num_zones(param->model_type, param->emis_type, get_iongrad_type(inp_param));

  return param;
}

/**
 * @brief get a new XILLVER PARAMETER STRUCTURE and initialize it with DEFAULT VALUES
 */
xillParam* get_xill_params(const ModelParams& inp_param) {
  auto *param = new xillParam;

  param->model_type = convertModelType(inp_param.get_model_name());
  param->prim_type = convertPrimSpecType(inp_param.primeSpec());

  // these parameters have to be given for any xillver parameter structure
  try {
    param->afe = (is_co_model(param->model_type))
                 ? inp_param.get_par(XPar::a_co)  // special definition of the xillver-co table
                 : inp_param.get_par(XPar::afe);
    param->incl = inp_param.get_par(XPar::incl);
    param->z = inp_param.get_par(XPar::z);

  } catch (ParamInputException &e) {
    throw ParamInputException("get_xill_params: model evaluation failed due to missing xillver parameters");
  }

  // important default values
  param->ect = (inp_param.primeSpec() == T_PrimSpec::Nthcomp)
               ? inp_param.get_otherwise_default(XPar::kte, 0)  // TODO: make kTe own parameter
               : inp_param.get_otherwise_default(XPar::ecut, 300);
  param->lxi = inp_param.get_otherwise_default(XPar::logxi, 0);  // default value for CO table
  param->dens = inp_param.get_otherwise_default(XPar::logn,               // CO-table has logN=17
                                                     is_co_model(param->model_type) ? 17 : 15);
  param->iongrad_index = inp_param.get_otherwise_default(XPar::iongrad_index, 0);
  param->boost = inp_param.get_otherwise_default(XPar::boost, -1);

  // those values should never be used, unless it is set by the model
  param->gam = inp_param.get_otherwise_default(XPar::gamma, 0);
  param->refl_frac = inp_param.get_otherwise_default(XPar::refl_frac, 0);
  param->frac_pl_bb = inp_param.get_otherwise_default(XPar::frac_pl_bb, 0);
  param->kTbb = inp_param.get_otherwise_default(XPar::ktbb, 0);

  param->interpret_reflfrac_as_boost =
      static_cast<int>(lround(inp_param.get_otherwise_default(XPar::switch_switch_reflfrac_boost, 0)));

  // to be deleted, only for testing
  param->shiftTmaxRRet = inp_param.get_otherwise_default(XPar::shifttmaxrrad, 0.0);

  return param;
}


int get_iongrad_type(const ModelParams &params) {
  if (params.get_model_name() == ModelName::relxilllpAlpha || params.get_model_name() == ModelName::relxillAlpha){
    return ION_GRAD_TYPE_ALPHA;
  } else {
    return static_cast<int>(lround(params.get_otherwise_default(XPar::switch_iongrad_type, 0)));
  }
}

