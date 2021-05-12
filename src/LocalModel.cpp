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

#include "LocalModel.h"
#include "XspecSpectrum.h"

#include <stdexcept>
#include <iostream>

/**
 * @brief get a new RELATIVISITC PARAMETER STRUCTURE and initialize it with DEFAULT VALUES
 */
relParam *LocalModel::get_rel_params() {
  auto *param = new relParam;

  param->model_type = convertModelType(m_name);
  param->emis_type = convertIrradType(m_info.irradiation());

  param->a = m_model_params[XPar::a];
  param->incl = m_model_params[XPar::incl] * M_PI / 180;  // conversion to rad is heritage from the old code
  param->emis1 = m_model_params[XPar::index1];
  param->emis2 = m_model_params[XPar::index2];
  param->rbr = m_model_params[XPar::rbr];
  param->rin = m_model_params[XPar::rin];
  param->rout = m_model_params[XPar::rout];
  param->lineE = m_model_params[XPar::linee];
  param->z = m_model_params[XPar::z];
  param->height = m_model_params[XPar::h];
  param->gamma = m_model_params[XPar::gamma];
  param->beta = m_model_params[XPar::beta];
  param->htop = m_model_params[XPar::htop];
  param->limb = static_cast<int>(lround(m_model_params[XPar::limb]));
  param->return_rad = static_cast<int>(lround(m_model_params[XPar::return_rad]));
  param->return_rad_flux_correction_factor = 1.0; // needs to be calculated in the code
  param->xillver_gshift_corr_fac = 1.0; // needs to be calculated in the code

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
 * @brief get a new XILLVER PARAMETER STRUCTURE and initialize it with DEFAULT VALUES
 */
xillParam *LocalModel::get_xill_params() {
  auto *param = new xillParam;

  param->model_type = convertModelType(m_name);
  param->prim_type = convertPrimSpecType(m_info.primeSpec());

  param->gam = m_model_params[XPar::gamma];
  param->afe = m_model_params[XPar::afe];
  param->lxi = m_model_params[XPar::logxi];
  param->ect = (m_info.primeSpec() == T_PrimSpec::Nthcomp) ? m_model_params[XPar::kte]
                                                           : m_model_params[XPar::ecut];  // TODO: make kTe own parameter
  param->dens = m_model_params[XPar::logn];
  param->incl = m_model_params[XPar::incl];
  param->z = m_model_params[XPar::z];
  param->refl_frac = m_model_params[XPar::refl_frac];
  param->fixReflFrac = static_cast<int>(lround(m_model_params[XPar::switch_fixreflfrac]));
  param->frac_pl_bb = m_model_params[XPar::frac_pl_bb];
  param->kTbb = m_model_params[XPar::ktbb];
  param->ion_grad_type = static_cast<int>(lround(m_model_params[XPar::switch_ion_grad_type]));
  param->ion_grad_index = m_model_params[XPar::xi_index];

  // special definition of the xillver-co table
  if (is_co_model(param->model_type)) {
    param->dens = 17;
    param->afe = m_model_params[XPar::a_co];
    param->lxi = 0.0;
  }

  return param;
}

/*
 * @brief calculate line model
 */
void LocalModel::line_model(const XspecSpectrum &spectrum) {

  relParam *rel_param = LocalModel::get_rel_params();

  // relline_base calculates the line for 1keV -> shift the energy grid accordingly
  spectrum.shift_energy_grid_1keV(rel_param->lineE);

  int status = EXIT_SUCCESS;

  relline_base(spectrum.energy(), spectrum.flux(), spectrum.num_flux_bins(), rel_param, &status);

  if (status != EXIT_SUCCESS) {
    throw std::exception();
  }

  delete rel_param;
}

/*
 * @brief calculate relxill model
 */
void LocalModel::relxill_model(const XspecSpectrum &spectrum) {

  relParam *rel_param = LocalModel::get_rel_params();
  xillParam *xill_param = LocalModel::get_xill_params();

  int status = EXIT_SUCCESS;
  relxill_kernel(spectrum.energy(), spectrum.flux(), spectrum.num_flux_bins(), xill_param, rel_param, &status);

  if (status != EXIT_SUCCESS) {
    throw std::exception();
  }

  delete rel_param;
  delete xill_param;
}

/**
 * @brief calculate convolution of a given spectrum
 */
void LocalModel::conv_model(const XspecSpectrum &spectrum) {

  if (calcSum(spectrum.flux(), spectrum.num_flux_bins()) <= 0.0) {
    throw ModelEvalFailed("input flux for convolution model needs to be >0");
  }

  relParam *rel_param = LocalModel::get_rel_params();

  int status = EXIT_SUCCESS;
  relconv_kernel(spectrum.energy(), spectrum.flux(), spectrum.num_flux_bins(), rel_param, &status);

  if (status != EXIT_SUCCESS) {
    throw ModelEvalFailed("executing convolution failed");
  }

  delete rel_param;
}


/**
 * @brief calculate xillver model
 */
void LocalModel::xillver_model(const XspecSpectrum &spectrum) {

  xillParam *xill_param = LocalModel::get_xill_params();

  int status = EXIT_SUCCESS;
  xillver_base(spectrum.energy(), spectrum.num_flux_bins(), spectrum.flux(), xill_param, &status);

  if (status != EXIT_SUCCESS) {
    throw std::exception();
  }

  delete xill_param;
}


/**
 * Wrapper function, which can be called from any C-type local model function as
 * required by Xspec
 * @param model_name: unique name of the model
 * @param parameter_values: input parameter_values array (size can be determined from the model definition for the model_name)
 * @param xspec_energy[num_flux_bins]: input energy grid
 * @param xspec_flux[num_flux_bins]: - flux array (already allocated), used to return the calculated values
 *                    - for convolution models this is also the input flux
 * @param num_flux_bins
 */
void xspec_C_wrapper_eval_model(ModelName model_name,
                                const double *parameter_values,
                                double *xspec_flux,
                                int num_flux_bins,
                                const double *xspec_energy) {

  try {
    LocalModel local_model{parameter_values, model_name};

    XspecSpectrum spectrum{xspec_energy, xspec_flux, static_cast<size_t>(num_flux_bins)};
    local_model.eval_model(spectrum);

  } catch (ModelNotFound &e) {
    std::cout << e.what();
  } catch (ModelEvalFailed &e) {
    std::cout << e.what();
    // TODO: what should we do if the evaluation fails? return zeros?
  }

}
