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

#include "test_relxill.h"

#include "writeOutfiles.h"




void set_std_param_relline_lp(double *inp_par) {
  inp_par[0] = 6.4;
  inp_par[1] = 3.0;
  inp_par[2] = 0.998;
  inp_par[3] = 30.0;
  inp_par[4] = -1.;
  inp_par[5] = 400.;
  inp_par[6] = 0.0;  // redshift
  inp_par[7] = 0.0;
  inp_par[8] = 2.0;  // gamma
}

relParam *init_par_relline_lp(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELLINELP, EMIS_TYPE_LP, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_RELLINELP);

  param->lineE = inp_par[0];
  param->height = inp_par[1];
  param->a = inp_par[2];
  param->incl = inp_par[3] * M_PI / 180;
  param->rin = inp_par[4];
  param->rout = inp_par[5];
  param->z = inp_par[6];
  param->limb = (int) (inp_par[7] + 0.5);
  param->gamma = inp_par[8];

  param->beta = 0.0;
  param->return_rad_flux_correction_factor = 1.0;
  param->xillver_gshift_corr_fac = 1.0;

  check_parameter_bounds(param, status);
  CHECK_STATUS_RET(*status, NULL);

  return param;
}


relParam *get_std_param_rellinelp(int *status) {
  int n_param = NUM_PARAM_RELLINELP;
  double inp_par[NUM_PARAM_RELLINELP];
  set_std_param_relline_lp(inp_par);
  return init_par_relline_lp(inp_par, n_param, status);
}


void set_std_param_xillver(double *inp_par) {
  inp_par[0] = 2.1;    // Gamma
  inp_par[1] = 1.0;    // Afe
  inp_par[2] = 300.0;  // Ecut
  inp_par[3] = 0.0;    // logxi
  inp_par[4] = 0.;     // redshift
  inp_par[5] = 45.0;   // inclination
  inp_par[6] = -1.0;   // refl. frac.
}

void set_std_param_xillver_dens(double *inp_par) {
  inp_par[0] = 2.1;    // Gamma
  inp_par[1] = 1.0;    // Afe
  inp_par[2] = 15.0;  // logN
  inp_par[3] = 0.0;    // logxi
  inp_par[4] = 0.;     // redshift
  inp_par[5] = 45.0;   // inclination
  inp_par[6] = -1.0;   // refl. frac.
}

void set_std_param_xillverns(double *inp_par) {
  inp_par[0] = 2.0;    // kTB
  inp_par[1] = 1.0;    // Afe
  inp_par[2] = 15.0;  // logN
  inp_par[3] = 1.0;    // logxi
  inp_par[4] = 0.;     // redshift
  inp_par[5] = 45.0;   // inclination
  inp_par[6] = -1.0;   // refl. frac.
}

void set_std_param_xillverco(double *inp_par) {
  inp_par[0] = 2.1;    // Gamma
  inp_par[1] = 5.0;    // A_CO
  inp_par[2] = 0.1;   // kTbb
  inp_par[3] = 0.1;   // frac_pl_bb
  inp_par[4] = 300.0;  // Ecut
  inp_par[5] = 0.;     // redshift
  inp_par[6] = 45.0;   // inclination
  inp_par[7] = -1.0;   // refl. frac.
}

void set_std_param_xillver_dens_nthcomp(double *inp_par) {
  inp_par[0] = 2.1;    // Gamma
  inp_par[1] = 1.0;    // Afe
  inp_par[2] = 60.0;  // kTe
  inp_par[3] = 15;  // logN
  inp_par[4] = 0.0;    // logxi
  inp_par[5] = 0.;     // redshift
  inp_par[6] = 45.0;   // inclination
  inp_par[7] = -1.0;   // refl. frac.
}

void set_std_param_xillver_nthcomp(double *inp_par) {
  inp_par[0] = 2.1;    // Gamma
  inp_par[1] = 1.0;    // Afe
  inp_par[2] = 60.0;  // kTe
  inp_par[3] = 0.0;    // logxi
  inp_par[4] = 0.;     // redshift
  inp_par[5] = 45.0;   // inclination
  inp_par[6] = -1.0;   // refl. frac.
}

void set_std_param_relline(double *inp_par) {
  inp_par[0] = 6.4;
  inp_par[1] = 3.0;
  inp_par[2] = 3.0;
  inp_par[3] = 15.0;
  inp_par[4] = 0.998;
  inp_par[5] = 30.0;
  inp_par[6] = -1.1;
  inp_par[7] = 1000.;
  inp_par[8] = 0.0;
  inp_par[9] = 1.0;
}

void set_std_param_relxilllpdens_nthcomp(double *inp_par) {
  inp_par[0] = 0.998; // a
  inp_par[1] = 60.0;  // incl
  inp_par[2] = -1.0;  // rin
  inp_par[3] = 400.;  // rout
  inp_par[4] = 6.0;   // height
  inp_par[5] = 0.0;   // htop
  inp_par[6] = 0.0;   // beta
  inp_par[7] = 2.1;   // pl Index
  inp_par[8] = 3.1;   // logxi
  inp_par[9] = 1.0;   // Afe
  inp_par[10] = 100.0; // kTe
  inp_par[11] = 15.0; // logN
  inp_par[12] = 3.0;   // refl_frac
  inp_par[13] = 0.0;   // fixReflFrac
  inp_par[14] = 0.0;    // redshift
}

xillParam *init_par_xillver(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVER, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVER);

  param->gam = inp_par[0];
  param->afe = inp_par[1];
  param->ect = inp_par[2];
  param->lxi = inp_par[3];
  param->dens = 15; // logN
  param->z = inp_par[4];
  param->incl = inp_par[5]; // is given in degrees !!
  param->refl_frac = inp_par[6];

  return param;
}

xillParam *init_par_xillver_nthcomp(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVER_NTHCOMP, PRIM_SPEC_NTHCOMP, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVER_NTHCOMP);

  param->gam = inp_par[0];
  param->afe = inp_par[1];
  param->ect = inp_par[2];   // is kTe internally
  param->lxi = inp_par[3];
  param->dens = 15; // logN
  param->z = inp_par[4];
  param->incl = inp_par[5]; // is given in degrees !!
  param->refl_frac = inp_par[6];

  return param;
}

xillParam *init_par_xillver_ns(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVERNS, PRIM_SPEC_BB, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVERNS);

  param->kTbb = inp_par[0];
  param->afe = inp_par[1];
  param->ect = 0.0;        // Ecut does not make sense for a BB spectrum
  param->lxi = inp_par[3];
  param->dens = inp_par[2]; // logN
  param->z = inp_par[4];
  param->incl = inp_par[5]; // is given in degrees !!
  param->refl_frac = inp_par[6];

  // TODO: check parameter bounds here as well
  /*	check_parameter_bounds_xillver(param,status);
      CHECK_STATUS_RET(*status,NULL); */
  
  return param;
}

xillParam *init_par_xillver_co(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVERCO, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVERCO);

  param->gam = inp_par[0];
  param->afe = inp_par[1]; // this is A_CO here for the xillverCO model
  param->kTbb = inp_par[2]; //
  param->frac_pl_bb = inp_par[3]; //
  param->ect = inp_par[4];
  param->z = inp_par[5];
  param->incl = inp_par[6]; // is given in degrees !!
  param->refl_frac = inp_par[7];

  param->dens = 17.0;
  param->lxi = 0.0;       // interestingly this model does not have an ionization

  // TODO: check parameter bounds here as well
  /*	check_parameter_bounds_xillver(param,status);
      CHECK_STATUS_RET(*status,NULL); */

  return param;
}

xillParam *init_par_xillver_dens(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVERDENS, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVERDENS);

  param->gam = inp_par[0];
  param->afe = inp_par[1];
  param->ect = 300.0;
  param->dens = inp_par[2];
  param->lxi = inp_par[3];
  param->z = inp_par[4];
  param->incl = inp_par[5]; // is given in degrees !!
  param->refl_frac = inp_par[6];

  // TODO: check parameter bounds here as well
  /*	check_parameter_bounds_xillver(param,status);
      CHECK_STATUS_RET(*status,NULL); */

  return param;
}

xillParam *init_par_xillver_dens_nthcomp(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVERDENS_NTHCOMP, PRIM_SPEC_NTHCOMP, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVERDENS_NTHCOMP);

  param->gam = inp_par[0];
  param->afe = inp_par[1];
  param->ect = inp_par[2];
  param->dens = inp_par[3];
  param->lxi = inp_par[4];
  param->z = inp_par[5];
  param->incl = inp_par[6]; // is given in degrees !!
  param->refl_frac = inp_par[7];

  // TODO: check parameter bounds here as well
  /*	check_parameter_bounds_xillver(param,status);
      CHECK_STATUS_RET(*status,NULL); */

  return param;
}

void init_par_relxilllp_dens_nthcomp(relParam **rel_param,
                                     xillParam **xill_param,
                                     const double *inp_par,
                                     const int n_parameter,
                                     int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILLLP, EMIS_TYPE_LP, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLLPDENS_NTHCOMP, PRIM_SPEC_NTHCOMP, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILLLPDENS_NTHCOMP);

  param->a = inp_par[0];
  param->incl = inp_par[1] * M_PI / 180;
  param->rin = inp_par[2];
  param->rout = inp_par[3];
  param->height = inp_par[4];
  param->htop = inp_par[5];
  param->beta = inp_par[6];

  param->gamma = inp_par[7];
  xparam->gam = param->gamma;
  xparam->lxi = inp_par[8];
  xparam->afe = inp_par[9];
  xparam->ect = inp_par[10];
  xparam->dens = inp_par[11];

  xparam->refl_frac = inp_par[12];
  xparam->fixReflFrac = (int) (inp_par[13] + 0.5); // make sure there is no problem with integer conversion

  param->z = inp_par[14];
  xparam->z = param->z;

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

}

relParam *init_par_relline(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELLINE, EMIS_TYPE_BKN, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_RELLINE);

  param->lineE = inp_par[0];
  param->emis1 = inp_par[1];
  param->emis2 = inp_par[2];
  param->rbr = inp_par[3];
  param->a = inp_par[4];
  param->incl = inp_par[5] * M_PI / 180;
  param->rin = inp_par[6];
  param->rout = inp_par[7];
  param->z = inp_par[8];
  param->limb = (int) (inp_par[9] + 0.5);

  check_parameter_bounds(param, status);
  CHECK_STATUS_RET(*status, NULL);

  return param;
}

xillParam *get_std_param_xillver(int *status) {
  int n_param = NUM_PARAM_XILLVER;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
  set_std_param_xillver(inp_par);
  return init_par_xillver(inp_par, n_param, status);
}

xillParam *get_std_param_xillver_dens(int *status) {
  int n_param = NUM_PARAM_XILLVERDENS;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
  set_std_param_xillver_dens(inp_par);
  return init_par_xillver_dens(inp_par, n_param, status);
}

xillParam *get_std_param_xillver_co(int *status) {
  int n_param = NUM_PARAM_XILLVERCO;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
  set_std_param_xillverco(inp_par);
  return init_par_xillver_co(inp_par, n_param, status);
}

xillParam *get_std_param_xillver_ns(int *status) {
  int n_param = NUM_PARAM_XILLVERNS;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
  set_std_param_xillverns(inp_par);
  return init_par_xillver_ns(inp_par, n_param, status);
}

xillParam *get_std_param_xillver_nthcomp(int *status) {
  int n_param = NUM_PARAM_XILLVER_NTHCOMP;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
  set_std_param_xillver_nthcomp(inp_par);
  return init_par_xillver_nthcomp(inp_par,
                                  n_param,
                                  status);  // this is bad design, but kTe and Ecut are identical in the code
}

xillParam *get_std_param_xillver_dens_nthcomp(int *status) {
  int n_param = NUM_PARAM_XILLVERDENS_NTHCOMP;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
  set_std_param_xillver_dens_nthcomp(inp_par);
  return init_par_xillver_dens_nthcomp(inp_par, n_param, status);
}

void get_std_param_relxilllpDCp(relParam **rel_param, xillParam **xill_param, int *status) {
  int n_param = NUM_PARAM_RELXILLLPDENS_NTHCOMP;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_VOID_STATUS(inp_par, status)

  set_std_param_relxilllpdens_nthcomp(inp_par);
  init_par_relxilllp_dens_nthcomp(rel_param, xill_param, inp_par, n_param, status);

}

relParam *get_std_param_relline(int *status) {
  int n_param = NUM_PARAM_RELLINE;
  double inp_par[NUM_PARAM_RELLINE];
  set_std_param_relline(inp_par);
  return init_par_relline(inp_par, n_param, status);
}


xillSpec *get_std_xill_spec(int *status) {
  xillParam *xill_param = get_std_param_xillver(status);
  xillSpec *xill_spec = get_xillver_spectra(xill_param, status);
  return xill_spec;
}

void get_RelProfileConstEmisZones(rel_spec **p_rel_profile, relParam **p_rel_param, int nzones, int *status) {

  // relline is per default re-normalized, but we need to test the physical normalization here
  putenv("RELLINE_PHYSICAL_NORM=1");

  *p_rel_param = get_std_param_relline(status);
  (*p_rel_param)->num_zones = nzones;
  // set emissivity to constant
  (*p_rel_param)->emis1 = 0.0;
  (*p_rel_param)->emis2 = 0.0;

  int n_ener;
  double *ener;
  get_std_relxill_energy_grid(&n_ener, &ener, status);

  *p_rel_profile = relbase(ener, n_ener, *p_rel_param, NULL, status);
  putenv("RELLINE_PHYSICAL_NORM=0");

}

rel_spec *get_stdRelProfile(int *status) {

  relParam *rel_param = get_std_param_relline(status);

  int n_ener;
  double *ener;
  get_std_relxill_energy_grid(&n_ener, &ener, status);

  return relbase(ener, n_ener, rel_param, NULL, status);
}

void init_std_relXill_spec(rel_spec **rel_profile, double **xill_spec_output, int *status) {

  CHECK_STATUS_VOID(*status);

  *rel_profile = get_stdRelProfile(status);

  double *xill_flux = malloc(sizeof(double) * (*rel_profile)->n_ener);
  CHECK_MALLOC_VOID_STATUS(xill_flux, status)

  xillSpec *xill_spec_table = get_std_xill_spec(status);
  rebin_spectrum((*rel_profile)->ener,
                 xill_flux,
                 (*rel_profile)->n_ener,
                 xill_spec_table->ener,
                 xill_spec_table->flu[0],
                 xill_spec_table->n_ener);
  *xill_spec_output = xill_flux;

}

double calc_FluxInStdBand(const double *flux, double *ener, const int n) {

  const double ELO_REF_RELPROFILE = 0.1;
  const double EHI_REF_RELPROFILE = 1.0;

  return calcSumInEnergyBand(flux, n, ener, ELO_REF_RELPROFILE, EHI_REF_RELPROFILE);
}

double calc_RelatFluxInStdBand(const rel_spec *spec) {

  assert(spec->n_zones == 1);
  return calc_FluxInStdBand(spec->flux[0], spec->ener, spec->n_ener);
}

double calc_XillverFluxInStdBand(const double *xillverSpec, double *ener, int n) {

  return calc_FluxInStdBand(xillverSpec, ener, n);
}

void compareReferenceFlux(double flux, double refFlux, int *status) {

  const double PREC = 1e-6;

  CHECK_STATUS_VOID(*status);

  if (fabs(flux - refFlux) > PREC) {
    RELXILL_ERROR("failed comparing the calculated flux to the reference flux", status);
    printf("  expecting a flux of %e, but calculated %e \n", refFlux, flux);
  }

}

int test_stdEvaluationRelatFlux(const rel_spec *rel_profile) {

  PRINT_RELXILL_TEST_MSG_DEFAULT();

  int status = EXIT_SUCCESS;

  const double ReferenceRelatStdFlux = 8.371512e-01;
  double relatFlux = calc_RelatFluxInStdBand(rel_profile);
  compareReferenceFlux(relatFlux, ReferenceRelatStdFlux, &status);

  print_relxill_test_result(status);
  return status;
}

int test_stdEvaluationXillverFlux(const double *xill_spec, const rel_spec *rel_profile) {

  PRINT_RELXILL_TEST_MSG_DEFAULT();

  int status = EXIT_SUCCESS;

  const double ReferenceXillverStdFlux = 1.802954e+01;
  double xillFlux = calc_XillverFluxInStdBand(xill_spec, rel_profile->ener, rel_profile->n_ener);
  compareReferenceFlux(xillFlux, ReferenceXillverStdFlux, &status);

  print_relxill_test_result(status);
  return status;
}

void test_stdEvaluationFluxes(int *status) {

  CHECK_STATUS_VOID(*status);

  PRINT_RELXILL_TEST_MSG("  : \n");

  rel_spec *rel_profile = NULL;
  double *xill_spec = NULL;
  init_std_relXill_spec(&rel_profile, &xill_spec, status);

  int relStatus = test_stdEvaluationRelatFlux(rel_profile);

  int xillStatus = test_stdEvaluationXillverFlux(xill_spec, rel_profile);

  if ((relStatus != EXIT_SUCCESS) || (xillStatus != EXIT_SUCCESS)) {
    *status = EXIT_FAILURE;
  }
}

