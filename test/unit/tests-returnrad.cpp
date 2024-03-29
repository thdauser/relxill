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

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "XspecSpectrum.h"
#include "common-functions.h"
#include "Rellp.h"
#include "IonGradient.h"
#include "Xillspec.h"
#include "Relphysics.h"
#include "Relreturn_Datastruct.h"
#include "Relreturn_Table.h"

extern "C" {
#include "relutility.h"
#include "writeOutfiles.h"
#include "xilltable.h"
}

#include <vector>

#define PREC 1e-6


static int is_grid_consistent(const double* rlo, const double* rhi, int nrad){
  for (int ii=1; ii<nrad; ii++){
    if (fabs(rhi[ii-1]-rlo[ii])>PREC){
      printf(" *** error: lower boundary rlo[ii]-rhi[ii-1] = %e - %e = %e  for ii=%i \n",
             rlo[ii], rhi[ii-1],
             fabs(rhi[ii-1]-rlo[ii]), ii);
      return 0;
    }
  }
  return 1;
}

static void test_table_radial_grid(tabulatedReturnFractions *dat) {

  REQUIRE(dat->nrad == RETURNRAD_TABLE_NR);
  REQUIRE(dat->ng == RETURNRAD_TABLE_NG);

  REQUIRE(dat->rlo[0] > 1.0);
  REQUIRE(dat->rhi[dat->nrad - 1] <= 1000.1); // don't care about border effects so add the 0.1

  REQUIRE(is_grid_consistent(dat->rlo, dat->rhi, dat->nrad)==1);
}

static int test_table_fracE_norm(const tabulatedReturnFractions *dat) {

  int status = EXIT_SUCCESS;
  double kSumfrac;
  double kSumfrac_ref = 1.0;

  for (int jj=0; jj<dat->nrad; jj++) {

    kSumfrac=0.0;
    for (int ii = 0; ii < dat->nrad; ii++) {

      kSumfrac += dat->frac_e[ii][jj];
    }
    if (fabs(kSumfrac - kSumfrac_ref) > PREC) {
      RELXILL_ERROR("testing the normalization of FRAC_E failed (return radiation table)", &status);
      printf(" expecting a normalization of %e, but found %e\n", kSumfrac_ref, kSumfrac);
      break;
    }
  }
  return status;
}

static int test_table_fracG_norm(const tabulatedReturnFractions *dat) {

  int status = EXIT_SUCCESS;
  double kSumfrac;
  double kSumfrac_ref = 1.0;

  for (int ii=0; ii<dat->nrad; ii++) {
    for (int jj=0; jj<dat->nrad; jj++) {

      kSumfrac = 0.0;
      for (int kk = 0; kk < dat->ng; kk++) {
        kSumfrac += dat->frac_g[ii][jj][kk];
      }

      if (fabs(kSumfrac-kSumfrac_ref)>PREC){
        RELXILL_ERROR("testing the normalization of FRAC_G failed (return radiation table)",&status);
        printf (" [%02i][%02i] expecting a normalization of %e, but found %e\n", ii,jj,kSumfrac_ref, kSumfrac);
        break;
      }
    }
  }
  return status;
}


static void test_table_fractions_normalization(tabulatedReturnFractions* dat){
  REQUIRE(dat->frac_e[0][0]!= 0);
//  REQUIRE(test_table_fracE_norm(dat) == EXIT_SUCCESS);
  REQUIRE(test_table_fracG_norm(dat) == EXIT_SUCCESS);
}



// ------- //
TEST_CASE(" Returning Radiation Table ", "[table]") {

  int status = EXIT_SUCCESS;
  returnTable *tab = get_returnrad_table(&status);
  REQUIRE(tab != NULL);
  REQUIRE(status == EXIT_SUCCESS);

  DYNAMIC_SECTION("testing radial grid") {
    for (int ii = 0; ii < tab->nspin; ii++) {
      test_table_radial_grid(tab->retFrac[ii]);
      REQUIRE(status == EXIT_SUCCESS);
    }
  }


  DYNAMIC_SECTION("testing fractions normalization") {
    for (int ii = 0; ii < tab->nspin; ii++) {
      test_table_fractions_normalization(tab->retFrac[ii]);
    }
  }

}





static void invert_emis_profile(emisProfile* emis){
  invertArray(emis->re, emis->nr);
  invertArray(emis->emis, emis->nr);

  if(emis->del_inc!= nullptr){
    invertArray(emis->del_inc, emis->nr);
  }

  if(emis->del_emit!= nullptr){
    invertArray(emis->del_emit, emis->nr);
  }
}

static void write_emis_profile(const std::string& fname, emisProfile* emis_profile){
  write_data_to_file( fname.c_str() , emis_profile->re, emis_profile->emis, emis_profile->nr);
}


emisProfile* get_test_emis_rrad(double rin, double rout, double spin, int* status){

  returningFractions *rf = get_rrad_fractions(spin,rin , rout, status);
  double gamma = 2;

  emisProfile* emis = new_emisProfile(rf->rad, rf->nrad, status);
  get_emis_bkn(emis->emis, emis->re, emis->nr,3.0,3.0,emis->re[0]);
  for (int ii=0; ii<rf->nrad; ii++){
    emis->emis[ii] /= emis->emis[0]; //
  }

  emisProfile* emis_return = calc_rrad_emis_corona(rf, nullptr, emis, gamma, status);

  free_returningFractions(&rf);
  free_emisProfile(emis);

  return emis_return;
}





// ------- //
TEST_CASE(" Changing number of radial bins if Rin is increased", "[returnrad]") {

  int status = EXIT_SUCCESS;

  double spin = 0.9;
  double Rout = 1000;

  double Rin = kerr_rms(spin);
  returningFractions *dat = get_rrad_fractions(spin, Rin, Rout, &status);
  REQUIRE(status==EXIT_SUCCESS);
  int nrad_rms = dat->nrad;

  free_returningFractions(&dat);

  Rin *= 2;
  dat = get_rrad_fractions(spin, Rin, Rout, &status);
  REQUIRE(status==EXIT_SUCCESS);
  int nrad_rfac2 = dat->nrad;

  free_returningFractions(&dat);

  REQUIRE(nrad_rms > nrad_rfac2);

}


// ------- //
TEST_CASE(" Rebining the return rad emissivity profile", "[returnrad]") {

  int status = EXIT_SUCCESS;

  double precRebinCoarse = 0.05;

  LocalModel lmod(ModelName::relline_lp);
  lmod.set_par(XPar::switch_switch_returnrad, 0);  // need this for the comparison (which is without return_rad)
  relParam* rel_param = lmod.get_rel_params();

  RelSysPar *sysPar = get_system_parameters(rel_param, &status);

 returningFractions *dat = get_rrad_fractions(rel_param->a, rel_param->rin, rel_param->rout, &status);
 invertArray(dat->rad, dat->nrad);
 emisProfile *emisCoarse = calc_emis_profile(dat->rad, dat->nrad, rel_param, &status);
 invert_emis_profile(emisCoarse);

  emisProfile *emisFineReference = calc_emis_profile(sysPar->re, sysPar->nr, rel_param, &status);
  emisProfile *emisRebin = new_emisProfile(sysPar->re, sysPar->nr, &status);
  rebin_emisprofile_on_radial_grid(emisRebin, emisCoarse, &status);

  REQUIRE(status==EXIT_SUCCESS);

  write_emis_profile("test_emisCoarse.dat", emisCoarse);
  write_emis_profile("test_emisFineReference.dat", emisFineReference);
  write_emis_profile("test_emisRebin.dat", emisRebin);

  for (int ii = 20; ii < emisRebin->nr;  ii++) {  // last bin in coarse grid deviates, so skip, as it is not the interpolation
    if (is_debug_run()) {
      printf(" rad: %.3e :  %e (ref=%e, ratio=%e) \n",
             emisRebin->re[ii], emisRebin->emis[ii], emisFineReference->emis[ii],
             emisRebin->emis[ii] / emisFineReference->emis[ii]);
    }
    REQUIRE(fabs(emisRebin->emis[ii] / emisFineReference->emis[ii] - 1) < precRebinCoarse);
  }

}

// ------- //
TEST_CASE(" Line profile for Returning Radiation ", "[returnrad]") {


  int status = EXIT_SUCCESS;

  LocalModel lmod(ModelName::relline_lp);
  relParam* rel_param = lmod.get_rel_params();
  rel_param->a = 0.998;
  rel_param->height = 5.0;
  rel_param->rout = 1000.0;

  const char* env_name_outfiles = "RELXILL_WRITE_OUTFILES";

  TestSpectrum *spec = getNewSpec(0.05, 10, 1000, &status);

  rel_param->return_rad=0;
  relline_spec_multizone *rel_profile_norrad = relbase(spec->ener, spec->nbins, rel_param, &status);

  setenv(env_name_outfiles, "1", 1);
  rel_param->return_rad = 1;
  relline_spec_multizone *rel_profile = relbase(spec->ener, spec->nbins, rel_param, &status);
  setenv(env_name_outfiles, "0", 1);

  REQUIRE(status==EXIT_SUCCESS);
  REQUIRE(rel_profile->n_zones == 1);
  REQUIRE(calcSum(rel_profile->flux[0], rel_profile->n_ener) > 1e-8);

  double abs_diff_profiles[rel_profile->n_ener];
  for (int ii = 0; ii < rel_profile->n_ener; ii++) {
    abs_diff_profiles[ii] = abs(rel_profile_norrad->flux[0][ii] - rel_profile->flux[0][ii]);
  }

  INFO(" require that the line profile with and without return radiation differs ");
  CHECK_THAT(calcSum(abs_diff_profiles, rel_profile->n_ener),
             !Catch::Matchers::WithinAbs(0.0, 1e-2)
  );

  free_Spectrum(&spec);

}


// ------- //
TEST_CASE(" Write emissivity profile", "[returnrad-emis]") {

  int status = EXIT_SUCCESS;

  LocalModel lmod(ModelName::relline_lp);
  lmod.set_par(XPar::switch_switch_returnrad, -1);  // need this for the comparison (which is without return_rad)
  relParam *rel_param = lmod.get_rel_params();
  RelSysPar *sysPar = get_system_parameters(rel_param, &status);
  emisProfile *emis_profile = calc_emis_profile(sysPar->re, sysPar->nr, rel_param, &status);
  write_emis_profile("__output_emis_profile_rrad.dat", emis_profile);

  REQUIRE(status == EXIT_SUCCESS);

}





// ------- //
TEST_CASE(" Increasing Rin has to reduce the flux at the next zone (in radius)", "[returnrad]") {

  int status = EXIT_SUCCESS;
  double spin = 0.99;
  double rin = kerr_rms(spin);
  double rin_grid = rin * 1.001;
  double rout = 1000;

  returningFractions *rf0 = get_rrad_fractions(spin, rin, rout, &status);
  returningFractions *rf_grid = get_rrad_fractions(spin, rin_grid, rout, &status);

  REQUIRE(rf0 != NULL);
  REQUIRE(rf_grid != NULL);

  INFO("changing Rin slightly has to reduce the returning incident fraction onto the next zone ");
  REQUIRE(abs(rf_grid->tf_r[1][0] - rf0->tf_r[1][0]) > 1e-8); // tf_r[ind_ri][ind_re]

  INFO("but not change for any other emitting radius");
  REQUIRE(abs(rf_grid->tf_r[1][1] - rf0->tf_r[1][1]) < 1e-8);

}

/*
// ------- //
TEST_CASE(" Interpolation of Returning Radiation Fractions", "[returnrad]") {

  int status = EXIT_SUCCESS;
  double spin = 0.995;
  int ind_rad_test = 2;
  double rout = 1000;

  returningFractions *rf = get_rrad_fractions(spin, kerr_rms(spin) , rout, &status);
  double rad_grid = rf->rad[ind_rad_test];

  double rad_lo = rad_grid*0.99999;
  double rad_hi  = rad_grid*1.00001;

  emisProfile *emis_lo = get_test_emis_rrad(rad_lo, rout, spin, &status);
  emisProfile *emis_0 = get_test_emis_rrad(rad_grid, rout, spin, &status);
  emisProfile *emis_hi = get_test_emis_rrad(rad_hi, rout, spin, &status);

  double sum_lo = calcSum(emis_lo->emis, emis_lo->nr);
  double sum_0 = calcSum(emis_0->emis, emis_0->nr);
  double sum_hi = calcSum(emis_hi->emis, emis_hi->nr);

  REQUIRE(sum_0 > 0);
  REQUIRE(sum_lo / sum_0 > 1.0001);
  REQUIRE(sum_0 / sum_hi > 1.0001);

  // absolute difference for a tiny change around the grid-point less than 1%
  REQUIRE(sum_lo - sum_0 < 0.01);
  REQUIRE(sum_0 - sum_hi < 0.01);

  // absolute difference for a tiny change around the grid-point less than 1%
  REQUIRE(sum_lo / sum_0 < 1.01);
  REQUIRE(sum_0 / sum_hi < 1.01);

  // fractional change in both direction of the grid point should be fairly similar
  REQUIRE( abs( sum_lo/sum_0 - sum_0/sum_hi) < 1e-4 );

} */


// ------- //
TEST_CASE(" Changing Rin should result in a change in line shape", "[returnrad-single]") {

  double spin = 0.995;
  double rin_scaled = kerr_rms(0.99025);
  double rin_scaled2 = kerr_rms(0.990);

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  const char* env_outfiles = "RELXILL_WRITE_OUTFILES";
  setenv(env_outfiles, "1", 1);

  LocalModel local_model{ModelName::relxilllp};
  local_model.set_par(XPar::a, spin);
  local_model.set_par(XPar::switch_switch_returnrad, 1);
  local_model.set_par(XPar::rout,1000.0);


  local_model.set_par(XPar::rin, rin_scaled2);
  local_model.eval_model(spec);

  double sum1 = sum_flux(spec.flux, spec.num_flux_bins());

  local_model.set_par(XPar::rin, rin_scaled);
  local_model.eval_model(spec);
  double sum2 = sum_flux(spec.flux, spec.num_flux_bins());

  REQUIRE(sum1 > 1e-6);
  REQUIRE(abs(sum1 - sum2) > 1e-3);

}

TEST_CASE("Test Flux Correction Factor", "[returnrad]") {
  //

  int status = EXIT_SUCCESS;

  LocalModel lmod{ModelName::relxilllp};
  lmod.set_par(XPar::switch_switch_returnrad, 1);

  lmod.set_par(XPar::logxi, 3.0);
  double flux_corr1=0;
  xillSpec* xill_spec = get_xillver_spectra(lmod.get_xill_params(), &status);
  get_xillver_fluxcorrection_factors(xill_spec, &flux_corr1, nullptr,
                                     get_xilltab_param(lmod.get_xill_params(), &status),
                                     &status);
  REQUIRE(flux_corr1 > 0.5);
  REQUIRE(flux_corr1 < 1.0); // only works for logxi<=4

  lmod.set_par(XPar::logxi, 0.0);
  xillSpec* xill_spec2 = get_xillver_spectra(lmod.get_xill_params(), &status);
  double flux_corr2=0;
  get_xillver_fluxcorrection_factors(xill_spec2, &flux_corr2, nullptr,
                                     get_xilltab_param(lmod.get_xill_params(), &status) ,
                                     &status);

  REQUIRE(flux_corr1 > 2 * flux_corr2);  // for low instead high ionization it is much lower, only works if logxi_hi>3.5

  free_xill_spec(xill_spec);
  free_xill_spec(xill_spec2);

}


TEST_CASE("Test Gshift Correction Factor", "[returnrad]") {
  int status = EXIT_SUCCESS;

  LocalModel lmod{ModelName::relxilllp};
  lmod.set_par(XPar::switch_switch_returnrad, 1);

  double gshift_ref_value = 1.5; // correction factor calculated for this value
  double gam = 2.0;

  lmod.set_par(XPar::gamma, gam);
  lmod.set_par(XPar::logxi, 3.0);
  double gshift_flux_corr1;
  xillSpec* xill_spec = get_xillver_spectra(lmod.get_xill_params(), &status);
  get_xillver_fluxcorrection_factors(xill_spec, nullptr, &gshift_flux_corr1,
                                     get_xilltab_param(lmod.get_xill_params(), &status), &status);
  REQUIRE(gshift_flux_corr1 > 1);
  REQUIRE(gshift_flux_corr1 < pow(gshift_ref_value,gam+0.5)/pow(gshift_ref_value,gam));

  lmod.set_par(XPar::logxi, 0.0);
  double gshift_flux_corr2;
  xillSpec* xill_spec2 = get_xillver_spectra(lmod.get_xill_params(), &status);
  get_xillver_fluxcorrection_factors(xill_spec2, nullptr, &gshift_flux_corr2,
                                     get_xilltab_param(lmod.get_xill_params(), &status), &status);

  REQUIRE(gshift_flux_corr2 < 1 );
  REQUIRE(gshift_flux_corr1 > 1.5 * gshift_flux_corr2);  // for low instead high ionization it is much lower, only works if logxi_hi>3.5

  free_xill_spec(xill_spec);
  free_xill_spec(xill_spec2);

}


TEST_CASE("Test Gshift Linear Interpolation of the Correction Factor", "[returnrad]") {

  double gamma = 2.0;

  double xill_gshift_fac = 1.2;
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 1.5, gamma) > pow(1.5, gamma));
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) > pow(0.3, gamma));
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 1.5, gamma) > 1);
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) < 1);
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) > 0);

  xill_gshift_fac = 0.9;
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 1.5, gamma) < pow(1.5, gamma));
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) < pow(0.3, gamma));
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 1.5, gamma) > 1);
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) < 1);
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) > 0);

}

TEST_CASE("Write Flux Correction Factor", "[returnrad]") {
  //

  int status = EXIT_SUCCESS;

  LocalModel lmod{ModelName::relxilllp};
  lmod.set_par(XPar::switch_switch_returnrad, 1);

  const int n_logxi = 100;
  double logxi[n_logxi] = {0};
  double flux_corr[n_logxi] = {0};
  double gshift_corr[n_logxi] = {0};

  double logxi_min = 0.0;
  double logxi_max = 4.0;

  const int n_gam = 3;
  double values_gamma[n_gam] = {1.5, 2.0, 2.5};

  char buffer[100];
  for (double gam : values_gamma) {
    lmod.set_par(XPar::gamma, gam);
    for (int ii = 0; ii < n_logxi; ii++) {
      logxi[ii] = (logxi_max - logxi_min) * (ii * 1.0 / (n_logxi - 1)) + logxi_min;
      lmod.set_par(XPar::logxi, logxi[ii]);
      xillSpec* xill_spec = get_xillver_spectra(lmod.get_xill_params(), &status);
      get_xillver_fluxcorrection_factors(xill_spec, &(flux_corr[ii]), &(gshift_corr[ii]),
                                         get_xilltab_param(lmod.get_xill_params(), &status), &status);
      free_xill_spec(xill_spec);
    }
    int retval = sprintf(buffer, "test-flux-corr-returnrad_gam%.1f.dat",gam);
    write_data_to_file(buffer, logxi, flux_corr, n_logxi);
    retval = sprintf(buffer, "test-gshift-corr-returnrad_gam%.1f.dat",gam);
    write_data_to_file(buffer, logxi, gshift_corr, n_logxi);
  }

}


// ------- //
TEST_CASE(" Test relxilllp including returning radiation", "[returnrad]") {

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllp};

  local_model.set_par(XPar::switch_switch_returnrad, 0);
  local_model.eval_model(spec);
  double model_no_rrad_flux = sum_flux(spec.flux, spec.num_flux_bins());


  local_model.set_par(XPar::switch_switch_returnrad, 1);
  local_model.eval_model(spec);
  double model_rrad_flux = sum_flux(spec.flux, spec.num_flux_bins());

  string fname = "test-spec-relxilllpret.dat";
  save_xillver_spectrum(spec.energy, spec.flux, spec.num_flux_bins(), fname.data());

  REQUIRE(model_rrad_flux > model_no_rrad_flux);


}



// ------- //
TEST_CASE(" Test return rad ENV variable", "[returnrad]") {

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllp};

  local_model.set_par(XPar::switch_switch_returnrad, 0);
  local_model.eval_model(spec);
  double model_no_rrad_flux = sum_flux(spec.flux, spec.num_flux_bins());


  const char* env = "RELXILL_RETURNRAD_SWITCH";
  setenv(env, "1", 1);
  local_model.eval_model(spec);
  double model_rrad_flux = sum_flux(spec.flux, spec.num_flux_bins());

  unsetenv(env);

  // env variable can not change if parameter value is given
  REQUIRE( fabs(model_rrad_flux - model_no_rrad_flux) <1e-8);
}

static void require_all_values_the_same(const double* val, const int n){
  double avg = calcSum(val, n) / n;
  REQUIRE( fabs( avg - val[0]) < 1e-6);
}

static void require_values_differ(const double* val, const int n){
  double avg = calcSum(val, n) / n;
  REQUIRE( fabs( (avg - val[0]) / avg) > 1e-2);
  REQUIRE( fabs( (avg - val[n-1]) / avg) > 1e-2);
}

TEST_CASE("Calculation of Correction Factors", "[returnrad]"){

  int status = EXIT_SUCCESS;

  LocalModel lmod{ModelName::relxilllp};
  lmod.set_par(XPar::switch_switch_returnrad, 1);

  relParam* rel_param = lmod.get_rel_params();

  const auto n_zones = rel_param->num_zones;
  double logxi_min = 0.0;
  double logxi_max = 4.0;

  auto radial_grid = RadialGrid(rel_param->rin, rel_param->rout, n_zones, rel_param->height);

  xillTableParam* xill_table_param[n_zones];
  xillSpec* xill_spec[n_zones];

  // TEST contant spectrum over the disk
  xillParam* xill_param = lmod.get_xill_params();
  for (int ii = 0; ii < n_zones; ii++) {
    xill_table_param[ii] = get_xilltab_param(xill_param, &status);
    xill_spec[ii] = get_xillver_spectra_table(xill_table_param[ii], &status);
  }

  rel_param->rrad_corr_factors =
      calc_rrad_corr_factors(xill_spec, radial_grid, xill_table_param, &status);

  require_all_values_the_same(rel_param->rrad_corr_factors->corrfac_flux, n_zones);
  require_all_values_the_same(rel_param->rrad_corr_factors->corrfac_gshift, n_zones);

  returningFractions *ret_fractions = get_rrad_fractions(rel_param->a, rel_param->rin, rel_param->rout, &status);

  rradCorrFactors* tab_corrfactors =
      rebin_corrfactors_to_rradtable_grid(rel_param->rrad_corr_factors, ret_fractions, &status);

  require_all_values_the_same(tab_corrfactors->corrfac_flux, tab_corrfactors->n_zones);
  require_all_values_the_same(tab_corrfactors->corrfac_gshift, tab_corrfactors->n_zones);

  // as all values are the same, re-binning should not change the actual value
  REQUIRE(rel_param->rrad_corr_factors->corrfac_flux[0] == tab_corrfactors->corrfac_flux[0]);
  REQUIRE(rel_param->rrad_corr_factors->corrfac_gshift[0] == tab_corrfactors->corrfac_gshift[0]);


  // TEST ionization gradient
  xillSpec* xill_spec_grad[n_zones];
  for (int ii = 0; ii < n_zones; ii++) {
    xill_table_param[ii]->lxi = (logxi_max - logxi_min) * (ii * 1.0 / (n_zones - 1)) + logxi_min;
    xill_spec_grad[ii] = get_xillver_spectra_table(xill_table_param[ii], &status);
  }

  rradCorrFactors* grad_corrfactors =
      calc_rrad_corr_factors(xill_spec_grad, radial_grid, xill_table_param, &status);

  require_values_differ(grad_corrfactors->corrfac_flux, n_zones);
  require_values_differ(grad_corrfactors->corrfac_gshift, n_zones);

  rradCorrFactors* tab_grad_corrfactors =
      rebin_corrfactors_to_rradtable_grid(grad_corrfactors, ret_fractions, &status);
  require_values_differ(tab_grad_corrfactors->corrfac_flux, tab_grad_corrfactors->n_zones);
  require_values_differ(tab_grad_corrfactors->corrfac_gshift, tab_grad_corrfactors->n_zones);

}

TEST_CASE(" Test that no warning is printed", "[returnrad]") {

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllpCp};

  local_model.set_par(XPar::a, 0.8);
  local_model.set_par(XPar::switch_switch_reflfrac_boost, 1);
  local_model.set_par(XPar::refl_frac, 1);
  local_model.set_par(XPar::afe, 10.0);
  local_model.set_par(XPar::h, 6.);
  local_model.set_par(XPar::switch_iongrad_type,2);
  local_model.set_par(XPar::switch_switch_returnrad, 1);

//  relParam* rel_param = local_model.get_rel_params();
//  xillParam * xill_param = local_model.get_xill_params();

  local_model.eval_model(spec);
  double model1 = sum_flux(spec.flux, spec.num_flux_bins());

}

TEST_CASE(" Caching of spectrum with return radiation correction factors", "[returnrad]"){

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllp};

  double a_default = 0.998;
  local_model.set_par(XPar::a, a_default);
  local_model.set_par(XPar::h, 3.);
  local_model.set_par(XPar::switch_switch_reflfrac_boost, 1);
  local_model.set_par(XPar::refl_frac, 1);


  local_model.eval_model(spec);
  double model1 = sum_flux(spec.flux, spec.num_flux_bins());

  // reset relxill total cache by this
  local_model.set_par(XPar::a, a_default*0.99);
  local_model.eval_model(spec);

  // evaluate wih the same values again, should give the same result
  local_model.set_par(XPar::a, a_default);
  local_model.eval_model(spec);
  double model2 = sum_flux(spec.flux, spec.num_flux_bins());
  REQUIRE( fabs(model1 - model2)/model1 < 1e-6);


}


TEST_CASE(" Test that returnrad is also implemented for low spin","[returnrad2]"){

  LocalModel local_model{ModelName::relxilllpCp};

  local_model.set_par(XPar::a, 0.0 );
  local_model.set_par(XPar::switch_switch_reflfrac_boost, 1);
  local_model.set_par(XPar::refl_frac, -1);
  local_model.set_par(XPar::afe, 1.0);
  local_model.set_par(XPar::h, 3.);
  local_model.set_par(XPar::rin, 10);
  local_model.set_par(XPar::switch_iongrad_type,0);

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  local_model.set_par(XPar::switch_switch_returnrad, 1);
  local_model.eval_model(spec);
  double model1 = sum_flux(spec.flux, spec.num_flux_bins());

  local_model.set_par(XPar::switch_switch_returnrad, 0);
  local_model.eval_model(spec);
  double model2 = sum_flux(spec.flux, spec.num_flux_bins());

  REQUIRE( fabs(model1 - model2)/model1 > 1e-3 );


}