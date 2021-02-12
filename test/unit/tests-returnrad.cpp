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
//#include "LocalModel.h"
//#include "XspecSpectrum.h"
//#include "common-functions.h"
extern "C" {
#include "relreturn.h"
#include "relutility.h"
#include "test_relxill.h"
#include "writeOutfiles.h"
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

static void test_table_radial_grid(returnFracData *dat) {

  REQUIRE(dat->nrad == RETURNRAD_TABLE_NR);
  REQUIRE(dat->ng == RETURNRAD_TABLE_NG);

  REQUIRE(dat->rlo[0] > 1.0);
  REQUIRE(dat->rhi[dat->nrad - 1] <= 1000.1); // don't care about border effects so add the 0.1

  REQUIRE(is_grid_consistent(dat->rlo, dat->rhi, dat->nrad)==1);
}

static int test_table_fracE_norm(const returnFracData *dat) {

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

static int test_table_fracG_norm(const returnFracData *dat) {

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


static void test_table_fractions_normalization(returnFracData* dat){
  REQUIRE(dat->frac_e[0][0]!= 0);
  REQUIRE(test_table_fracE_norm(dat) == EXIT_SUCCESS);
  REQUIRE(test_table_fracG_norm(dat) == EXIT_SUCCESS);
}



// ------- //
TEST_CASE(" Returning Radiation Table ", "[table]") {

  int status = EXIT_SUCCESS;
  returnTable *tab = get_returnRadTable(&status);
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




// ------- //
TEST_CASE(" Changing number of radial bins if Rin is increased", "[returnrad]") {

  int status = EXIT_SUCCESS;

  double spin = 0.9;
  double Rout = 1000;

  double Rin = kerr_rms(spin);
  returningFractions *dat = get_rr_fractions(spin, Rin, Rout, &status);
  REQUIRE(status==EXIT_SUCCESS);
  int nrad_rms = dat->nrad;

  Rin *= 2;
  dat = get_rr_fractions(spin, Rin, Rout, &status);
  REQUIRE(status==EXIT_SUCCESS);
  int nrad_rfac2 = dat->nrad;

  REQUIRE(nrad_rms > nrad_rfac2);

}



// ------- //
TEST_CASE(" Return Fraction interpolation for different spins", "[returnrad]") {

  int status = EXIT_SUCCESS;

  int nspin = 2;
  double spin[] = { 0.86,0.861 };

  double Rout = 1000;

  INFO(" this test is missing any useful actions \n");
  for (int ii=0; ii<nspin; ii++){

    double Rin = kerr_rms(spin[ii]);
    returningFractions *dat = get_rr_fractions(spin[ii], Rin, Rout, &status);

    REQUIRE(status==EXIT_SUCCESS);
    REQUIRE(dat != nullptr);

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

// ------- //
TEST_CASE(" Rebining the return rad emissivity profile", "[returnrad]") {

  int status = EXIT_SUCCESS;

  double precRebinCoarse = 0.05;

  relParam* rel_param =  get_std_param_rellinelp(&status);

  RelSysPar *sysPar = get_system_parameters(rel_param, &status);

 returningFractions *dat = get_rr_fractions(rel_param->a, rel_param->rin, rel_param->rout, &status);
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

  relParam *rel_param = get_std_param_rellinelp(&status);
  rel_param->a = 0.998;
  rel_param->height = 5.0;
  rel_param->rout = 1000.0;

  const char* env_name_outfiles = "RELXILL_WRITE_OUTFILES";

  Spectrum *spec = getNewSpec(0.05, 10, 1000, &status);

  rel_param->return_rad=0;
  rel_spec *rel_profile_norrad = relbase(spec->ener, spec->nbins, rel_param, nullptr, &status);

  setenv(env_name_outfiles, "1", 1);
  rel_param->return_rad = 1;
  rel_spec *rel_profile = relbase(spec->ener, spec->nbins, rel_param, nullptr, &status);
  setenv(env_name_outfiles, "0", 1);

  REQUIRE(status==EXIT_SUCCESS);
  REQUIRE(rel_profile->n_zones == 1);
  REQUIRE(calcSum(rel_profile->flux[0], rel_profile->n_ener) > 1e-8);


  double abs_diff_profiles[rel_profile->n_ener];
  for(int ii=0; ii<rel_profile->n_ener; ii++){
    abs_diff_profiles[ii] = abs(rel_profile_norrad->flux[0][ii] - rel_profile->flux[0][ii] ) ;
  }

  INFO(" require that the line profile with and without return radiation differs ");
  CHECK_THAT(calcSum(abs_diff_profiles, rel_profile->n_ener),
              ! Catch::Matchers::WithinAbs(0.0, 1e-2)
              );


}
