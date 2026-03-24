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

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "XspecSpectrum.h"
#include "Relreturn_BlackBody.h"
#include "tests-bbody-returnrad.h"

extern "C" {
#include "relutility.h"
}

#if !defined(RELXILL_SOURCE_DIR)
#warning: RELXILL_SOURCE_DIR is not defined
#endif

#define PREC 1e-6
#define PREC_REFDATA 1e-3
#define REFDATA_DIR "test/unit/refdata/"  // relativ to RELXILL_SOURCE_DIR (set when compiling)



// Function definitions
static double* integSpecArea(double** spec, int n_ener, double spin, const returnSpec2D *returnSpec, int* status);

static void fits_write_spec(const char* fname, const double* ener, double* spec, int n_ener, int* status){

  CHECK_STATUS_VOID(*status);

  char filename[300];
  if (fname[0]!='!'){
    sprintf(filename,"!%s", fname);
  } else {
    sprintf(filename,"%s",fname);
  }

  // open the fits file
  fitsfile *fptr;
  if (fits_create_file(&fptr, filename, status)) {
    CHECK_STATUS_VOID(*status);
    printf("  *** error : creating file %s failed\n", fname);
    return;
  }

  char tf1[3] = "1D";
  char f1[10] = "bin_lo";
  char f2[10] = "bin_hi";
  char f3[10] = "flux";

  const int tfields = 3;
  char *ttype[] = { f1, f2, f3 };
  char *tform[] = {tf1, tf1, tf1};

  if (fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, nullptr, "spec", status)) {
    CHECK_STATUS_VOID(*status);
  }

  double bin_lo[n_ener];
  double bin_hi[n_ener];

  for (int ii=0; ii<n_ener; ii++){
    bin_lo[ii] = ener[ii];
    bin_hi[ii] = ener[ii+1];
  }

  int firstrow = 1;  /* first row in table to write   */
  int firstelem = 1;  /* first element in row  */

  fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, n_ener, bin_lo, status);
  fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, n_ener, bin_hi, status);
  fits_write_col(fptr, TDOUBLE, 3, firstrow, firstelem, n_ener, spec, status);

  relxill_check_fits_error(status);

  if (fptr != nullptr) { fits_close_file(fptr, status); }



}



static refSpecData *new_refSpecData(int nener, int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  auto *spec = new refSpecData;
  CHECK_MALLOC_RET_STATUS(spec, status, nullptr);

  spec->nener = nener;

  spec->elo = new double[nener];
  spec->ehi = new double[nener];
  spec->flux = new double[nener];

  return spec;
}

static void free_refSpecData(refSpecData *spec) {
  delete[] spec->elo;
  delete[] spec->ehi;
  delete[] spec->flux;
  delete spec;
}

static refSpecData *fits_read_refdata_spec(char *fname, int *status) {

  // open the reference data
  fitsfile *fptr = nullptr;
  (fits_open_table(&fptr, fname, READONLY, status));
  relxill_check_fits_error(status);
  CHECK_STATUS_RET(*status, nullptr);

  long nener;

  // get the number of rows
  if (fits_get_num_rows(fptr, &nener, status)) return nullptr;

  auto spec = new_refSpecData((int) nener, status);
  CHECK_MALLOC_RET_STATUS(spec, status, nullptr);

  int anynul = 0;
  double nullptrval = 0.0;
  auto nelem = (LONGLONG) nener;
  fits_read_col(fptr, TDOUBLE, 1, 1, 1, nelem, &nullptrval, spec->elo, &anynul, status);
  fits_read_col(fptr, TDOUBLE, 2, 1, 1, nelem, &nullptrval, spec->ehi, &anynul, status);
  fits_read_col(fptr, TDOUBLE, 3, 1, 1, nelem, &nullptrval, spec->flux, &anynul, status);

  if (*status != EXIT_SUCCESS) {
    printf(" *** error *** initializing of the RETURN RADIATION table %s failed \n", fname);
    free_refSpecData(spec);
  }

  if (fptr != nullptr) {
    fits_close_file(fptr, status);
  }

  return spec;
}

static int compare_refSpecDat(const refSpecData *spec, const double *ener0, const double *flux0, int nener0) {

  int status = EXIT_SUCCESS;

  if (spec->nener != nener0) {
    status = EXIT_FAILURE;
    printf(" *** error : spectrum has %i bins, but expected to have %i bins\n", spec->nener, nener0);
    return EXIT_FAILURE;
  }

  for (int ii = 0; ii < nener0; ii++) {

    if ((spec->flux[ii] / flux0[ii]) > PREC_REFDATA && (fabs(spec->flux[ii] - flux0[ii]) > PREC_REFDATA)) {
      status = EXIT_FAILURE;
      printf(" *** error : spectrum flux in bin [%.3f, %.3f] is %e, but expected to be %e \n",
             ener0[ii], ener0[ii + 1], flux0[ii], spec->flux[ii]);
    }

  }

  return status;
}

static int compareWithRefdata(const char *fname_ref, const double *ener0, const double *flux0, int nener0, int *status) {

  PRINT_RELXILL_TEST_MSG(fname_ref);

  char fullname_ref[1000];
  sprintf(fullname_ref, "%s/%s/%s", RELXILL_SOURCE_DIR, REFDATA_DIR, fname_ref);
  refSpecData *spec = fits_read_refdata_spec(fullname_ref, status);

  if (spec == nullptr){
    return -1;
  }

  int status_comp = compare_refSpecDat(spec, ener0, flux0, nener0);
  free_refSpecData(spec);
  return status_comp;
}



static double* integSpecArea(double** spec, int n_ener, double spin, const returnSpec2D *returnSpec, int* status) {

  auto *areaIntegSpec = (double*) malloc(sizeof(double)*n_ener);
  CHECK_MALLOC_RET_STATUS(areaIntegSpec, status, areaIntegSpec)

  for (int jj=0; jj<n_ener; jj++) {
    areaIntegSpec[jj] = 0.0;
  }
  for (int ii=0; ii<returnSpec->nrad; ii++){
    double prop_area = calc_proper_area_ring(returnSpec->rlo[ii],returnSpec->rhi[ii],spin);
    for (int jj=0; jj<n_ener; jj++) {
      areaIntegSpec[jj] += spec[ii][jj] * prop_area;
    }
  }
  return areaIntegSpec;
}


static int checkTemperatureProfile(double* temperature, int n){

  if (temperature[0] < 0.0) {
    printf("\n *** error: lowest bin of temperature profile is %f  and NOT >=0 as expected \n",
           temperature[0]);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

void fitsAddKeyRadialZone(char* fname, const returnSpec2D* returnSpec, int izone, int* status){

  fitsfile* fptr = nullptr;
  fits_open_table(&fptr, fname, READWRITE, status);

  double rlo = returnSpec->rlo[izone];
  double rhi = returnSpec->rhi[izone];
  fits_update_key(fptr, TDOUBLE, "rlo",&(rlo),"radial zone grid", status);
  fits_update_key(fptr, TDOUBLE, "rhi",&(rhi),"radial zone grid", status);
  fits_update_key(fptr, TINT, "izone",&izone,"radial zone grid", status);

  fits_close_file(fptr, status);


  //relxill_check_fits_error(status);
}


//  ======= TEST FUNCTIONS ========   //

// ========== //
TEST_CASE(" Check Temperature Profile", "[bbret]"){

  int status = EXIT_SUCCESS;

  std::vector<double> spinArray = { 0.85,0.851,0.86,0.87,0.89,0.9 };

  double Rout = 1000;

  returningFractions* dat = nullptr;
  for (double spin : spinArray){
    double Rin = kerr_rms(spin);
    dat = get_rrad_fractions(spin, Rin, Rout, &status);

   VecD temperature = getTemperatureProfileDiskZones(dat, Rin, 1.0, &status);

    REQUIRE( checkTemperatureProfile(temperature.data(), dat->nrad) == EXIT_SUCCESS);

    delete_returningFractions(&dat);
  }

}

// ========== //
TEST_CASE(" DiskBB Spectrum ", "[bbret]"){

  int status = EXIT_SUCCESS;
  setenv("DEBUG_RELXILL","1",1);

  //create an energy grid
  int n_ener = 50;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 20.0);

  // call the relline model
  double photar[n_ener];

  double Tin = 1.0;
  double spin = 0.998; // simply determines the radial grid here

  spec_diskbb(ener, photar, n_ener, Tin, spin, &status);

  fits_write_spec("!testrr-spec-diskbb.fits", ener, photar, n_ener, &status);

  REQUIRE( compareWithRefdata("refvalues_diskbb.fits", ener, photar, n_ener, &status) == EXIT_SUCCESS );
  REQUIRE( status == EXIT_SUCCESS );

  unsetenv("DEBUG_RELXILL");
}

// ========== //
TEST_CASE(" Return Radiation Black Body Spectrum", "[bbret]"){

  int status = EXIT_SUCCESS;

  /* create an energy grid */
  int n_ener = 200;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 50.0);

  /* call the relline model */
  double photar[n_ener];
  double photar0[n_ener];

  double Tin = 1.0;
  double spin = 0.998; // simply determines the radial grid here

  double Rin = kerr_rms(spin);
  double Rout = RMAX_RELRET;

  returnSpec2D *returnSpec = spec_returnrad_blackbody(ener, n_ener, Tin, Rin, Rout, spin, &status);

  assert(returnSpec != nullptr);

  double* photar_areaInteg = integSpecArea(returnSpec->specRet, n_ener, spin, returnSpec, &status);
  double* photar0_areaInteg = integSpecArea(returnSpec->specPri, n_ener, spin, returnSpec, &status);

  fits_rr_write_2Dspec("!testrr-rframe-rr-bbody.fits", returnSpec->specRet, ener, n_ener,
                       returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, &status);
  fits_rr_write_2Dspec("!testrr-rframe-prim-bbody.fits", returnSpec->specPri, ener, n_ener,
                       returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, &status);

  // refdata currently missing
  // REQUIRE( compareWithRefdata("refvalues_rrbbody_ret.fits", ener, photar_areaInteg, n_ener, &status) == EXIT_SUCCESS);
  // REQUIRE( compareWithRefdata("refvalues_rrbbody_pri.fits", ener, photar0_areaInteg, n_ener, &status ) == EXIT_SUCCESS);

  free(photar_areaInteg);
  free(photar0_areaInteg);

  delete_returnSpec2D(&returnSpec);

}

TEST_CASE(" Check loading of relxillBB table","[bbret]")
{
  int status = EXIT_SUCCESS;

  LocalModel lmod(ModelName::relxillBB);
  auto xill_params = lmod.xill_params();
  auto xill_table_params = get_xilltab_param(xill_params, &status);
  auto xill_spec = get_xillver_spectra_table(xill_table_params, &status);

  for (size_t ii = 0; ii < xill_spec->n_incl; ii++)
  {
    REQUIRE(calcSum(xill_spec->flu[ii], xill_spec->n_ener) > 1e-8);
  }

  free_xill_spec(xill_spec);
  free(xill_table_params);
}




TEST_CASE("Check if Second Evaluation returns the identical results", "[bbret]"){

  LocalModel lmod(ModelName::relxillBB);

  DefaultSpec default_spec1{};
  XspecSpectrum spec1 = default_spec1.get_xspec_spectrum();
  lmod.eval_model(spec1);

  DefaultSpec default_spec0{};
  XspecSpectrum spec0 = default_spec0.get_xspec_spectrum();
  lmod.eval_model(spec0);

  auto sum0 = calcSum(spec0.flux,spec0.num_flux_bins());
  REQUIRE( std::isnormal (sum0));
  auto sum1 = calcSum(spec1.flux,spec0.num_flux_bins());
    REQUIRE(fabs(sum0 - sum1) < PREC);
}

void test_boost_parameter(LocalModel& lmod) {
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  lmod.set_par(XPar::boost, 1);
  lmod.eval_model(spec);
  double sum_total = calcSum(spec.flux, spec.num_flux_bins());
  REQUIRE(sum_total > 1e-8);

  lmod.set_par(XPar::boost,-1);
  lmod.eval_model(spec);
  double sum_refl = calcSum(spec.flux, spec.num_flux_bins());
  REQUIRE(sum_refl > 1e-8);

  lmod.set_par(XPar::boost, 0);
  lmod.eval_model(spec);
  double sum_direct = calcSum(spec.flux, spec.num_flux_bins());
  REQUIRE(sum_direct > 1e-8);


  // check that boost parmameter works: refl+direct has to give total flux
  printf(" total:%.3e direct:%.3e refl:%.3e (diff:%.3e) \n",
         sum_total,
         sum_direct,
         sum_refl,
         fabs(sum_total - sum_direct - sum_refl) / sum_total);
  REQUIRE( fabs(sum_total - sum_direct - sum_refl)/sum_total < 1e-5);
}


TEST_CASE(" Evaluate RelxillBB (only black body)","[bbret-fail]") {

  int status = EXIT_SUCCESS;
  setenv("RELXILL_WRITE_OUTFILES","1",1);

  LocalModel lmod(ModelName::relxillBB);
  setenv("RELXILL_BBRET_NOREFL","1",1);
  test_boost_parameter(lmod);
  setenv("RELXILL_BBRET_NOREFL", "0", 1);
  test_boost_parameter(lmod);

  LocalModel lmod2(ModelName::relxillBBxill);
  test_boost_parameter(lmod2);

  unsetenv("RELXILL_WRITE_OUTFILES");
  unsetenv("RELXILL_BBRET_NOREFL");
}

TEST_CASE("Test different spin values", "[bbret-speed]"){

 std::vector<double> spinArray{0.35,0.71,0.801,0.99,0.998};
  LocalModel lmod(ModelName::relxillBB);

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  for ( double spin : spinArray) {
    lmod.set_par(XPar::a, spin);
    lmod.eval_model(spec);

    // std::cout << "spin " << spin << " : " << calcSum(spec.flux(), spec.num_flux_bins()) << std::endl;

    REQUIRE(calcSum(spec.flux, spec.num_flux_bins()) > 1e-8);
  }

}
