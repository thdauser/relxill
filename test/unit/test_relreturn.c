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

    Copyright 2020 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "test_relreturn.h"

#include "relbase.h"
#include "relreturn.h"
#include "test_relxill.h"
#include "relreturn_corona.h"
#include "writeOutfiles.h"

#ifndef RELXILL_SOURCE_DIR
#error: RELXILL_SOURCE_DIR is not defined
#endif

#define PREC 1e-6
#define PREC_REFDATA 1e-3
#define REFDATA_DIR "test/returnrad/refdata/"  // relativ to RELXILL_SOURCE_DIR (set when compiling)



// Function definitions
static double* integSpecArea(double** spec, int n_ener, double spin, const returnSpec2D *returnSpec, int* status);

static void fits_write_spec(char* fname, const double* ener, double* spec, int n_ener, int* status){

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
    relxill_check_fits_error(status);
    printf("  *** error : creating file %s failed\n", fname);
    CHECK_STATUS_VOID(*status);
  }

  const int tfields = 3;
  char *ttype[] = {"bin_lo", "bin_hi", "flux"};
  char *tform[] = {"1D","1D","1D"};

  if (fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, NULL, "spec", status)) {
    relxill_check_fits_error(status);
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

  if (fptr != NULL) { fits_close_file(fptr, status); }



}

static void testReturnRadTableRadialGrid(returnFracData *dat, const int *status) {

  CHECK_STATUS_VOID(*status);

  assert(dat->nrad == RETURNRAD_TABLE_NR);
  assert(dat->ng == RETURNRAD_TABLE_NG);

  assert(dat->rlo[0] > 1.0);
  assert(dat->rhi[dat->nrad - 1] <= 1000.1); // don't care about border effects so add the 0.1

  for (int ii=1; ii<dat->nrad; ii++){
    assert(fabs(dat->rhi[ii-1]-dat->rlo[ii])<1e-6);
  }

}

void test_norm_frac_e(const returnFracData *dat, int *status) {

  double kSumfrac;
  double kSumfrac_ref = 1.0;

  for (int jj=0; jj<dat->nrad; jj++) {

    kSumfrac=0.0;
    for (int ii = 0; ii < dat->nrad; ii++) {

      kSumfrac += dat->frac_e[ii][jj];
    }
    if (fabs(kSumfrac - kSumfrac_ref) > PREC) {
      RELXILL_ERROR("testing the normalization of FRAC_E failed (return radiation table)", status);
      printf(" expecting a normalization of %e, but found %e\n", kSumfrac_ref, kSumfrac);
    }
  }
}

void test_norm_frac_g(const returnFracData *dat, int *status) {

  double kSumfrac;
  double kSumfrac_ref = 1.0;

  for (int ii=0; ii<dat->nrad; ii++) {
    for (int jj=0; jj<dat->nrad; jj++) {

      kSumfrac = 0.0;
      for (int kk = 0; kk < dat->ng; kk++) {
        kSumfrac += dat->frac_g[ii][jj][kk];
      }

      if (fabs(kSumfrac-kSumfrac_ref)>PREC){
        RELXILL_ERROR("testing the normalization of FRAC_G failed (return radiation table)",status);
        printf (" [%02i][%02i] expecting a normalization of %e, but found %e\n", ii,jj,kSumfrac_ref, kSumfrac);
      }
    }
  }
}

static void testReturnRadTableFractionNormalization(returnFracData* dat, int* status){

  assert(dat->frac_e[0][0]!= 0);

  test_norm_frac_e(dat, status);
  test_norm_frac_g(dat, status);

}

static refSpecData *new_refSpecData(int nener, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  refSpecData *spec = malloc(sizeof(refSpecData));
  CHECK_MALLOC_RET_STATUS(spec, status, spec)

  spec->nener = nener;

  spec->elo = (double *) malloc(sizeof(double) * nener);
  CHECK_MALLOC_RET_STATUS(spec->elo, status, spec)
  spec->ehi = (double *) malloc(sizeof(double) * nener);
  CHECK_MALLOC_RET_STATUS(spec->ehi, status, spec)
  spec->flux = (double *) malloc(sizeof(double) * nener);
  CHECK_MALLOC_RET_STATUS(spec->flux, status, spec)

  return spec;
}

static void free_refSpecData(refSpecData *spec) {
  if (spec != NULL) {
    free(spec->elo);
    free(spec->ehi);
    free(spec->flux);
  }
}

static refSpecData *fits_read_refdata_spec(char *fname, int *status) {

  // open the reference data
  fitsfile *fptr = NULL;
  (fits_open_table(&fptr, fname, READONLY, status));
  relxill_check_fits_error(status);
  CHECK_STATUS_RET(*status, NULL);

  long nener;

  // get the number of rows
  if (fits_get_num_rows(fptr, &nener, status)) return NULL;

  refSpecData *spec = new_refSpecData((int) nener, status);

  int anynul = 0;
  double nullval = 0.0;
  LONGLONG nelem = (LONGLONG) nener;
  fits_read_col(fptr, TDOUBLE, 1, 1, 1, nelem, &nullval, spec->elo, &anynul, status);
  fits_read_col(fptr, TDOUBLE, 2, 1, 1, nelem, &nullval, spec->ehi, &anynul, status);
  fits_read_col(fptr, TDOUBLE, 3, 1, 1, nelem, &nullval, spec->flux, &anynul, status);

  if (*status != EXIT_SUCCESS) {
    printf(" *** error *** initializing of the RETURN RADIATION table %s failed \n", fname);
    free_refSpecData(spec);
  }

  if (fptr != NULL) {
    fits_close_file(fptr, status);
  }

  return spec;
}

static void compare_refSpecDat(refSpecData *spec, double *ener0, double *flux0, int nener0, int *status) {

  CHECK_STATUS_VOID(*status);

  if (spec->nener != nener0) {
    *status = EXIT_FAILURE;
    printf(" *** error : spectrum has %i bins, but expected to have %i bins\n", spec->nener, nener0);
    return;
  }

  for (int ii = 0; ii < nener0; ii++) {

    if ((spec->flux[ii] / flux0[ii]) > PREC_REFDATA && (fabs(spec->flux[ii] - flux0[ii]) > PREC_REFDATA)) {
      *status = EXIT_FAILURE;
      printf(" *** error : spectrum flux in bin [%.3f, %.3f] is %e, but expected to be %e \n",
             ener0[ii], ener0[ii + 1], flux0[ii], spec->flux[ii]);
    }

  }

}

static void compareWithRefdata(char *fname_ref, double *ener0, double *flux0, int nener0, int *status) {

  PRINT_RELXILL_TEST_MSG(fname_ref);

  char fullname_ref[1000];
  sprintf(fullname_ref, "%s/%s/%s", RELXILL_SOURCE_DIR, REFDATA_DIR, fname_ref);
  refSpecData *spec = fits_read_refdata_spec(fullname_ref, status);

  compare_refSpecDat(spec, ener0, flux0, nener0, status);

  print_relxill_test_result(*status);

}

static double* integSpecArea(double** spec, int n_ener, double spin, const returnSpec2D *returnSpec, int* status) {

  double* areaIntegSpec = malloc(sizeof(double)*n_ener);
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

  fitsfile* fptr = NULL;
  fits_open_table(&fptr, fname, READWRITE, status);

  fits_update_key(fptr, TDOUBLE, "rlo",&(returnSpec->rlo[izone]),"radial zone grid", status);
  fits_update_key(fptr, TDOUBLE, "rhi",&(returnSpec->rhi[izone]),"radial zone grid", status);
  fits_update_key(fptr, TINT, "izone",&izone,"radial zone grid", status);

  fits_close_file(fptr, status);


  //relxill_check_fits_error(status);
}

void writeZoneXillverAllIncidentReturnSpec(double Tshift, int indZone, int n_ener, double *ener, xillParam *xill_param,
                                           const returnSpec2D *returnSpec, rel_spec *rel_profile, int *status) {


  char* fname_base = "testrr-spec-rframe-bbody";
  char fname_zone[200];
  char fnameSingleBuffer[200];

  sprintf(fname_zone,"%s-izone%02i-Tshift%.2f-%%s.fits",fname_base, indZone, Tshift);


  double Tin = xill_param->kTbb;
  xill_param->kTbb *= Tshift;


  // -0- direct radiation
  sprintf(fnameSingleBuffer,fname_zone, "direct");

  double out_direct_flux[n_ener];
  getZoneDirectPrimaryFlux(xill_param, returnSpec, out_direct_flux, indZone);
  fits_write_spec(fnameSingleBuffer,ener, out_direct_flux, n_ener, status);

  fitsAddKeyRadialZone(fnameSingleBuffer, returnSpec, indZone, status);
  CHECK_STATUS_VOID(*status);

  // -1- incident radiation
  sprintf(fnameSingleBuffer,fname_zone, "incident");

  double out_return_flux[n_ener];
  getZoneIncidentReturnFlux(xill_param, returnSpec, out_return_flux, indZone);
  fits_write_spec(fnameSingleBuffer,ener, out_return_flux, n_ener, status);

  fitsAddKeyRadialZone(fnameSingleBuffer, returnSpec, indZone, status);
  CHECK_STATUS_VOID(*status);

  // -2- reflected radiation
  sprintf(fnameSingleBuffer,fname_zone, "reflect");

  double xillver_out[n_ener];
  getZoneReflectedReturnFluxDiskframe(xill_param, rel_profile, returnSpec, xillver_out, indZone, status);
  fits_write_spec(fnameSingleBuffer,ener, xillver_out, n_ener, status);
  fitsAddKeyRadialZone(fnameSingleBuffer, returnSpec, indZone, status);
  CHECK_STATUS_VOID(*status);


  // -3- primary radiation causing the reflection
  sprintf(fnameSingleBuffer,fname_zone, "reflectPrim");

  double* xillver_prim_out = scaledXillverPrimaryBBodyHighener(xill_param->kTbb, returnSpec->specRet[indZone],
                                                               returnSpec->ener, returnSpec->n_ener, status);
  fits_write_spec(fnameSingleBuffer,ener, xillver_prim_out, n_ener, status);
  fitsAddKeyRadialZone(fnameSingleBuffer, returnSpec, indZone, status);
  CHECK_STATUS_VOID(*status);
  free(xillver_prim_out);

  // reset temperature
  xill_param->kTbb = Tin;
}


void get_std_param_relxill_bbret(relParam** p_rel_param, xillParam** p_xill_param, int* status) {

  CHECK_STATUS_VOID(*status);

  // set the standard parameters
  int n_parameter = NUM_PARAM_RELXILLBBRET;
  double parameter[n_parameter];
  set_std_param_relxill_bbret(parameter);
  init_par_relxill_bbret(p_rel_param, p_xill_param, parameter, n_parameter, status);

}


//  ======= TEST FUNCTIONS ========   //


static void testReturnfractionsRadialGridWhenChangingRin(int* status){

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();


  double spin =0.9;
  double Rout = 1000;

  double Rin = kerr_rms(spin);
  returnFracIpol *dat = get_rr_fractions(spin, Rin, Rout, status);
  int nrad_rms  = dat->nrad;

  Rin *= 2;
  dat = get_rr_fractions(spin, Rin, Rout, status);
  int nrad_rfac2 = dat->nrad;

  if (nrad_rms>nrad_rfac2){
    *status = EXIT_SUCCESS;
  } else {
    *status = EXIT_FAILURE;
  }

    print_relxill_test_result(*status);

}

static void testReturnfractionsInterpolationForDifferentSpins(int* status){

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();


  int nspin = 1;
  double spin[] = { 0.86 };

  double Rout = 1000;

  for (int ii=0; ii<nspin; ii++){

    double Rin = kerr_rms(spin[ii]);
    returnFracIpol *dat = get_rr_fractions(spin[ii], Rin, Rout, status);

    assert(dat != NULL);

    printf(" this test is missing any useful actions \n");
    // TODO: Test missing here
  }

  print_relxill_test_result(*status);

}

static void testTemperatureProfile(int* status){

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();


  int nspin = 6;
  double spin[] = { 0.85,0.851,0.86,0.87,0.89,0.9 };

  double Rout = 1000;

  returnFracIpol* dat = NULL;
  for (int ii=0; ii<nspin; ii++){
    double Rin = kerr_rms(spin[ii]);
    dat = get_rr_fractions(spin[ii], Rin, Rout, status);

    double* temperature = getTemperatureProfileDiskZones(dat, Rin, 1.0, status);

    *status = checkTemperatureProfile(temperature, dat->nrad);
    CHECK_STATUS_BREAK(*status);

    free(temperature);
  }

  print_relxill_test_result(*status);

}

static void testReturnRadTableNormlizationAndRadialGrid(returnTable* tab, int* status){

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  for (int ii=0; ii<tab->nspin; ii++) {

    testReturnRadTableRadialGrid(tab->retFrac[ii], status);
    if (*status!=EXIT_SUCCESS){
      printf(" *** FAILED: testing radial grid\n");
      return;
    }

    testReturnRadTableFractionNormalization(tab->retFrac[ii], status);
    if (*status != EXIT_SUCCESS) {
      printf(" *** FAILED: testing FRAC_E normalization\n");
      return;
    }

  }

  print_relxill_test_result(*status);

}

static void testDiskbbSpec(int *status) {

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  /* create an energy grid */
  int n_ener = 50;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 20.0);

  /* call the relline model */
  double photar[n_ener];

  double Tin = 1.0;
  double spin = 0.998; // simply determines the radial grid here

  spec_diskbb(ener, photar, n_ener, Tin, spin, status);

  // compareWithRefdata("refvalues_diskbb.fits", ener, photar, n_ener, status);

  print_relxill_test_result(*status);
}

static void testReturnRadBBodySpec(int* status){

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG(": \n");

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

  returnSpec2D *returnSpec = spec_returnrad_blackbody(ener, photar, photar0, n_ener, Tin, Rin, Rout, spin, status);

  assert(returnSpec != NULL);

  double* photar_areaInteg = integSpecArea(returnSpec->specRet, n_ener, spin, returnSpec, status);
  double* photar0_areaInteg = integSpecArea(returnSpec->specPri, n_ener, spin, returnSpec, status);

  fits_rr_write_2Dspec("!testrr-rframe-rr-bbody.fits", returnSpec->specRet, ener, n_ener,
                       returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, NULL, status);
  fits_rr_write_2Dspec("!testrr-rframe-prim-bbody.fits", returnSpec->specPri, ener, n_ener,
                       returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, NULL, status);

  compareWithRefdata("refvalues_rrbbody_ret.fits", ener, photar_areaInteg, n_ener, status);
  compareWithRefdata("refvalues_rrbbody_pri.fits", ener, photar0_areaInteg, n_ener, status);

  free(photar_areaInteg);
  free(photar0_areaInteg);


}

static void testSingleZoneRframeReflectSpectrum(int* status){

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();


  // get a standard grid for the convolution (is rebinned later to the input grid)
  int n_ener;
  double *ener;
  get_std_relxill_energy_grid(&n_ener, &ener, status);

  // set the standard parameters
  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;
  get_std_param_relxill_bbret(&rel_param, &xill_param, status);
  xill_param->refl_frac = -1.0;

  xillTable *xill_tab = NULL;
  get_init_xillver_table(&xill_tab, xill_param, status);

  returnSpec2D *returnSpec = spec_returnrad_blackbody(ener, NULL, NULL, n_ener, xill_param->kTbb, rel_param->rin,
                                                      rel_param->rout, rel_param->a, status);

  double *radialGrid = getRadialGridFromReturntab(returnSpec, status);
  rel_spec *rel_profile = relbase_multizone(ener, n_ener, rel_param, xill_tab, radialGrid, returnSpec->nrad, status);

  enum{ nzones = 4};
  enum{ nshift = 2};
  int indZone[nzones] = {10,20,30,40};
  double Tshift[nshift] = {1.0,1.4};

  for (int ii=0; ii<nzones; ii++) {
    for (int jj=0; jj<nshift; jj++) {
      writeZoneXillverAllIncidentReturnSpec(Tshift[jj], indZone[ii], n_ener, ener, xill_param, returnSpec, rel_profile, status);
    }
  }

  printf(" [ !! missing comparison !! ]  ");

  print_relxill_test_result(*status);

}

static void test_rr_bbody_lmod(int* status){

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  double Tin = 1.0;
  double spin = 0.86;


  /* create an energy grid */
  int n_ener = 1000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 50.0);

  // set the standard parameters
  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;
  get_std_param_relxill_bbret(&rel_param, &xill_param, status);

  xill_param->kTbb = Tin;
  rel_param->a = spin;
  rel_param->rin = kerr_rms(spin);
  rel_param->incl = 60.0/180.0*M_PI;

  // create space for the output spectrum
  double photar[n_ener];


  xill_param->refl_frac = -1.0;
  relxill_bb_kernel(ener, photar, n_ener, xill_param, rel_param, 1, status);
  fits_write_spec("!testrr-spec-rr-bbody.fits",ener, photar, n_ener, status);


  xill_param->refl_frac = 0.0;
  relxill_bb_kernel(ener, photar, n_ener, xill_param, rel_param, 0, status);
  fits_write_spec("!testrr-spec-prim-bbody.fits",ener, photar, n_ener, status);

  xill_param->refl_frac = -1.0;
  relxill_bb_kernel(ener, photar, n_ener, xill_param, rel_param, 0, status);
  fits_write_spec("!testrr-spec-rr-bbody-reflect.fits",ener, photar, n_ener, status);

  relxill_bb_kernel(ener, photar, n_ener, xill_param, rel_param, 0, status);
  fits_write_spec("!testrr-spec-rr-bbody-reflect.fits",ener, photar, n_ener, status);


  print_relxill_test_result(*status);
}


static void testSecondEvalReturnsIdenticalResults(int* status) {

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  /* set the parameters */
  enum {n_param = NUM_PARAM_RELXILLBBRET};
  double inp_par[n_param];
  set_std_param_relxill_bbret(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  Spectrum *spec0  = getNewSpec(0.1, 10, 1000, status);
  Spectrum *spec1 = getNewSpec(0.1, 10, 1000, status);
  CHECK_RELXILL_DEFAULT_ERROR(status);

  inp_par[9] = 1.0;

  tdrelxillbbret(spec0->ener, spec0->nbins, spec0->flux, inp_par, n_param, status);
  CHECK_RELXILL_DEFAULT_ERROR(status);

  inp_par[9] = 1.0;

 // tdrelxillbbret(spec0->ener, spec0->nbins, spec0->flux, inp_par, n_param, status);
  CHECK_RELXILL_DEFAULT_ERROR(status);

  inp_par[9] = 1.0;

  tdrelxillbbret(spec1->ener, spec1->nbins, spec1->flux, inp_par, n_param, status);
  CHECK_RELXILL_DEFAULT_ERROR(status);


  for (int ii=0; ii<spec0->nbins; ii++){
    if ( fabs(spec0->flux[ii] - spec1->flux[ii]) > PREC){
      *status = EXIT_FAILURE;
      CHECK_RELXILL_DEFAULT_ERROR(status);
      printf("    energ bin [%e,%e]  :  flux0=%e not equal to flux1=%e \n",
          spec0->ener[ii],spec0->ener[ii+1], spec0->flux[ii], spec1->flux[ii]);
      return;
    }
  }


  free_Spectrum(&spec0);
  free_Spectrum(&spec1);

  print_relxill_test_result(*status);

}

static void testEvaluateRelxillbbret(int* status) {

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG("(switching on DEBUG=1)\n");
  putenv("DEBUG_RELXILL=1");

  /* set the parameters */
  enum {n_param = NUM_PARAM_RELXILLBBRET};
  double inp_par[n_param];
  set_std_param_relxill_bbret(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */

  Spectrum* spec = getNewSpec(0.1,100,1000, status);
  CHECK_RELXILL_DEFAULT_ERROR(status);

  tdrelxillbbret(spec->ener, spec->nbins, spec->flux, inp_par, n_param, status);
  CHECK_RELXILL_DEFAULT_ERROR(status);

  fits_write_spec("!testrr-spec-relxillbbret.fits",spec->ener, spec->flux, spec->nbins, status);

  putenv("DEBUG_RELXILL=0");

  print_relxill_test_result(*status);

}

static void testBBretDifferentSpinValues(int *status) {

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  Spectrum *spec = getNewSpec(0.1, 100, 1000, status);
  CHECK_RELXILL_DEFAULT_ERROR(status);


  /* set the parameters */
  enum { n_param = NUM_PARAM_RELXILLBBRET };
  double inp_par[n_param];
  set_std_param_relxill_bbret(inp_par);
  CHECK_STATUS_VOID(*status);


  int nspin = 5;
  double spinArray[] = {0.35,0.71,0.801,0.99,0.998};

  for (int ii=0; ii< nspin; ii++){
    inp_par[0] = spinArray[ii];
    printf("\n    -> testing spin=%.3e ", spinArray[ii]);
    tdrelxillbbret(spec->ener, spec->nbins, spec->flux, inp_par, n_param, status);

    if (calcSum(spec->flux, spec->nbins) < 0) {
      printf("\n  *** error: testing spin=%.3e failed (summed output flux = %e ) \n",
             spinArray[ii], calcSum(spec->flux, spec->nbins));
      *status = EXIT_FAILURE;
      break;
    }
  }

  print_relxill_test_result(*status);

}

void save_emisProfile(char *fname, emisProfile *emis) {
  write_data_to_file(fname, emis->re, emis->emis, emis->nr);
}


/* ******* test_coronaRet ******** */

static void set_std_param_relxilllpRet(double *inp_par) {
  inp_par[0] = 3;   // height
  inp_par[1] = 0.998; // a
  inp_par[2] = 60.0;  // incl
  inp_par[3] = -1.0;  // rin
  inp_par[4] = 1000.;  // rout
  inp_par[5] = 0.0;    // redshift
  inp_par[6] = 2.1;   // pl Index
  inp_par[7] = 0.0;   // logxi
  inp_par[8] = 1.0;   // Afe
  inp_par[9] = 300.0; // Ecut
  inp_par[10] = 3.0;   // refl_frac
  inp_par[11] = 0.0;   // fixReflFrac
  inp_par[12] = 1;
}

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


static void testRebinEmisProfiles(int *status) {

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  double precRebinCoarse = 0.05;

  relParam* rel_param =  get_std_param_rellinelp(status);

  RelSysPar *sysPar = get_system_parameters(rel_param, status);

  returnFracIpol *dat = get_rr_fractions(rel_param->a, rel_param->rin, rel_param->rout, status);
  double rmean[dat->nrad]; // descending grid
  for (int ii = 0; ii < dat->nrad; ii++) {
    rmean[dat->nrad - ii - 1] = 0.5 * (dat->rlo[ii] + dat->rhi[ii]);
  }
  emisProfile *emisCoarse = calc_emis_profile(rmean, dat->nrad, rel_param, status);
  // invertArray(emisCoarse->re, emisCoarse->nr);
  // invertArray(emisCoarse->emis, emisCoarse->nr);

  emisProfile *emisFineReference = calc_emis_profile(sysPar->re, sysPar->nr, rel_param, status);
  emisProfile *emisRebin = new_emisProfile(sysPar->re, sysPar->nr, status);
  interpolEmisProfile(emisRebin, emisCoarse, status);

  if (*status == EXIT_SUCCESS) {

    save_emisProfile("test_emisCoarse.dat", emisCoarse);
    save_emisProfile("test_emisFineReference.dat", emisFineReference);
    save_emisProfile("test_emisRebin.dat", emisRebin);

    for (int ii = 20; ii < emisRebin->nr;
         ii++) {  // last bin in coarse grid deviates, so skip, as it is not the interpolation
      if (is_debug_run()) {
        printf(" rad: %.3e :  %e (ref=%e, ratio=%e) \n",
               emisRebin->re[ii], emisRebin->emis[ii], emisFineReference->emis[ii],
               emisRebin->emis[ii] / emisFineReference->emis[ii]);
      }
      assert(fabs(emisRebin->emis[ii] / emisFineReference->emis[ii] - 1) < precRebinCoarse);
    }

  }

  print_relxill_test_result(*status);

}



static void returnEmisProfileLoaded(int *status) {

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  relParam *rel_param = get_std_param_rellinelp(status);

  rel_param->return_rad = 0;
  RelSysPar *sysPar = get_system_parameters(rel_param, status);
  assert(sysPar->emisReturn == NULL);

  rel_param->return_rad = 1;
  sysPar = get_system_parameters(rel_param, status);
  assert(sysPar->emisReturn != NULL);

  print_relxill_test_result(*status);

}

static void returnEmisLineProfileCalculation(int *status) {

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  relParam *rel_param = get_std_param_rellinelp(status);
  rel_param->a = 0.998;
  rel_param->height = 5;

  //rel_param->return_rad=0;

  Spectrum *spec = getNewSpec(0.05, 10, 1000, status);

  rel_param->return_rad = 1;
  rel_spec *rel_profile = relbase(spec->ener, spec->nbins, rel_param, NULL, status);
  CHECK_STATUS_VOID(*status);

  assert(rel_profile->n_zones == 1);
  assert(calcSum(rel_profile->flux[0], rel_profile->n_ener) > 1e-8);

  print_relxill_test_result(*status);

}

static void evalLmodRelxilllpRet(int *status) {

  CHECK_STATUS_VOID(*status);
  PRINT_RELXILL_TEST_MSG_DEFAULT();

  enum { n_param = NUM_PARAM_RELXILLLPRET };
  double inp_par[n_param];
  set_std_param_relxilllpRet(inp_par);
  CHECK_STATUS_VOID(*status);

  Spectrum *spec = getNewSpec(0.1, 100, 1000, status);
  CHECK_RELXILL_DEFAULT_ERROR(status);

  tdrelxilllpret(spec->ener, spec->nbins, spec->flux, inp_par, n_param, status);
  CHECK_RELXILL_DEFAULT_ERROR(status);

  fits_write_spec("!testrr-spec-relxilllpret.fits", spec->ener, spec->flux, spec->nbins, status);

  print_relxill_test_result(*status);

}


//  ======= MAIN ========   //

void test_general(int *status) {

  returnTable *tab = get_returnRadTable(status);

  testReturnRadTableNormlizationAndRadialGrid(tab, status);

  testReturnfractionsRadialGridWhenChangingRin(status);

  testReturnfractionsInterpolationForDifferentSpins(status);

}

void test_bbodyRet(int *status) {

  testTemperatureProfile(status);
  testDiskbbSpec(status);
  testReturnRadBBodySpec(status);
  testSingleZoneRframeReflectSpectrum(status);
  test_rr_bbody_lmod(status);
  testSecondEvalReturnsIdenticalResults(status);
  testEvaluateRelxillbbret(status);
  testBBretDifferentSpinValues(status);

}

void test_coronaRet(int *status) {

  testRebinEmisProfiles(status);
  //evalLmodRelxilllpRet(status);
  //returnEmisProfileLoaded(status);
  returnEmisLineProfileCalculation(status);

}

void test_relreturn(void) {

  int status = EXIT_SUCCESS;

  print_version_number();
  printf("\n### Testing RETURN RADIATION ###\n");

  test_general(&status);

  test_coronaRet(&status);

  //  test_bbodyRet(&status);

}

int main(int argc, char *argv[]) {

  test_relreturn();

}
