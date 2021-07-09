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

#include "relreturn.h"

#include "common.h"
#include "relbase.h"
#include "relutility.h"
#include "relphysics.h"
#include "relreturn_datastruct.h"
#include "relmodels.h"
//#include "relxill.h"

#define LIM_GFAC_RR_BBODY 0.001 // difference between gmin and gmax, above which the energy shift is taken into account

// Function definitions
static double **get_returnrad_specs(double *ener_inp,
                                    int nener_inp,
                                    returningFractions *dat,
                                    const double *temperature,
                                    int *status);
static double **get_bbody_specs(double *ener, int nener, returningFractions *dat, double *temperature, int *status);

void normalizeFluxRrad(int nrad, int nener, const double *ener, double **spec);
double getEmissivityNormFactor(returningFractions *dat, int nener, const double *ener, double **spec);

// Program Code

double* getTemperatureProfileDiskZones(returningFractions* dat, double Rin, double Tin, int* status){
  return get_tprofile(dat->rlo, dat->rhi, dat->nrad, Rin, Tin, TPROFILE_ALPHA, status);
}

returnSpec2D *spec_returnrad_blackbody(double *ener, double *spec, double *spec_prim, int nener,
                                       double Tin, double Rin, double Rout, double spin, int *status) {
/*
 * MAIN FUNCTION: returns the returning black body emission in the frame of the disk
 * output: returnSpec (containing primary and returning spectrum)
 * input:
 */

  CHECK_STATUS_RET(*status, NULL);

  // 1 - get the fractions from the table (plus interpolation to current parameters)
  returningFractions *dat = get_rrad_fractions(spin, Rin, Rout, status);

  // 2 - temperature profile of the whole disk
  double *temperature = getTemperatureProfileDiskZones(dat, Rin, Tin, status);

  // 3 - get primary and returning 2D spectrum
  double **spec_rr_zones = get_returnrad_specs(ener, nener, dat, temperature, status);
  double **spec_prim_zones = get_bbody_specs(ener, nener, dat, temperature, status);

  // normalizeSpectraOverTheFullDisk(spec_prim_zones, spec_rr_zones, ener, nener, dat);

  returnSpec2D *returnSpec = getReturnradOutputStructure(dat, spec_rr_zones, spec_prim_zones, ener, nener, status);

  // get the primary spectrum for the given radial grid and temperature profile
  if (spec!=NULL) { // TODO: decide if we need this spec
    sum_2Dspec(spec, spec_rr_zones, nener, dat->nrad, status);
  }
  if (spec_prim!=NULL) { // TODO: decide if we need this spec
    sum_2Dspec(spec_prim, spec_prim_zones, nener, dat->nrad, status);
  }

  if (is_debug_run()) {
    fits_rr_write_2Dspec("!debug-testrr-rframe-rr-bbody.fits", returnSpec->specRet, ener, nener,
                         dat->rlo, dat->rhi, dat->nrad, dat, status);
    fits_rr_write_2Dspec("!debug-testrr-rframe-prim-bbody.fits", returnSpec->specPri, ener, nener,
                         dat->rlo, dat->rhi, dat->nrad, dat, status);
  }

  free(temperature);

  return returnSpec;

}

/*** Routines for the Black Body Case ***/

static void calc_rr_bbspec_gzone(double *ener, int nener, double *spec, double temp,
                                 double *gfac, int ng, const double *frac_g) {
  // loop over all gfac and calculate the spectrum for each
  double spec_g[nener];
  for (int kk = 0; kk < ng; kk++) {
    bbody_spec(ener, nener, spec_g, temp, gfac[kk]);
    for (int jj = 0; jj < nener; jj++) {
      spec[jj] += spec_g[jj] * frac_g[kk];
    }
  }

}


static void calc_rr_bbspec_ring(double* ener, double* spec, int nener, int irad, const double* temp, returningFractions* dat, const int* status){

  CHECK_STATUS_VOID(*status);

  for (int jj=0; jj<nener; jj++){
    spec[jj] = 0.0;
  }

  double gfac[dat->tabData->ng];
  double spec_r[nener];
  for (int ii=0; ii<dat->nrad; ii++){  // loop over all radial zones

    for (int jj=0; jj<nener; jj++){
      spec_r[jj] = 0.0;
    }

    get_gfac_grid(gfac, dat->tabData->gmin[irad][ii], dat->tabData->gmax[irad][ii], dat->tabData->ng);

    // only calculate it for a significant redshift, otherwise return BBODY without gshift
    if ( fabs(dat->tabData->gmax[irad][ii] - dat->tabData->gmin[irad][ii]) > LIM_GFAC_RR_BBODY ) {
      calc_rr_bbspec_gzone(ener, nener, spec_r, temp[ii], gfac, dat->tabData->ng, dat->tabData->frac_g[irad][ii]);
    } else {
      bbody_spec(ener, nener, spec_r, temp[ii], dat->tabData->gmin[irad][ii]); // caveat, gmean might not be 1.0
    }

    // apply the correct fractions and add it to the zone output spectrum
    for (int jj = 0; jj < nener; jj++) {
      spec[jj] += spec_r[jj] * dat->tabData->f_ret[ii] * dat->frac_i[irad][ii];
    }

  }
}



static double **get_returnrad_specs(double *ener_inp, int nener_inp, returningFractions *dat,
                                    const double *temperature, int *status) {
/* returns: - 2D-spectral array, unit is cts/bin  (xspec standard)
 *          - each zone is integrated over the GR proper area
 *          - normalization
 * input:   - temperature profile (in keV)
 *          - interpolated return radiation fractions
 */


  CHECK_STATUS_RET(*status, NULL);

  int nener;
  double *ener;
  get_std_bbody_energy_grid(&nener, &ener, status);

  double **spec_zones = new_specZonesArr(nener_inp, dat->nrad, status);

  assert(dat->frac_i[0] != NULL);

  double spec[nener];
  for (int ii = 0; ii < dat->nrad; ii++) {
    calc_rr_bbspec_ring(ener, spec, nener, ii, temperature, dat, status);  // return ph / cm²/s/keV ( not bin integ.)

    rebin_mean_flux(ener_inp, spec_zones[ii], nener_inp, ener, spec, nener, status);
  }

  normalizeFluxRrad(dat->nrad, nener_inp, ener_inp, spec_zones);


  return spec_zones;
}


static double **get_bbody_specs(double *ener, int nener, returningFractions *dat, double *temperature, int *status) {
/* returns: - 2D-spectral array, unit is cts/bin  (xspec standard)
 *          - outside the energy band EMIN_XILLVER to EMAX_XILLVER the function is 0
 * input: temperature profile (in keV)
 */

  CHECK_STATUS_RET(*status, NULL);

  double **spec_array = new_specZonesArr(nener, dat->nrad, status);

  for (int ii = 0; ii < dat->nrad; ii++) {

    bbody_spec(ener, nener, spec_array[ii], temperature[ii], 1.0);
    for (int jj = 0; jj < nener; jj++) {
      if (ener[jj] < EMIN_XILLVER || ener[jj + 1] > EMAX_XILLVER) {
        spec_array[ii][jj] = 0.0;
      }
    }
  }

  normalizeFluxRrad(dat->nrad, nener, ener, spec_array);

  return spec_array;
}


void normalizeFluxRrad(int nrad, int nener, const double *ener, double **spec) {

  for (int ii = 0; ii < nrad; ii++) {
    for (int jj = 0; jj < nener; jj++) {
      spec[ii][jj] *= (ener[jj + 1] - ener[jj]);  // now make it cts/bin
    }
  }
}


double getEmissivityNormFactor(returningFractions *dat, int nener, const double *ener, double **spec) {
/*
 *  Input: bin integrated spectrum  [cts /bin/cm^2]
 */

double sumRadius[dat->nrad];
  double sumTotalSpec = 0.0;

  for (int ii = 0; ii < dat->nrad; ii++) {
    sumRadius[ii] = 0.0;
    for (int jj = 0; jj < nener; jj++) {
      sumRadius[ii] += spec[ii][jj] *  dat->proper_area_ring[ii]; // integrate spectrum over the whole disk
    }

    sumTotalSpec += sumRadius[ii];
  }

  return sumTotalSpec;
}


void fits_rr_write_2Dspec(char *fname, double **spec_arr, double *ener, int nener,
                          double* rlo, double* rhi, int nrad, returningFractions *dat, int *status) {

  CHECK_STATUS_VOID(*status);

  // open the fits file
  fitsfile *fptr;
  if (fits_create_file(&fptr, fname, status)) {
    relxill_check_fits_error(status);
    printf("   creating file %s failed\n", fname);
    CHECK_STATUS_VOID(*status);
  }


  int n1 = nrad;
  int n2 = nener;
  char dim1[10];
  char dim2[10];
  sprintf(dim1, "%iD", (int) 1);
  sprintf(dim2, "%iD", (int) n2);

  if (dat!=NULL) {
    const int tfields = 6;
    char *ttype[] = {"rlo", "rhi", "ener", "spec", "fraci", "fret"};
    char *tform[] = {dim1, dim1, dim2, dim2, "50D", dim1};
    fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, NULL, "2Dspec", status);

  } else {
    const int tfields = 4;
    char *ttype[] = {"rlo", "rhi", "ener", "spec"};
    char *tform[] = {dim1, dim1, dim2, dim2};
    fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, NULL, "2Dspec", status);
  }

   relxill_check_fits_error(status);
   CHECK_STATUS_VOID(*status);

  int firstrow = 1;  /* first row in table to write   */
  int firstelem = 1;  /* first element in row  */

  fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, n1, rlo, status);
  fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, n1, rhi, status);
  for (int ii = 0; ii < n1; ii++) {
    fits_write_col(fptr, TDOUBLE, 3, firstrow + ii, firstelem, n2, ener, status);
    fits_write_col(fptr, TDOUBLE, 4, firstrow + ii, firstelem, n2, spec_arr[ii], status);
  }

  if (dat!=NULL) {
    for (int ii = 0; ii < n1; ii++) {
      fits_write_col(fptr, TDOUBLE, 5, firstrow + ii, firstelem, n1, dat->frac_i[ii], status);
    }
    fits_write_col(fptr, TDOUBLE, 6, firstrow, firstelem, n1, dat->tabData->f_ret, status);
    if( fits_write_key(fptr, TDOUBLE, "SPIN", &(dat->a),NULL, status) ) {}
  }


  relxill_check_fits_error(status);

  if (fptr != NULL) { fits_close_file(fptr, status); }

}


void spec_diskbb(double* ener, double* spec, int n, double Tin,  double spin, int* status) {

  CHECK_STATUS_VOID(*status);

  returningFractions *dat = get_rrad_fractions(spin, -1, RMAX_RELRET, status);

  double *temperature =
      get_tprofile(dat->rlo, dat->rhi, dat->nrad, 0.5 * (dat->rlo[0] + dat->rhi[0]), Tin, TPROFILE_DISKBB, status);

  double **spec_arr = get_bbody_specs(ener, n, dat, temperature, status);  // spectra integrated over the ring area

  /* need area correction for diskbb */
  for (int jj = 0; jj < n; jj++) {
    for (int ii = 0; ii < dat->nrad; ii++) {
      // need to divide by the GR area and multiply by the normal ring area to mimick a diskbb spectrum
      spec_arr[ii][jj] *=
          (M_PI * (pow(dat->rhi[ii], 2) - pow(dat->rlo[ii], 2))) / dat->proper_area_ring[ii];
    }
  }

  if (is_debug_run()) {
    fits_rr_write_2Dspec("!debug-testrr-diskbb.fits", spec_arr, ener, n, dat->rlo, dat->rhi, dat->nrad, dat, status);
  }

  sum_2Dspec(spec, spec_arr, n, dat->nrad, status);

  free_2d(&spec_arr, dat->nrad);
  free(temperature);
}



double *getRadialGridFromReturntab(returnSpec2D *spec, int* status) {

  double* rgrid = malloc(sizeof(double)*(spec->nrad+1)); // we use n+1 grid points
  CHECK_MALLOC_RET_STATUS(rgrid, status, rgrid)

  for (int ii=0; ii<spec->nrad; ii++){
    rgrid[ii] = spec->rlo[ii];
  }
  rgrid[spec->nrad] = spec->rhi[spec->nrad-1];

  return rgrid;
}


static double getXillverNormFactorFromPrimarySpectrum(double* spec, double* ener, int n_ener, int* status){

  CHECK_STATUS_RET(*status,0.0);

  EnerGrid* egrid = get_stdXillverEnergygrid(status);

  double* xillverInputSpec = malloc(sizeof(double)*egrid->nbins);
  CHECK_MALLOC_RET_STATUS(xillverInputSpec, status, 0.0)

  rebin_spectrum(egrid->ener, xillverInputSpec, egrid->nbins, ener, spec, n_ener);

  // divide by the primary normalization factor, to get the scaling of the xillver reflection spectrum
  double normFactorXill = calcNormWrtXillverTableSpec(xillverInputSpec, egrid->ener, egrid->nbins, status);

  free(xillverInputSpec);

  return normFactorXill;
}

double calcXillverNormfacRetrad2BoodyAtHighenergy(double kTbb, double * spec_in, double* spec_bb, double* ener, int n_ener, int* status){
  /* calculate the normalization factor between the primary returning radiation and the black body (xillver primary
   * spectra) at high energies
   */

  const double ELO_LIMIT_HIGHENER = 4*kTbb;

  double normReturnSpec = calcSumInEnergyBand(spec_in,n_ener, ener, ELO_LIMIT_HIGHENER , EMAX_XILLVER_NORMALIZATION);

  double normBbodySpec = calcSumInEnergyBand(spec_bb,n_ener, ener, ELO_LIMIT_HIGHENER, EMAX_XILLVER_NORMALIZATION);

  return normReturnSpec / normBbodySpec;
}

void getZoneReflectedReturnFluxDiskframe(xillParam *xill_param, rel_spec* rel_profile, const returnSpec2D *returnSpec,
                                         double *xill_flux_returnrad, int izone, int* status) {

  // assumption: we use Tin for all zones
  getNormalizedXillverSpec(xill_flux_returnrad, returnSpec->ener, returnSpec->n_ener, xill_param,
                           rel_profile->rel_cosne->dist[izone], status);


  double xillverReflectionNormFactor = getXillverNormFactorFromPrimarySpectrum(returnSpec->specRet[izone], returnSpec->ener, returnSpec->n_ener, status);
  CHECK_STATUS_VOID(*status);

  double* xillver_prim_out = getXillverPrimaryBBodyNormalized(xill_param->kTbb, returnSpec->specRet[izone],
      returnSpec->ener, returnSpec->n_ener, status);

  double normfacMatchAtHighEnergies = calcXillverNormfacRetrad2BoodyAtHighenergy(xill_param->kTbb,
                                                                                 returnSpec->specRet[izone], xillver_prim_out,
                                                                                 returnSpec->ener,
                                                                                 returnSpec->n_ener, status);
  free(xillver_prim_out);

  for (int jj = 0; jj < returnSpec->n_ener; jj++) {
    xill_flux_returnrad[jj] *= fabs(xill_param->refl_frac) ;
    xill_flux_returnrad[jj] *= xillverReflectionNormFactor * normfacMatchAtHighEnergies;

    if (xill_param->refl_frac >= 0) {
      xill_flux_returnrad[jj] += returnSpec->specPri[izone][jj];
    }
  }

}

void getZoneIncidentReturnFlux(xillParam *xill_param, const returnSpec2D *returnSpec, double *returnFlux, int ii) {

  for (int jj = 0; jj < returnSpec->n_ener; jj++) {
    returnFlux[jj] = returnSpec->specRet[ii][jj] * fabs(xill_param->refl_frac);
    if (xill_param->refl_frac >= 0) {
      returnFlux[jj] += returnSpec->specPri[ii][jj];
    }
  }

}

void getZoneDirectPrimaryFlux(xillParam *xill_param, const returnSpec2D *returnSpec, double *returnFlux, int ii) {

  for (int jj = 0; jj < returnSpec->n_ener; jj++) {
    returnFlux[jj] =  returnSpec->specPri[ii][jj];
  }

}


double* getXillverPrimaryBBodyNormalized(double kTbb, double* spec_in, double* ener, int n_ener, int* status){

  double* spec_out = malloc(sizeof(double)*n_ener);
  CHECK_MALLOC_RET_STATUS(spec_out, status, spec_out)

  double xillverReflectionNormFactor = getXillverNormFactorFromPrimarySpectrum(spec_in, ener, n_ener, status);

  bbody_spec(ener, n_ener, spec_out, kTbb, 1.0);

  for (int jj = 0; jj < n_ener; jj++) {
      spec_out[jj] *= (ener[jj + 1] - ener[jj]);
  }

  double bbodyNormFactor = getXillverNormFactorFromPrimarySpectrum(spec_out, ener, n_ener, status);

  for (int jj = 0; jj < n_ener; jj++) {
    spec_out[jj] *= xillverReflectionNormFactor / bbodyNormFactor;
  }

  return spec_out;
}

double* scaledXillverPrimaryBBodyHighener(double kTbb, double* spec_in, double* ener, int n_ener, int* status){

  double* xill_out_prim = getXillverPrimaryBBodyNormalized(kTbb, spec_in,ener, n_ener, status);

  double normFac = calcXillverNormfacRetrad2BoodyAtHighenergy(kTbb, spec_in, xill_out_prim, ener, n_ener, status);

  for (int jj = 0; jj < n_ener; jj++) {
    xill_out_prim[jj] *= normFac;
  }

  return xill_out_prim;
}


static void setLowValuesToZero(double* spec, int n, double lowestModelValue){

  double maxVal =  0.0;

  for (int ii=0; ii<n; ii++){
    maxVal = fmax(maxVal, spec[ii]);
  }

  for (int ii=0; ii<n; ii++) {
    if (spec[ii] < maxVal*lowestModelValue){
      spec[ii] = 0.0;
    }
  }
}

static void setValuesOutsideToZero(double* spec, const double* ener, int n){

  for (int ii=0; ii<n; ii++) {
    if (ener[ii] < EMIN_XILLVER || ener[ii+1]>EMAX_XILLVER){
      spec[ii] = 0.0;
    }
  }
}


void relxill_bb_kernel(double *ener_inp, double *spec_inp, int n_ener_inp, xillParam *xill_param, relParam *rel_param,
    int noXillverReflection, int *status) {

  CHECK_STATUS_VOID(*status);
  assert(xill_param->model_type == MOD_TYPE_RELXILLBBRET);

  // get a standard grid for the convolution (is rebinned later to the input grid)
  int n_ener;
  double *ener;
  get_std_relxill_energy_grid(&n_ener, &ener, status);

  xillTable *xill_tab = NULL;
  get_init_xillver_table(&xill_tab, xill_param, status);

  returnSpec2D *returnSpec = spec_returnrad_blackbody(ener, NULL, NULL, n_ener, xill_param->kTbb, rel_param->rin,
                                                      rel_param->rout, rel_param->a, status);

  double *radialGrid = getRadialGridFromReturntab(returnSpec, status);
  rel_spec *rel_profile = relbase_multizone(ener, n_ener, rel_param, xill_tab, radialGrid, returnSpec->nrad, status);

  // ========== //
  double single_spec_inp[n_ener_inp];

  double **spec_conv_out = malloc(sizeof(double *) * returnSpec->nrad);
  CHECK_MALLOC_VOID_STATUS(spec_conv_out, status)

  double **xillver_out = malloc(sizeof(double *) * returnSpec->nrad);
  CHECK_MALLOC_VOID_STATUS(spec_conv_out, status)

  double **xillver_prim_out = malloc(sizeof(double *) * returnSpec->nrad);
  CHECK_MALLOC_VOID_STATUS(spec_conv_out, status)
  // ========== //

  double Tin = xill_param ->kTbb;

  specCache* spec_cache =  init_global_specCache(status);
  CHECK_STATUS_VOID(*status);
  setArrayToZero(spec_inp, n_ener_inp);
  for (int ii = 0; ii < rel_profile->n_zones; ii++) {
    assert(returnSpec->n_ener==n_ener);

    xillver_out[ii] = malloc(sizeof(double *) * n_ener);
    CHECK_MALLOC_VOID_STATUS(xillver_out[ii],status)

    xill_param->kTbb=Tin*xill_param->shiftTmaxRRet;  // currently set for testing
    if (noXillverReflection){
      getZoneIncidentReturnFlux(xill_param, returnSpec, xillver_out[ii], ii);
    } else {
      getZoneReflectedReturnFluxDiskframe(xill_param, rel_profile, returnSpec, xillver_out[ii], ii, status);
    }

    xillver_prim_out[ii] = scaledXillverPrimaryBBodyHighener(xill_param->kTbb, returnSpec->specRet[ii],
                                                            returnSpec->ener, returnSpec->n_ener, status);

    spec_conv_out[ii] = malloc(sizeof(double *) * n_ener);
    CHECK_MALLOC_VOID_STATUS(spec_conv_out[ii],status)
    convolveSpectrumFFTNormalized(ener, xillver_out[ii], rel_profile->flux[ii], spec_conv_out[ii], n_ener,
        1, 1, ii, spec_cache, status);


    rebin_spectrum(ener_inp, single_spec_inp, n_ener_inp, ener, spec_conv_out[ii], n_ener);

    for (int jj = 0; jj < n_ener_inp; jj++) {
      spec_inp[jj] += single_spec_inp[jj];
    }

  }

  // clean spectrum
  double lowestModelValue = 1e-8;
  setLowValuesToZero(spec_inp, n_ener_inp, lowestModelValue);
  setValuesOutsideToZero(spec_inp, ener_inp, n_ener_inp);

  // reset Tin parameter to be safe
  xill_param->kTbb = Tin;

  if (is_debug_run()) {
    char *fname = "!debug-testrr-bbody-obs-reflect.fits";
    if (fabs(xill_param->refl_frac) < 1e-8) {
      fname = "!debug-testrr-bbody-obs-primary.fits";
    }
    if (noXillverReflection){
      fname = "!debug-testrr-bbody-obs-mirror.fits";
    }
    fits_rr_write_2Dspec(fname, spec_conv_out, ener, n_ener,
                         returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, NULL, status);


    fits_rr_write_2Dspec("!debug-testrr-bbody-rframe-xillverRefl.fits", xillver_out, ener, n_ener,
                         returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, NULL, status);

    fits_rr_write_2Dspec("!debug-testrr-bbody-rframe-xillverPrim.fits", xillver_prim_out, ener, n_ener,
                         returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, NULL, status);


    fits_rr_write_2Dspec("!debug-testrr-bbody-rframe-specRet.fits", returnSpec->specRet, ener, n_ener,
                         returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, NULL, status);
    fits_rr_write_2Dspec("!debug-testrr-bbody-rframe-specPri.fits", returnSpec->specPri, ener, n_ener,
                         returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, NULL, status);

  }
  free_2d(&spec_conv_out, returnSpec->nrad);
  free_2d(&xillver_prim_out, returnSpec->nrad);

}

void set_std_param_relxill_bbret(double *inp_par) {
  inp_par[0] = 0.998;
  inp_par[1] = 60.0;
  inp_par[2] = -1.0;
  inp_par[3] = 1000.;
  inp_par[4] = 0.0;   // redshift
  inp_par[5] = 1.0;   // kTbb
  inp_par[6] = 2.0;   // logxi
  inp_par[7] = 1.0;   // Afe
  inp_par[8] = 15.0; // logN
  inp_par[9] = 1.0;   // refl_frac
  inp_par[10] = 1;   // fixReflFrac
  inp_par[11] = 1.4;   // shiftTmaxRRad
}



/** shift the spectrum such that we can calculate the line for 1 keV **/
double *shift_energ_spec_1keV(const double *ener, const int n_ener, double line_energ, double z, int *status) {

  double *ener1keV = (double *) malloc((n_ener + 1) * sizeof(double));
  CHECK_MALLOC_RET_STATUS(ener1keV, status, NULL)

  int ii;
  for (ii = 0; ii <= n_ener; ii++) {
    ener1keV[ii] = ener[ii] * (1 + z) / line_energ;
  }
  return ener1keV;
}

/** RELXILL MODEL FUNCTION for the BB returning radiation **/
// DEPRECATED!!!!!!!!!!!!!
void tdrelxillbbret(const double *ener0,
                    const int n_ener0,
                    double *photar,
                    const double *parameter,
                    const int n_parameter,
                    int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

 // init_par_relxill_bbret(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  double *ener = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);

  relxill_bb_kernel(ener, photar, n_ener0, xill_param, rel_param, 0, status);
  CHECK_STATUS_VOID(*status);

  free(ener);
  free(xill_param);
  free(rel_param);

}

