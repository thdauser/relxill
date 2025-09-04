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

#include <iostream>
#include "Relreturn_BlackBody.h"
#include "Relbase.h"
#include "Relphysics.h"
#include "Relreturn_Datastruct.h"

extern "C" {
#include "xilltable.h"
#include "common.h"
#include "relutility.h"
}

#define LIM_GFAC_RR_BBODY 0.001 // difference between gmin and gmax, above which the energy shift is taken into account

// Function definitions
static double **get_returnrad_specs(double *ener_inp,
                                    int nener_inp,
                                    returningFractions *dat,
                                    const double *temperature,
                                    int *status);
static double **get_bbody_specs(double *ener, int nener, returningFractions *dat, double *temperature, int *status);

void normalizeFluxRrad(int nrad, int nener, const double *ener, double **spec);

// Program Code




void fits_rr_write_2Dspec(const char *fname, double **spec_arr, double *ener, int nener,
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

  if (dat!=nullptr) {
    const int tfields = 6;
    char *ttype[] = { (char*) "rlo", (char*) "rhi", (char*) "ener", (char*) "spec", (char*) "fraci", (char*) "fret"};
    char *tform[] = {dim1, dim1, dim2, dim2, (char*) "50D", dim1};
    fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, nullptr, "2Dspec", status);

  } else {
    const int tfields = 4;
    char *ttype[] = {(char*) "rlo", (char*) "rhi", (char*) "ener", (char*) "spec"};
    char *tform[] = {dim1, dim1, dim2, dim2};
    fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, nullptr, "2Dspec", status);
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

  if (dat!=nullptr) {
    for (int ii = 0; ii < n1; ii++) {
      fits_write_col(fptr, TDOUBLE, 5, firstrow + ii, firstelem, n1, dat->tf_r[ii], status);
    }
    fits_write_col(fptr, TDOUBLE, 6, firstrow, firstelem, n1, dat->tabData->f_ret, status);
    if( fits_write_key(fptr, TDOUBLE, "SPIN", &(dat->a),nullptr, status) ) {}
  }


  relxill_check_fits_error(status);

  if (fptr != nullptr) { fits_close_file(fptr, status); }

}


double* getTemperatureProfileDiskZones(returningFractions* dat, double Rin, double Tin, int* status){
  return get_tprofile(dat->rlo, dat->rhi, dat->nrad, Rin, Tin, TPROFILE_ALPHA, status);
}


/**
 * Computes the returning black body emission spectrum for the given radial zones.
 * The primary (emitted) and (incident) returning radiation are computed and returned
 * in a 2D structure.
 *
 * @param ener Array of energy values for the spectral grid.
 * @param spec Pointer to the array to receive the summed returning radiation spectrum.
 *           If nullptr, it is not computed.
 * @param spec_prim Pointer to the array to receive the summed primary radiation spectrum.
 *           If nullptr, it is not computed.
 * @param nener Number of energy bins in the spectral grid.
 * @param Tin Inner disk temperature in keV.
 * @param Rin Inner radius of the disk in gravitational radii.
 * @param Rout Outer radius of the disk in gravitational radii.
 * @param spin Dimensionless black hole spin parameter.
 * @param status Pointer to the integer status flag for error reporting and propagation.
 *           If an error occurs, the function sets this flag and may return a nullptr.
 * @return Pointer to a returnSpec2D structure that contains the computed primary and
 *         returning radiation spectra over the disk zones. If an error occurs,
 *         returns nullptr.
 *
 * Steps:
 * 1. Interpolates returning radiation fractions based on spin, Rin, Rout.
 * 2. Calculates the temperature profile across the disk zones using the provided
 *    inner disk temperature and inner radius.
 * 3. Computes the 2D spectral zones for the returning radiation and primary
 *    blackbody radiation.
 * 4. Constructs and returns the result structure, normalizing the spectra if needed.
 * 5. Optionally writes debugging FITS files if debugging mode is enabled.
 */
returnSpec2D *spec_returnrad_blackbody(double *ener, double *spec, double *spec_prim, int nener,
                                       double Tin, double Rin, double Rout, double spin, int *status) {

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
  if (spec!=nullptr) { // TODO: decide if we need this spec
    sum_2Dspec(spec, spec_rr_zones, nener, dat->nrad, status);
  }
  if (spec_prim!=nullptr) { // TODO: decide if we need this spec
    sum_2Dspec(spec_prim, spec_prim_zones, nener, dat->nrad, status);
  }

  if ( is_debug_run() ) {
    fits_rr_write_2Dspec("!debug-testrr-rframe-rr-bbody.fits", returnSpec->specRet, ener, nener,
                         dat->rlo, dat->rhi, dat->nrad, dat, status);
    fits_rr_write_2Dspec("!debug-testrr-rframe-prim-bbody.fits", returnSpec->specPri, ener, nener,
                         dat->rlo, dat->rhi, dat->nrad, dat, status);
  }

  free(temperature);

  return returnSpec;

}


/*
 * Helper Routines for Black Body returning radiation calculations
 */


static void calc_rr_bbspec_gzone(double *ener, int nener, double *spec, double temp,
                                 double *gfac, int ng, const double *frac_g) {
  // loop over all gfac and calculate the spectrum for each
  auto spec_g = new double[nener];
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

  auto gfac = new double[dat->tabData->ng];
  auto spec_r = new double[nener];
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
      spec[jj] += spec_r[jj] * dat->tabData->f_ret[ii] * dat->tf_r[irad][ii];
    }

  }
}



/**
 * Computes the returning radiation spectra for a given set of input photon energy grid,
 * considering the radial distribution of returning fractions and the provided temperature profile.
 *
 * @param ener_inp Pointer to the array of input photon energies (in keV).
 * @param nener_inp Number of input energy bins corresponding to `ener_inp`.
 * @param dat Pointer to the `returningFractions` structure, which holds radial properties,
 *            fractional distribution of radiation, and tabulated data for calculations.
 * @param temperature Pointer to the array of temperatures (in keV) defined per radial zone.
 * @param status Pointer to the status code variable. If non-zero, the function will return early with NULL.
 *
 * @return A dynamically allocated 2D array containing the computed spectra for all radial zones
 *         (units: cts/bin, standard for xspec). Each row in the array corresponds to a different
 *         radial zone, and each column represents a rebinned photon energy point. The function
 *         normalizes the flux over the radial zones after computation.
 *
 * The function operates as follows:
 * - Calls `get_std_bbody_energy_grid` to retrieve the standard blackbody energy grid configuration,
 *   which is used as the base energy grid.
 * - Allocates memory for the output 2D spectra array via `new_specZonesArr`.
 * - Iteratively calculates the returning radiation spectra for each radial zone using `calc_rr_bbspec_ring`.
 * - Rebinned spectra are computed for input energy bins using `rebin_mean_flux`.
 * - After iterating through all radial zones, the spectra data is normalized using `normalizeFluxRrad`
 *   to produce a final output in the desired format.
 *
 * The caller is responsible for deallocating the returned 2D array to avoid memory leaks.
 */
static double **get_returnrad_specs(double *ener_inp, int nener_inp, returningFractions *dat,
                                    const double *temperature, int *status) {
/* returns: - 2D-spectral array, unit is cts/bin  (xspec standard)
 *          - normalization
 * input:   - temperature profile (in keV)
 *          - interpolated return radiation fractions
 */


  CHECK_STATUS_RET(*status, NULL);

  int nener;
  double *ener;
  get_std_bbody_energy_grid(&nener, &ener, status);

  double **spec_zones = new_specZonesArr(nener_inp, dat->nrad, status);

  assert(dat->tf_r[0] != nullptr);

  auto spec = new double[nener];
  for (int ii = 0; ii < dat->nrad; ii++) {
    calc_rr_bbspec_ring(ener, spec, nener, ii, temperature, dat, status);  // return ph / cm²/s/keV ( not bin integ.)

    rebin_mean_flux(ener_inp, spec_zones[ii], nener_inp, ener, spec, nener, status);
  }

  normalizeFluxRrad(dat->nrad, nener_inp, ener_inp, spec_zones);


  return spec_zones;
}


/**
 * Computes a 2D spectral array based on a given blackbody temperature profile
 * and energy bins, normalizing the flux for radial zones.
 *
 * @param ener Pointer to an array containing the energy grid in keV.
 * @param nener The number of energy grid points, including boundaries.
 * @param dat Pointer to a returningFractions structure containing radial data
 *            and fractional photon information.
 * @param temperature Pointer to an array of temperatures (in keV) for each radial zone.
 * @param status Pointer to an integer representing the status, used to handle errors.
 *               If the status is non-zero, the computation is not performed.
 *
 * @return A dynamically allocated **2D spectral array** with dimensions [dat->nrad][nener],
 *         where each row corresponds to an energy spectrum for a radial zone.
 *         The spectra are normalized to counts/bin using the provided energy grid.
 *         Returns `NULL` if an error is encountered (status becomes non-zero).
 *
 * @note The spectral values are set to zero outside the predefined energy band
 *       [EMIN_XILLVER, EMAX_XILLVER].
 * @note Memory for the 2D spectral array is allocated internally. The caller is
 *       responsible for releasing this memory using the appropriate deallocation function.
 */
static double **get_bbody_specs(double *ener, int nener, returningFractions *dat, double *temperature, int *status) {

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



/*
 *  Input: bin integrated spectrum  [cts /bin/cm^2]
double getEmissivityNormFactor(returningFractions *dat, int nener, const double *ener, double **spec) {

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
 */

/**
 * diskbb spectrum in Xspec units [cts/bin]
 *
 * tested to give identical results to the Xspec implementation of diskbb
 * (deviation at <0.1keV are expected due to our outer radius only being 1000Rg)
 */
void spec_diskbb(double* ener, double* spec, int n, double Tin,  double spin, int* status) {

  CHECK_STATUS_VOID(*status);

  const double rin = kerr_rms(spin);
  returningFractions *dat = get_rrad_fractions(spin,  rin, RMAX_RELRET, status);

  double *temperature =
      get_tprofile(dat->rlo, dat->rhi, dat->nrad, 0.5 * (dat->rlo[0] + dat->rhi[0]), Tin, TPROFILE_DISKBB, status);

  double **spec_arr = get_bbody_specs(ener, n, dat, temperature, status);  // spectra integrated over the ring area

  /* need area correction for diskbb */
  assert(dat->rhi[0] > dat->rlo[0] );
  for (int jj = 0; jj < n; jj++) {
    for (int ii = 0; ii < dat->nrad; ii++) {
      // need to multiply by the normal ring area to mimick a diskbb spectrum
      spec_arr[ii][jj] *=
          (M_PI * (pow(dat->rhi[ii], 2) - pow(dat->rlo[ii], 2)));
    }
  }

  if (is_debug_run()) {
    fits_rr_write_2Dspec("!debug-testrr-spec-diskbb.fits", spec_arr, ener, n, dat->rlo, dat->rhi, dat->nrad, dat, status);
  }

  sum_2Dspec(spec, spec_arr, n, dat->nrad, status);

  free_2d(&spec_arr, dat->nrad);
  delete[]temperature;
}



double *getRadialGridFromReturntab(returnSpec2D *spec, int* status) {

  auto rgrid = new double[spec->nrad+1]; // we use n+1 grid points
  CHECK_MALLOC_RET_STATUS(rgrid, status, rgrid)

  for (int ii=0; ii<spec->nrad; ii++){
    rgrid[ii] = spec->rlo[ii];
  }
  rgrid[spec->nrad] = spec->rhi[spec->nrad-1];

  return rgrid;
}


static double getXillverNormFactorFromPrimarySpectrum(double* spec_in, double* ener, int n_ener, int* status){

  CHECK_STATUS_RET(*status,0.0);

  EnerGrid *egrid = get_coarse_xillver_energrid();

  auto xillverInputSpec = new double[egrid->nbins];
  CHECK_MALLOC_RET_STATUS(xillverInputSpec, status, 0.0)

  _rebin_spectrum(egrid->ener, xillverInputSpec, egrid->nbins, ener, spec_in, n_ener);

  // divide by the primary normalization factor to get the scaling of the xillver reflection spectrum
  double normFactorXill = calcXillverNormFromPrimarySpectrum(xillverInputSpec, egrid->ener, egrid->nbins, status);

  delete[] xillverInputSpec;

  return normFactorXill;
}


double normalizationFactorBBodyAtHighener(double kTbb, double *spec1, double *spec_bb, double *ener, int n_ener) {

  assert(spec1 != nullptr);
  assert(spec_bb != nullptr);

  const double ELO_LIMIT_HIGHENER = 4 * kTbb;

  double normSpecIn = calcSumInEnergyBand(spec1, n_ener, ener, ELO_LIMIT_HIGHENER, EMAX_XILLVER_NORMALIZATION);
  double normBbodySpec = calcSumInEnergyBand(spec_bb, n_ener, ener, ELO_LIMIT_HIGHENER, EMAX_XILLVER_NORMALIZATION);

  return normSpecIn / normBbodySpec;
}


double calcNormfacBBodyAtHighenergy(double kTbb, double *spec_in, double *ener, int n_ener, int *status) {

  assert(spec_in != nullptr);

  double *spec_bb = getXillverPrimaryBBody(kTbb, spec_in, ener, n_ener, status);
  if (*status != EXIT_SUCCESS) {
    return 0.0;
  }

  double norm_fac = normalizationFactorBBodyAtHighener(kTbb, spec_in, spec_bb, ener, n_ener);
  free(spec_bb);

  return norm_fac;
}


/**
 * Computes the reflected return flux for a specific zone of a multi-zone accretion disk model.
 *
 * This function calculates the zone-specific reflected return flux in the disk frame
 * by calculating the xillver reflection spectrum, which is normalized by a normalization
 * factor derived from the primary spectrum, which is then adjusted by another factor at high energies
 * according to the black body show which was actually used to calculate the loaded xillver spectrum.
 * The final flux is multiplied with the value of the of boost parameter (optionally add direct spectrum
 * if boost >=0).
 *
 * Assumptions:
 * - A single blackbody temperature (kTbb) is assumed for the reflection of all zones.
 *
 * Parameters:
 * @param xill_param          Pointer to the structure containing xillver model parameters (e.g., spectral parameters and disk settings).
 * @param rel_profile         Pointer to the multi-zone diskline profile, including relativistic correction data for all zones.
 * @param returnSpec          Pointer to the 2D return spectrum structure containing zone-specific primary and reflected flux data.
 * @param xill_flux_returnrad Array where the calculated reflected return flux for the specified zone is stored (length `returnSpec->n_ener`).
 * @param izone               Index of the zone for which the return flux is calculated.
 * @param status              Pointer to the status flag, used to report success or failure.
 */
void getZoneReflectedReturnFluxDiskframe(xillParam *xill_param, relline_spec_multizone* rel_profile, const returnSpec2D *returnSpec,
                                         double *xill_flux_returnrad, int izone, int* status) {

  getAnglecorrXillSpec(xill_flux_returnrad, returnSpec->ener, returnSpec->n_ener, xill_param,
                       rel_profile->rel_cosne->dist[izone], status);


  double xillverReflectionNormFactor
      = getXillverNormFactorFromPrimarySpectrum(returnSpec->specRet[izone], returnSpec->ener, returnSpec->n_ener, status);
  CHECK_STATUS_VOID(*status);

  double normfacMatchAtHighEnergies
      = calcNormfacBBodyAtHighenergy(xill_param->kTbb,
                                     returnSpec->specRet[izone],
                                     returnSpec->ener,
                                     returnSpec->n_ener,
                                     status);
  CHECK_STATUS_VOID(*status);

  for (int jj = 0; jj < returnSpec->n_ener; jj++) {
    xill_flux_returnrad[jj] *= fabs(xill_param->boost) ;
    xill_flux_returnrad[jj] *= xillverReflectionNormFactor * normfacMatchAtHighEnergies;

    if (xill_param->boost >= 0) {
      xill_flux_returnrad[jj] += returnSpec->specPri[izone][jj];
    }
  }

}

void getZoneIncidentReturnFlux(xillParam *xill_param, const returnSpec2D *returnSpec, double *returnFlux, int ii) {

  for (int jj = 0; jj < returnSpec->n_ener; jj++) {
    returnFlux[jj] = returnSpec->specRet[ii][jj] * fabs(xill_param->boost);
    if (xill_param->boost >= 0) {
      returnFlux[jj] += returnSpec->specPri[ii][jj];
    }
  }

}

/*
void getZoneDirectPrimaryFlux(xillParam *xill_param, const returnSpec2D *returnSpec, double *returnFlux, int ii) {

  for (int jj = 0; jj < returnSpec->n_ener; jj++) {
    returnFlux[jj] =  returnSpec->specPri[ii][jj];
  }

}*/



/**
 * @brief Generates a blackbody primary spectrum  of the given Xillver spectrum.
 *
 * This function calculates a blackbody spectrum and applies appropriate scaling
 * to normalize it based on the Xillver reflection normalization factor. It adjusts
 * the computed blackbody spectrum to match the normalization derived from the original
 * input spectrum. The resultant spectrum is scaled and returned to the caller.
 *
 * @param kTbb Temperature of the blackbody spectrum in keV.
 * @param spec_in Pointer to the array containing the input spectrum.
 * @param ener Pointer to the array containing the energy grid for the spectrum.
 * @param n_ener Number of elements in the energy grid and the spectrum arrays.
 * @param status Pointer to an integer used to track the success or failure of the operation.
 *
 * @return Pointer to an array containing the generated blackbody spectrum normalized
 *         to the Xillver primary spectrum. The caller assumes responsibility for deallocating
 *         the memory using delete[] when no longer needed. Returns nullptr if an error occurs.
 */
double* getXillverPrimaryBBody(double kTbb, double* spec_in, double* ener, int n_ener, int* status){

  auto spec_bb = new double[n_ener];
  CHECK_MALLOC_RET_STATUS(spec_bb, status, nullptr)

  double xillverReflectionNormFactor = getXillverNormFactorFromPrimarySpectrum(spec_in, ener, n_ener, status);

  bbody_spec(ener, n_ener, spec_bb, kTbb, 1.0);
  for (int jj = 0; jj < n_ener; jj++) {
    spec_bb[jj] *= (ener[jj + 1] - ener[jj]);
  }

  double bbodyNormFactor = getXillverNormFactorFromPrimarySpectrum(spec_bb, ener, n_ener, status);

  for (int jj = 0; jj < n_ener; jj++) {
    spec_bb[jj] *= xillverReflectionNormFactor / bbodyNormFactor;
  }

  // make sure memory is freed in case of an error
  if (*status != EXIT_SUCCESS) {
    delete[] spec_bb;
    return nullptr;
  }

  return spec_bb;
}


double* getXillverPimaryBBodyNormalizedAtHighener(double kTbb, double* spec_in, double* ener, int n_ener, int* status){

  double *spec_bbody_xillver = getXillverPrimaryBBody(kTbb, spec_in, ener, n_ener, status);
  if (*status != EXIT_SUCCESS || spec_bbody_xillver == nullptr) {
    return nullptr;
  }

  double normFac = normalizationFactorBBodyAtHighener(kTbb, spec_in, spec_bbody_xillver, ener, n_ener);

  for (int jj = 0; jj < n_ener; jj++) {
    spec_bbody_xillver[jj] *= normFac;
  }

  return spec_bbody_xillver;
}


static void setLowValuesToZero(double* spec, int n){

  double maxVal =  0.0;
  const double lowestModelValue = 1e-8;

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


static int should_noXillverRefl_calculated(){

  char* env = getenv("RELXILL_BBRET_NOREFL");

  if (  env!= nullptr &&  env[0]=='1'){
    return 1;
  } else {
    return 0;
  }

}


/**
 * @brief Computes the relativistic reflection spectrum using a blackbody as the primary radiation source.
 *
 * This function generates the relativistic convolution of the incident spectrum using the parameters of the
 * blackbody source and the relativistic disk. It combines contributions from multiple zones to produce the
 * full output spectrum.
 *
 * @param ener_inp Pointer to the input energy grid where the final output spectrum will be calculated.
 * @param specOutput Pointer to the output spectrum array to store the computed spectrum.
 * @param n_ener_inp Number of energy bins in the input energy grid.
 * @param xill_param Struct pointer that holds parameters related to the XILLVER reflection model, blackbody
 *        parameters, and source properties.
 * @param rel_param Struct pointer that holds parameters for the relativistic disk, emission, and system properties.
 * @param status Pointer to the status flag. Will be checked and updated during execution for error handling.
 *
 * The function performs the following major operations:
 * 1. Initialize and assert various parameters required for the computation.
 * 2. Generate a standard energy grid used for convolution.
 * 3. Load the XILLVER table data based on the model and primary source type.
 * 4. Compute the incident and reflected spectrum for each zone in the relativistic disk.
 * 5. Convolve the resulting spectrum using relativistic effects for each reflection zone.
 * 6. Rebin the convolved spectrum onto the input energy grid and combine contributions from all zones.
 * 7. Apply corrections to clean and format the output spectrum.
 * 8. Optionally write diagnostic FITS files for debugging purposes based on the control flags.
 *
 * Memory allocation and deallocation are handled for intermediate spectra within the function, ensuring proper cleanup
 * of resources to prevent memory leaks.
 *
 * Notes:
 * - The function assumes that the input and output spectra have been properly initialized by the caller.
 * - The input parameter xill_param->kTbb may be modified temporarily during computations but is restored
 *   to its original value upon completion.
 * - Diagnostic FITS files are written if diagnostic flags are set via shouldOutfilesBeWritten().
 */
void relxill_bb_kernel(double *ener_inp, double *specOutput, int n_ener_inp, xillParam *xill_param, relParam *rel_param,
                       int *status) {

  CHECK_STATUS_VOID(*status);
  assert(xill_param->model_type == MOD_TYPE_RELXILLBBRET);


  // get a standard grid for the convolution (is rebinned later to the input grid)
  auto ener_grid = get_relxill_conv_energy_grid();
  int n_ener = ener_grid->nbins;
  double *ener = ener_grid->ener;


  xillTable *xill_tab = nullptr;
  get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, status);

  returnSpec2D *returnSpec = spec_returnrad_blackbody(ener, nullptr, nullptr, n_ener, xill_param->kTbb, rel_param->rin,
                                                      rel_param->rout, rel_param->a, status);

  double *radialGrid = getRadialGridFromReturntab(returnSpec, status);
  RelSysPar *systemParameters = get_system_parameters(rel_param, status);
  relline_spec_multizone
      *relProfile = relbase_profile(ener,
                                    n_ener,
                                    rel_param,
                                    systemParameters,
                                    xill_tab,
                                    radialGrid,
                                    returnSpec->nrad,
                                    status);

  // ========== //
  auto specSingleZone = new double[n_ener_inp];
  auto specConvOutputZones = new double *[returnSpec->nrad];
  auto specXillverReflZones = new double *[returnSpec->nrad];
  auto specXillverPrimZones = new double *[returnSpec->nrad];
 // ========== //

  double Tin = xill_param ->kTbb;

  specCache* spec_cache =  init_global_specCache(status);
  CHECK_STATUS_VOID(*status);
  setArrayToZero(specOutput, n_ener_inp);

  for (int izone = 0; izone < relProfile->n_zones; izone++) {
    assert(returnSpec->n_ener==n_ener);

    specXillverReflZones[izone] = new double[n_ener];

    if ( should_noXillverRefl_calculated() ){
      getZoneIncidentReturnFlux(xill_param, returnSpec, specXillverReflZones[izone], izone);
    } else {
      xill_param->kTbb=Tin*xill_param->shiftTmaxRRet;  // currently set for testing
      getZoneReflectedReturnFluxDiskframe(xill_param,
                                          relProfile,
                                          returnSpec,
                                          specXillverReflZones[izone],
                                          izone,
                                          status);
    }

    // calculate xillver primary spectrum (only used for additional/debug output)
    specXillverPrimZones[izone] = getXillverPimaryBBodyNormalizedAtHighener(xill_param->kTbb, returnSpec->specRet[izone],
                                                                    returnSpec->ener, returnSpec->n_ener, status);
    //   calculated_combined_flux(specXillverReflZones[izone], specXillverPrimZones[izone], n_ener, xill_param->boost);

    specConvOutputZones[izone] = new double[n_ener];
    convolveSpectrumFFTNormalized(ener,
                                  specXillverReflZones[izone],
                                  relProfile->flux[izone],
                                  specConvOutputZones[izone],
                                  n_ener,
                                  1,
                                  1,
                                  izone,
                                  spec_cache,
                                  status);


    _rebin_spectrum(ener_inp, specSingleZone, n_ener_inp, ener, specConvOutputZones[izone], n_ener);

    // add spectra from the current zone to the output spectrum
    for (int jj = 0; jj < n_ener_inp; jj++) {
      specOutput[jj] += specSingleZone[jj];
    }

  }

  // clean spectrum
  setLowValuesToZero(specOutput, n_ener_inp);
  setValuesOutsideToZero(specOutput, ener_inp, n_ener_inp);

  // reset Tin parameter to be safe
  xill_param->kTbb = Tin;

  if ( shouldOutfilesBeWritten() ) {

  //  std::cout << " writing BBret diagnose outfiles " << std::endl;

    std::string fname = "!undefined.fits";

    if ( should_noXillverRefl_calculated() ){
      if (fabs(xill_param->boost) < 1e-8) {
        fname = "!debug-testrr-bbody-obs-mirror-primary.fits";
      } else if (xill_param->boost < 0 ){
        fname = "!debug-testrr-bbody-obs-mirror-refl.fits";
      } else {
        fname = "!debug-testrr-bbody-obs-mirror.fits";
      }

    } else {
      if (fabs(xill_param->boost) < 1e-8) {
        fname = "!debug-testrr-bbody-obs-primary.fits";
      } else if (xill_param->boost < 1) {
        fname = "!debug-testrr-bbody-obs-reflect.fits";
      } else {
        fname = "!debug-testrr-bbody-obs-total.fits";
      }
    }

    fits_rr_write_2Dspec(fname.c_str(), specConvOutputZones, ener, n_ener,
                         returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);


    fits_rr_write_2Dspec("!debug-testrr-bbody-rframe-xillverRefl.fits", specXillverReflZones, ener, n_ener,
                         returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);

    fits_rr_write_2Dspec("!debug-testrr-bbody-rframe-xillverPrim.fits", specXillverPrimZones, ener, n_ener,
                         returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);


    fits_rr_write_2Dspec("!debug-testrr-bbody-rframe-specRet.fits", returnSpec->specRet, ener, n_ener,
                         returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);
    fits_rr_write_2Dspec("!debug-testrr-bbody-rframe-specPri.fits", returnSpec->specPri, ener, n_ener,
                         returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);

  }
  free_2d(&specConvOutputZones, returnSpec->nrad);
  free_2d(&specXillverPrimZones, returnSpec->nrad);

}
