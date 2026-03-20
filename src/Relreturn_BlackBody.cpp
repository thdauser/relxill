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
#include "Rebin.h"
#include "Relbase.h"
#include "Relphysics.h"
#include "Relreturn_Datastruct.h"

extern "C" {
#include "xilltable_reading.h"
#include "common.h"
#include "relutility.h"
}

#define LIM_GFAC_RR_BBODY 0.001 // difference between gmin and gmax, above which the energy shift is taken into account



// Function definitions
static VecD_2D get_returnrad_specs(double* ener_inp,
                                   int nener_inp,
                                   returningFractions* dat,
                                   const VecD& temperature,
                                   int* status,
                                   int norder = 1);
static VecD_2D get_bbody_specs_binned(double* ener, int nener, const VecD &temperature, const int* status);

void bin_integrate_flux(int nener, const double* ener, VecD_2D& spec);

// Program Code


void fits_rr_write_2Dspec(const char* fname, double** spec_arr, const double* ener, int nener,
                          const VecD& rlo, const VecD& rhi, int nrad, returningFractions* dat, int* status)
{
    CHECK_STATUS_VOID(*status);

    // open the fits file
    fitsfile* fptr;
    if (fits_create_file(&fptr, fname, status))
    {
        relxill_check_fits_error(status);
        printf("   creating file %s failed\n", fname);
        CHECK_STATUS_VOID(*status);
    }


    int n1 = nrad;
    int n2 = nener;
    char dim1[10];
    char dim2[10];
    sprintf(dim1, "%iD", (int)1);
    sprintf(dim2, "%iD", (int)n2);

    if (dat != nullptr)
    {
        const int tfields = 6;
        char* ttype[] = {(char*)"rlo", (char*)"rhi", (char*)"ener", (char*)"spec", (char*)"fraci", (char*)"fret"};
        char* tform[] = {dim1, dim1, dim2, dim2, (char*)"50D", dim1};
        fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, nullptr, "2Dspec", status);
    }
    else
    {
        const int tfields = 4;
        char* ttype[] = {(char*)"rlo", (char*)"rhi", (char*)"ener", (char*)"spec"};
        char* tform[] = {dim1, dim1, dim2, dim2};
        fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, nullptr, "2Dspec", status);
    }

    relxill_check_fits_error(status);
    CHECK_STATUS_VOID(*status);

    int firstrow = 1; /* first row in table to write   */
    int firstelem = 1; /* first element in row  */

    // we use "const_cast" as with legacy cfitsio this is the only way (and assume those values are never changed)
    fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, n1, const_cast<double*>(rlo.data()), status);
    fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, n1, const_cast<double*>(rhi.data()), status);
    for (int ii = 0; ii < n1; ii++)
    {
        fits_write_col(fptr, TDOUBLE, 3, firstrow + ii, firstelem, n2, const_cast<double*>(ener), status);
        fits_write_col(fptr, TDOUBLE, 4, firstrow + ii, firstelem, n2, spec_arr[ii], status);
    }

    if (dat != nullptr)
    {
        for (int ii = 0; ii < n1; ii++)
        {
            fits_write_col(fptr, TDOUBLE, 5, firstrow + ii, firstelem, n1, dat->tf_r[ii].data(), status);
        }
        fits_write_col(fptr, TDOUBLE, 6, firstrow, firstelem, n1, dat->tabData->f_ret, status);
        if (fits_write_key(fptr, TDOUBLE, "SPIN", &(dat->a), nullptr, status))
        {
        }
    }


    relxill_check_fits_error(status);

    if (fptr != nullptr) { fits_close_file(fptr, status); }
}


void fits_rr_write_2Dspec(const char* fname, const std::vector<std::vector<double>>& spec_arr,
                          const double* ener, int nener,
                          const VecD& rlo, const VecD& rhi, int nrad, returningFractions* dat, int* status)
{
    // Convert to double** temporarily
    std::vector<double*> ptrs(spec_arr.size());
    for (size_t i = 0; i < spec_arr.size(); i++)
    {
        ptrs[i] = const_cast<double*>(spec_arr[i].data());
    }

    // Call the original implementation
    fits_rr_write_2Dspec(fname, ptrs.data(), ener, nener, rlo, rhi, nrad, dat, status);
}

VecD getTemperatureProfileDiskZones(returningFractions* dat, double Rin, double Tin, int* status)
{
    return get_tprofile(dat->rlo, dat->rhi, dat->nrad, Rin, Tin, TPROFILE_ALPHA, status);
}


/**
 * Computes a 2D spectral array based on a given blackbody temperature profile
 * and energy bins, normalizing the flux for radial zones.
 *
 * @param ener Pointer to an array containing the energy grid in keV.
 * @param nener The number of energy grid points, including boundaries.
 * @param temperature Pointer to an array of temperatures (in keV) for each radial zone.
 * @param status Pointer to an integer representing the status, used to handle errors.
 *               If the status is non-zero, the computation is not performed.
 *
 * @return A 2D vector with dimensions [dat->nrad][nener],
 *         where each row corresponds to an energy spectrum for a radial zone.
 *         The spectra are normalized to photons/s/keV/cm^2 using the provided energy grid.
 *         Returns empty VecD_2D if an error is encountered (status becomes non-zero).
 *
 * @note For the bin-intergrated version (as needed, e.g., by Xspec, see "get_bbody_spec_binned()")
 * @note The spectral values are set to zero outside the predefined energy band
 *       [EMIN_XILLVER, EMAX_XILLVER].
 * @note Memory for the 2D spectral array is allocated internally. The caller is
 *       responsible for releasing this memory using the appropriate deallocation function.
 */
static VecD_2D get_bbody_specs(double* ener, int nener, const VecD &temperature, const int* status)
{
    CHECK_STATUS_RET(*status, VecD_2D());

    const int nrad = static_cast<int>(temperature.size());
    VecD_2D spec_array(nrad, VecD(nener));

    for (int ii = 0; ii < nrad; ii++)
    {
        bbody_spec(ener, nener, spec_array[ii].data(), temperature[ii], 1.0);
        for (int jj = 0; jj < nener; jj++)
        {
            if (ener[jj] < EMIN_XILLVER || ener[jj + 1] > EMAX_XILLVER)
            {
                spec_array[ii][jj] = 0.0;
            }
        }
    }

    return spec_array;
}

/**
 * @brief Compute the binned blackbody return spectrum.
 *
 * Wrapper around `get_bbody_specs()` that additionally bin-integrates the
 * result on the provided energy grid.
 *
 * @param ener Energy bin boundaries.
 * @param nener Number of energy bins.
 * @param temperature Disk temperature profile.
 * @param status Status flag.
 *
 * @return Bin-integrated 2D spectrum, one row per radial zone.
 */
static VecD_2D get_bbody_specs_binned(double* ener, int nener, const VecD& temperature, const int* status) {

    auto spec_array =  get_bbody_specs(ener, nener, temperature, status);
    bin_integrate_flux(nener, ener, spec_array);
    return spec_array;
}


/**
 * @brief Calculates a 2D blackbody return radiation spectrum for relativistic disk models
 *
 * This function computes the returning radiation spectrum from a relativistic accretion disk
 * using blackbody emission. The spectrum is calculated across multiple radial zones and
 * energy bins, taking into account relativistic effects near a spinning black hole.
 *
 * @param ener Pointer to array of energy values [keV] (input energy grid)
 * @param nener Number of energy bins in the energy array
 * @param Tin Temperature at the inner disk radius [K]
 * @param Rin Inner disk radius [gravitational radii, r_g = GM/c²]
 * @param Rout Outer disk radius [gravitational radii]
 * @param spin Dimensionless black hole spin parameter (a* = a/M, range: -1 to 1)
 * @param status Pointer to status variable (0 = success, non-zero = error)
 *
 * @return Pointer to returnSpec2D structure containing:
 *         - 2D spectral array [n_radial_zones][n_energy_bins]
 *         - Radial grid information
 *         - Energy grid information
 *         Returns NULL on error
 *
 * @note The function assumes a standard temperature profile for the alpha diks with zero torque boundary
 * (see Poutanen+2007)
 * @note Memory for the returned structure must be freed by the caller
 * @note The spin parameter follows the convention where positive values indicate prograde rotation
 *
 * @see returnSpec2D for the structure definition
 * @see spec_diskbb for single-zone blackbody spectrum calculation
 */
returnSpec2D* spec_returnrad_blackbody(double* ener, int nener,
                                       double Tin, double Rin, double Rout, double spin, int* status)
{
    CHECK_STATUS_RET(*status, NULL);

    // 1 - get the fractions from the table (plus interpolation to current parameters)
    returningFractions* dat = get_rrad_fractions(spin, Rin, Rout, status);



    // 2 - temperature profile of the whole disk
    VecD temperature = getTemperatureProfileDiskZones(dat, Rin, Tin, status);

    // 3 - get primary and returning 2D spectrum (bin integrated spectra)
    VecD_2D spec_prim_zones = get_bbody_specs_binned(ener, nener, temperature, status);
    VecD_2D spec_rr_zones = get_returnrad_specs(ener, nener, dat, temperature, status, 1);
    VecD_2D spec_rr_zones_2ndOrder = get_returnrad_specs(ener, nener, dat, temperature, status, 2);

    if (is_debug_run())
    {
        for (int ii=0; ii<dat->nrad; ii++)
        {
            assert_double_array_isnormal(spec_rr_zones[ii].data(), nener);
            assert_double_array_isnormal(spec_prim_zones[ii].data(), nener);
        }
    }

    returnSpec2D* returnSpec = getReturnradOutputStructure(dat, spec_rr_zones, spec_prim_zones, ener, nener, status);
    delete_returningFractions(&dat);

    return returnSpec;
}


/*
 * Helper Routines for Black Body returning radiation calculations
 */


static void calc_rr_bbspec_gzone(double* ener, int nener, double* spec_ret, const VecD& spec_pri,
                                 double* gfac, int ng, const double* frac_g)
{
    // loop over all gfac and calculate the spectrum for each
    VecD spec_ret_g(nener, 0.0);
    setArrayToZero(spec_ret, nener);
    for (int kk = 0; kk < ng; kk++)
    {
        // we have to shift 1/gfac, as gfac=E_i/E_e and do not shift the spectrum, but shift the energy grid "backwards"
        // N_i(E_i) = N_e(E_e) = N_e(E_i / gfac)
        rebin_mean_flux_shifted(ener, nener, spec_ret_g, spec_pri, 1./gfac[kk]);
        for (int jj = 0; jj < nener; jj++)
        {
            spec_ret[jj] += spec_ret_g[jj] * frac_g[kk];
        }
    }
}


/**
 * @brief Compute the returning radiation blackbody spectrum for a single radial zone.
 *
 * Sums contributions from all radial zones that return photons to zone `irad`,
 * applying the appropriate energy shift (g-factor) and returning fractions.
 * The result is a spectrum in ph/s/keV/cm² (not bin-integrated).
 *
 * @param[out] spec     Output spectrum array (size nener), zeroed and filled on return.
 * @param      ener     Energy grid boundaries (size nener+1) [keV].
 * @param      nener    Number of energy bins.
 * @param      irad     Index of the receiving radial zone.
 * @param      spec_pri Primary blackbody spectra per zone (ph/s/keV/cm²).
 * @param      dat      Returning fractions data (g-factors, transfer fractions, f_ret).
 * @param      status   Error status flag.
 */
static void calc_rr_bbspec_ring(VecD& spec,
                                double* ener,
                                int nener,
                                int irad,
                                const VecD_2D& spec_pri,
                                returningFractions* dat,
                                const int* status)
{
    VecD spec_r(nener, 0.0);
    assert(spec.size() == nener);
    spec.assign(spec.size(), 0.0);

    for (int ii = 0; ii < dat->nrad; ii++)
    {
        auto gfac = get_gfac_grid(dat->tabData->gmin[irad][ii], dat->tabData->gmax[irad][ii], dat->tabData->ng);

        // only calculate it for a significant redshift, otherwise return BBODY without gshift
        if (fabs(dat->tabData->gmax[irad][ii] - dat->tabData->gmin[irad][ii]) > LIM_GFAC_RR_BBODY)
        {
            calc_rr_bbspec_gzone(ener, nener, spec_r.data(), spec_pri[ii], gfac.data(), dat->tabData->ng, dat->tabData->frac_g[irad][ii]);
        }
        else
        {
            // we have to shift 1/gfac, as gfac=E_i/E_e and do not shift the spectrum, but shift the energy grid "backwards"
            // N_i(E_i) = N_e(E_e) = N_e(E_i / gfac)
            rebin_mean_flux_shifted(ener, nener, spec_r, spec_pri[ii], 1./dat->tabData->gmin[irad][ii]);
        }

        // apply the correct fractions and add it to the zone output spectrum
        for (int jj = 0; jj < nener; jj++)
        {
            spec[jj] += spec_r[jj] * dat->tabData->f_ret[ii] * dat->tf_r[irad][ii];
        }
    }
    assert_double_array_isnormal(spec.data(), nener);
}


/**
 * @brief Compute returning radiation spectra from given primary spectra for all zones.
 *
 * Loops over all radial zones and calls `calc_rr_bbspec_ring` to apply the
 * returning fractions and energy shifts. Result is not bin-integrated.
 *
 * @param      dat      Returning fractions data (g-factors, transfer fractions, f_ret).
 * @param      status   Error status flag.
 * @param      nener    Number of energy bins.
 * @param      ener     Energy grid boundaries (size nener+1) [keV].
 * @param      spec_pri Primary spectra per zone (ph/s/keV/cm²).
 *
 * @return 2D spectrum [nrad][nener] in ph/s/keV/cm² (not bin-integrated).
 */
static VecD_2D calc_returnrad_from_input_spectra(returningFractions *dat, int *status, int nener, double *ener, const VecD_2D& spec_pri) {
    VecD_2D spec_zones(dat->nrad, VecD(nener));
    for (int ii = 0; ii < dat->nrad; ii++) {
        calc_rr_bbspec_ring(spec_zones[ii], ener, nener, ii, spec_pri, dat, status); // return ph / cm²/s/keV ( not bin integ.)
    }
    return spec_zones;
}


/**
 * @brief Compute bin-integrated returning radiation spectra for all radial zones.
 *
 * Builds blackbody primary spectra from the temperature profile on the internal
 * bbody energy grid, computes returning radiation via `calc_returnrad_from_input_spectra`,
 * then rebins onto the output grid and bin-integrates.
 *
 * @param      ener_inp    Output energy grid boundaries (size nener_inp+1) [keV].
 * @param      nener_inp   Number of output energy bins.
 * @param      dat         Returning fractions data (g-factors, transfer fractions, f_ret).
 * @param      temperature Disk temperature profile per radial zone [keV].
 * @param      status      Error status flag.
 * @param      norder     Number of returning radiation calculations (only this order is returned!)
 *
 * @return 2D spectrum [nrad][nener_inp] in cts/bin, or empty on error.
 */
static VecD_2D get_returnrad_specs(double* ener_inp,
                                   int nener_inp,
                                   returningFractions* dat,
                                   const VecD& temperature,
                                   int* status,
                                   int norder)
{
    CHECK_STATUS_RET(*status, VecD_2D());
    assert(!dat->tf_r[0].empty());

    int nener;
    double* ener;
    get_std_bbody_energy_grid(&nener, &ener, status);
    auto spec_pri = get_bbody_specs(ener, nener, temperature, status);

    std::vector<VecD_2D> spec_out(norder);
    // match the spectra to the output grid
    VecD_2D spec_sum_out(dat->nrad, VecD(nener));
    for (int iorder = 0; iorder < norder; iorder++) {
        spec_out[iorder] = calc_returnrad_from_input_spectra(dat, status, nener, ener,
            iorder == 0 ? spec_pri : spec_out[iorder - 1]);
        if (iorder == 0) {
            spec_sum_out = spec_out[iorder];
        } else {
            // Current implementation: do not add spectra, but give the last order only
            // add_2d(spec_sum_out, spec_out[iorder]);
            spec_sum_out =  spec_out[iorder];
        }
    }

    VecD_2D spec_zones_rebinned(dat->nrad, VecD(nener_inp));
    for (int ii = 0; ii < dat->nrad; ii++)
    {
        rebin_mean_flux(ener_inp, spec_zones_rebinned[ii].data(), nener_inp, ener, spec_sum_out[ii].data(), nener, status);
        assert_double_array_isnormal(spec_zones_rebinned[ii].data(),nener_inp);
    }
    bin_integrate_flux(nener_inp, ener_inp, spec_zones_rebinned);

    return spec_zones_rebinned;
}




/**
 * @brief Converts spectral values to counts per energy bin for multiple radial zones
 *
 * This function adjusts the spectral flux values for a given 2D spectra by integrating over
 * energy bins. For each radial zone, the spectral values are multiplied by the energy
 * bin width to produce counts per bin.
 *
 * @param nener Number of energy bins in the input energy array
 * @param ener Pointer to an array containing the boundaries of the energy bins [keV]
 *             (array size must be nener + 1)
 * @param spec 2D vector of spectral values, organized as spec[nrad][nener], containing
 *             the flux values to be converted (e.g., photons/cm²/s/keV).
 */
void bin_integrate_flux(int nener, const double* ener, VecD_2D& spec)
{
    const int nrad_spec = static_cast<int>(spec.size());
    for (int ii = 0; ii < nrad_spec; ii++)
    {
        for (int jj = 0; jj < nener; jj++)
        {
            spec[ii][jj] *= (ener[jj + 1] - ener[jj]); // now make it cts/bin
        }
    }
}

/**
 * diskbb spectrum in Xspec units [cts/bin]
 *
 * tested to give identical results to the Xspec implementation of diskbb
 * (deviation at <0.1keV are expected due to our outer radius only being 1000Rg)
 */
void spec_diskbb(double* ener, double* spec, int n, double Tin, double spin, int* status)
{
    CHECK_STATUS_VOID(*status);

    const double rin = kerr_rms(spin);
    returningFractions* dat = get_rrad_fractions(spin, rin, RMAX_RELRET, status);

    VecD temperature =
        get_tprofile(dat->rlo, dat->rhi, dat->nrad, 0.5 * (dat->rlo[0] + dat->rhi[0]), Tin, TPROFILE_DISKBB, status);

    VecD_2D spec_arr = get_bbody_specs_binned(ener, n, temperature, status);

    // need to multiply by the normal ring area to mimick a diskbb spectrum
    assert(dat->rhi[0] > dat->rlo[0]);
for (int jj = 0; jj < n; jj++)
    {
        for (int ii = 0; ii < dat->nrad; ii++)
        {

            spec_arr[ii][jj] *=
                (M_PI * (pow(dat->rhi[ii], 2) - pow(dat->rlo[ii], 2)));
        }
    }

    if (is_debug_run())
    {
        fits_rr_write_2Dspec("!debug-testrr-spec-diskbb.fits", spec_arr, ener, n, dat->rlo, dat->rhi, dat->nrad, dat,
                             status);
    }

    // Convert to double** temporarily for sum_2Dspec
    std::vector<double*> ptrs(dat->nrad);
    for (int ii = 0; ii < dat->nrad; ii++)
    {
        ptrs[ii] = spec_arr[ii].data();
    }
    sum_2Dspec(spec, ptrs.data(), n, dat->nrad, status);
}


std::vector<double> getRadialGridFromReturntab(returnSpec2D* spec)
{
    assert(spec != nullptr);
    assert(spec->nrad > 0);
    assert(!spec->rlo.empty());
    assert(!spec->rhi.empty());

    VecD rgrid;
    rgrid.reserve(spec->nrad + 1); // we use n+1 grid points

    // Add lower bounds
    for (int ii = 0; ii < spec->nrad; ii++)
    {
        rgrid.push_back(spec->rlo[ii]);
    }

    // Add the upper bound of the last zone
    rgrid.push_back(spec->rhi[spec->nrad - 1]);

    return rgrid;
}

static double getXillverNormFactorFromPrimarySpectrum(double* spec_in, const double* ener, int n_ener, int* status)
{
    CHECK_STATUS_RET(*status, 0.0);

    EnerGrid* egrid = get_coarse_xillver_energrid();

    auto xillverInputSpec = new double[egrid->nbins];
    CHECK_MALLOC_RET_STATUS(xillverInputSpec, status, 0.0)

    _rebin_spectrum(egrid->ener, xillverInputSpec, egrid->nbins, ener, spec_in, n_ener);

    // divide by the primary normalization factor to get the scaling of the xillver reflection spectrum
    double normFactorXill = calcXillverNormFromPrimarySpectrum(xillverInputSpec, egrid->ener, egrid->nbins, status);

    delete[] xillverInputSpec;

    return normFactorXill;
}


double normalizationFactorBBodyAtHighener(double kTbb, const double* spec1, const VecD& spec_bb, const double* ener, int n_ener)
{
    assert(spec1 != nullptr);
    assert(spec_bb.size() == n_ener);

    const double ELO_LIMIT_HIGHENER = 4 * kTbb;

    double normSpecIn = calcSumInEnergyBand(spec1, n_ener, ener, ELO_LIMIT_HIGHENER, EMAX_XILLVER_NORMALIZATION);
    double normBbodySpec = calcSumInEnergyBand(spec_bb.data(), n_ener, ener, ELO_LIMIT_HIGHENER, EMAX_XILLVER_NORMALIZATION);

    return normSpecIn / normBbodySpec;
}


double calcNormfacBBodyAtHighenergy(double kTbb, double* spec_in, const double* ener, int n_ener, int* status)
{
    assert(spec_in != nullptr);

    VecD spec_bb = getXillverPrimaryBBody(kTbb, spec_in, ener, n_ener, status);
    if (*status != EXIT_SUCCESS)
    {
        return 0.0;
    }

    double norm_fac = normalizationFactorBBodyAtHighener(kTbb, spec_in, spec_bb, ener, n_ener);

    assert( std::isnormal(norm_fac));

    return norm_fac;
}


/**
 * Computes the reflected return flux for a specific zone of a multi-zone accretion disk model.
 *
 * This function calculates the zone-specific reflected return flux in the disk frame
 * by calculating the xillver reflection spectrum, which is normalized by a normalization
 * factor derived from the primary spectrum, which is then adjusted by another factor at high energies
 * according to the black body which was actually used to calculate the tabulated xillver spectrum.
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
void calcZoneDirectAndReflectedReturnFlux(xillParam* xill_param, relline_spec_multizone* rel_profile,
                                         const returnSpec2D* returnSpec,
                                         VecD& xill_flux_returnrad, int izone, int* status)
{
    getAnglecorrXillSpec(xill_flux_returnrad, returnSpec->ener, returnSpec->n_ener, xill_param,
                         rel_profile->rel_cosne->dist[izone], status);

    double xillverReflectionNormFactor
        = getXillverNormFactorFromPrimarySpectrum(returnSpec->specRet[izone], returnSpec->ener, returnSpec->n_ener,
                                                  status);
    CHECK_STATUS_VOID(*status);

    double normfacMatchAtHighEnergies
        = calcNormfacBBodyAtHighenergy(xill_param->kTbb,
                                       returnSpec->specRet[izone],
                                       returnSpec->ener,
                                       returnSpec->n_ener,
                                       status);
    assert(std::isnormal(normfacMatchAtHighEnergies));
    CHECK_STATUS_VOID(*status);

    for (int jj = 0; jj < returnSpec->n_ener; jj++)
    {
        xill_flux_returnrad[jj] *= fabs(xill_param->boost);
        xill_flux_returnrad[jj] *= xillverReflectionNormFactor * normfacMatchAtHighEnergies;

        if (xill_param->boost >= 0)
        {
            xill_flux_returnrad[jj] += returnSpec->specPri[izone][jj];
        }
    }
}

void calcZoneDirectAndReturnFlux(xillParam* xill_param, const returnSpec2D* returnSpec, VecD& returnFlux, int ii)
{
    for (int jj = 0; jj < returnSpec->n_ener; jj++) {
        returnFlux[jj] = returnSpec->specRet[ii][jj] * fabs(xill_param->boost);
        if (xill_param->boost >= 0)
        {
            returnFlux[jj] += returnSpec->specPri[ii][jj];
        }
    }
}


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
VecD getXillverPrimaryBBody(double kTbb, double* spec_in, const double* ener, int n_ener, int* status)
{
    VecD spec_bb(n_ener);

    double xillverReflectionNormFactor = getXillverNormFactorFromPrimarySpectrum(spec_in, ener, n_ener, status);

    bbody_spec(ener, n_ener, spec_bb.data(), kTbb, 1.0);
    for (int jj = 0; jj < n_ener; jj++)
    {
        spec_bb[jj] *= (ener[jj + 1] - ener[jj]);
    }

    double bbodyNormFactor = getXillverNormFactorFromPrimarySpectrum(spec_bb.data(), ener, n_ener, status);

    for (int jj = 0; jj < n_ener; jj++)
    {
        spec_bb[jj] *= xillverReflectionNormFactor / bbodyNormFactor;
    }

    return spec_bb;
}


VecD getXillverPimaryBBodyNormalizedAtHighener(double kTbb, double* spec_in, const double* ener, int n_ener, int* status)
{
    VecD spec_bbody_xillver = getXillverPrimaryBBody(kTbb, spec_in, ener, n_ener, status);

    double normFac = normalizationFactorBBodyAtHighener(kTbb, spec_in, spec_bbody_xillver, ener, n_ener);

    for (int jj = 0; jj < n_ener; jj++)
    {
        spec_bbody_xillver[jj] *= normFac;
    }

    return spec_bbody_xillver;
}


static void setLowValuesToZero(double* spec, int n)
{
    double maxVal = 0.0;
    const double lowestModelValue = 1e-8;

    for (int ii = 0; ii < n; ii++)
    {
        maxVal = fmax(maxVal, spec[ii]);
    }

    for (int ii = 0; ii < n; ii++)
    {
        if (spec[ii] < maxVal * lowestModelValue)
        {
            spec[ii] = 0.0;
        }
    }
}

static void setValuesOutsideToZero(double* spec, const double* ener, int n)
{
    for (int ii = 0; ii < n; ii++)
    {
        if (ener[ii] < EMIN_XILLVER || ener[ii + 1] > EMAX_XILLVER)
        {
            spec[ii] = 0.0;
        }
    }
}


static int should_noXillverRefl_calculated()
{
    char* env = getenv("RELXILL_BBRET_NOREFL");

    if (env != nullptr && env[0] == '1')
    {
        return 1;
    }
    return 0;
}


void check_and_reset_value(double* ptr_val, double val_ref, const char* param_name)
{
    const double prec = 1e-6;
    if (fabs(*ptr_val-val_ref) > prec)
    {
        *ptr_val = val_ref;
        printf(" re-setting parameter %s to %.3e (only this value allowed)\n", param_name , val_ref);
    }
}

void check_parameter_bounds_relxillBB(xillParam* xill_param,relParam* rel_param)
{
    if (rel_param->model_type == MOD_TYPE_RELXILLBBRETXILL)
    {
        check_and_reset_value( &(rel_param->a), 0.998, "a");
        check_and_reset_value( &(rel_param->rin), kerr_rms(0.998), "Rin");
        check_and_reset_value( &(rel_param->rout), 1000, "Rout");
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
 * @param specOutput Pointer to the output spectrum array to store the computed spectrum. Any values written to this
 *        array will be set to zero
 * @param n_ener_inp Number of energy bins in the input energy grid.
 * @param xill_param_orig Struct pointer that holds parameters related to the XILLVER reflection model, blackbody
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
void relxill_bb_kernel(double* ener_inp, double* specOutput, int n_ener_inp, const xillParam* xill_param_orig, relParam* rel_param,
                       int* status)
{
    xillParam xill_param_struct = *xill_param_orig; // make a copy of the input parameters
    xillParam* xill_param = &xill_param_struct; // make a pointer to the copy

    check_parameter_bounds_relxillBB(xill_param, rel_param);

    CHECK_STATUS_VOID(*status);
    assert(xill_param->model_type == MOD_TYPE_RELXILLBBRET || xill_param->model_type == MOD_TYPE_RELXILLBBRETXILL);

    // get a standard grid for the convolution (is rebinned later to the input grid)
    auto ener_grid = get_relxill_conv_energy_grid();
    int n_ener = ener_grid->nbins;
    double* ener = ener_grid->ener;

    xillTable* xill_tab = nullptr;
    get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, status);

    returnSpec2D* returnSpec = spec_returnrad_blackbody(ener, n_ener, xill_param->kTbb,
                                                        rel_param->rin, rel_param->rout, rel_param->a,
                                                        status);
   CHECK_RELXILL_ERROR("calculating return rad spec, failed", status);

    auto radialGrid = getRadialGridFromReturntab(returnSpec);
    auto nrad_grid = returnSpec->nrad;
    assert(radialGrid.size()-1 == nrad_grid); // uses upper and lower bounds of the bins so nrad+1
    RelSysPar* systemParameters = get_system_parameters(rel_param, status);
    relline_spec_multizone
        * relProfile = relbase_profile(ener,
                                       n_ener,
                                       rel_param,
                                       systemParameters,
                                       xill_tab,
                                       radialGrid,
                                       nrad_grid,
                                       status);

    // ========== //
    VecD specSingleZone(n_ener_inp);
    VecD_2D specZonesConvOutput(returnSpec->nrad);
    VecD_2D specZonesEmittedCombined(returnSpec->nrad);
    VecD_2D specBbodyPrimZones(returnSpec->nrad);
    // ========== //

    double Tin = xill_param->kTbb;

    specCache* spec_cache = init_global_specCache(status);
    CHECK_STATUS_VOID(*status);
    setArrayToZero(specOutput, n_ener_inp);

    for (int izone = 0; izone < relProfile->n_zones; izone++)
    {
        assert(returnSpec->n_ener == n_ener);

        specZonesEmittedCombined[izone] = VecD(n_ener, 0.0);

        if (should_noXillverRefl_calculated())
        {
            // TODO: remove this once everything is working; only for testing
            calcZoneDirectAndReturnFlux(xill_param, returnSpec, specZonesEmittedCombined[izone], izone);
        }
        else
        {
            if (xill_param->model_type == MOD_TYPE_RELXILLBBRET)
            {
                xill_param->kTbb = Tin * xill_param->shiftTmaxRRet; // currently set for testing the xillverNS table
            }
            if (xill_param->model_type == MOD_TYPE_RELXILLBBRETXILL)
            {
                // very dangerous; need to check that it is within the table boundaries as well
                // (will only be used for the model with the xillverRTR table (MOD_TYPE_RELXILLBBRETXILL) )
                xill_param->radzone = (double)izone;
            }
            calcZoneDirectAndReflectedReturnFlux(xill_param, relProfile, returnSpec, specZonesEmittedCombined[izone],
                                                izone, status);
        }


        // calculate xillver primary spectrum (only used for additional/debug output),
        if (shouldOutfilesBeWritten())
        {
            specBbodyPrimZones[izone] =
                getXillverPimaryBBodyNormalizedAtHighener(xill_param->kTbb, returnSpec->specRet[izone],
                                                          returnSpec->ener,
                                                          returnSpec->n_ener, status);
        }

        specZonesConvOutput[izone] = VecD(n_ener, 0.0);
        convolveSpectrumFFTNormalized(ener,
                                      specZonesEmittedCombined[izone].data(),
                                      relProfile->flux[izone],
                                      specZonesConvOutput[izone].data(),
                                      n_ener,
                                      1,
                                      1,
                                      izone,
                                      spec_cache,
                                      status);


        _rebin_spectrum(ener_inp, specSingleZone.data(), n_ener_inp, ener, specZonesConvOutput[izone].data(), n_ener);

        // add spectra from the current zone to the output spectrum
        for (int jj = 0; jj < n_ener_inp; jj++)
        {
            specOutput[jj] += specSingleZone[jj];
        }
    }

    // clean spectrum
    setLowValuesToZero(specOutput, n_ener_inp);
    setValuesOutsideToZero(specOutput, ener_inp, n_ener_inp);

    // reset Tin parameter to be safe
    xill_param->kTbb = Tin;

    if (shouldOutfilesBeWritten())
    {
        std::string fname = "!undefined.fits";
        std::string fbase = "!__relxillOutput-rrad-bbody-";

        if (should_noXillverRefl_calculated())
        {
            if (fabs(xill_param->boost) < 1e-8)
            {
                fname = fbase + "obs-mirror-primary.fits";
            }
            else if (xill_param->boost < 0)
            {
                fname = fbase + "obs-mirror-refl.fits";
            }
            else
            {
                fname = fbase + "obs-mirror.fits";
            }
        }
        else
        {
            if (fabs(xill_param->boost) < 1e-8)
            {
                fname = fbase + "obs-primary.fits";
            }
            else if (xill_param->boost < 1)
            {
                fname = fbase + "obs-reflect.fits";
            }
            else
            {
                fname = fbase + "obs-total.fits";
            }
        }
        std::cout << " writing BBret diagnose outfile: " << fname.c_str() << std::endl;

        fits_rr_write_2Dspec(fname.c_str(), specZonesConvOutput, ener, n_ener,
                             returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);
        
        fname = fbase + "rframe-total.fits";
        fits_rr_write_2Dspec(fname.c_str(), specZonesEmittedCombined, ener, n_ener,
                             returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);

        fname = fbase + "rframe-bbodyRet.fits";
        fits_rr_write_2Dspec(fname.c_str(), returnSpec->specRet, ener, n_ener,
                             returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);

        fname = fbase + "rframe-bbodyEmit.fits";
        fits_rr_write_2Dspec(fname.c_str(), returnSpec->specPri, ener, n_ener,
                             returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);

    }

    delete_returnSpec2D(&returnSpec);
}

