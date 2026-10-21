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

#include <algorithm>
#include <iostream>
#include <numeric>
#include "Relreturn_BlackBody.h"
#include "ModelDefinition.h"
#include "Rebin.h"
#include "Relbase.h"
#include "Relphysics.h"
#include "Relreturn_Datastruct.h"

#include "writeOutfiles.h"

extern "C" {
#include "xilltable_reading.h"
#include "common.h"
#include "relutility.h"
}

#define LIM_GFAC_RR_BBODY 0.001 // difference between gmin and gmax, above which the energy shift is taken into account
const std::string fbase_rrad_bbody = "!__relxillOutput-rrad-bbody-";  // filename base for output files
const std::string fbase_rrad_bbody_dat = "__relxillOutput-rrad-bbody-"; // filename base for output files


// Function definitions
static VecD_2D get_returnrad_specs(double* ener_inp,
                                   int nener_inp,
                                   const returningFractions *dat,
                                   const VecD& temperature,
                                   int* status,
                                   int norder = 1);

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
    sprintf(dim1, "%iD", 1);
    sprintf(dim2, "%iD", n2);

    if (dat != nullptr)
    {
        constexpr int tfields = 6;
        char *ttype[] = {
            const_cast<char *>("rlo"), const_cast<char *>("rhi"), const_cast<char *>("ener"),
            const_cast<char *>("spec"), const_cast<char *>("fraci"), const_cast<char *>("fret")
        };
        char *tform[] = {dim1, dim1, dim2, dim2, const_cast<char *>("50D"), dim1};
        fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, nullptr, "2Dspec", status);
    } else {
        constexpr int tfields = 4;
        char *ttype[] = {
            const_cast<char *>("rlo"), const_cast<char *>("rhi"),
            const_cast<char *>("ener"), const_cast<char *>("spec")};
        char* tform[] = {dim1, dim1, dim2, dim2};
        fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, nullptr, "2Dspec", status);
    }

    relxill_check_fits_error(status);
    CHECK_STATUS_VOID(*status);

    int firstrow = 1; /* first row in table to write   */
    int firstelem = 1; /* first element in a row */

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

// temperature profile of the disk emission, selected with the environment variable RELXILL_RRAD_BBODY_TPROFILE
static int get_tprofile_type_from_env() {
    switch (get_from_env_rrad_bbody_tprofile()) {
        case 2: return TPROFILE_DISKBB;
        case 3: return TPROFILE_NOVIKOV_THORNE;
        default: return TPROFILE_ALPHA; // this is for a value of 1 as well
    }
}

// temperature of each radial zone (evaluated at its mean radius), as used for the returning radiation
VecD get_temperature_profile_disk_zones(const returningFractions *dat, double Rin, double Tin, int *status) {
    return get_tprofile(dat->rlo, dat->rhi, dat->nrad, Rin, Tin, dat->a, get_tprofile_type_from_env(), status);
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
 * @note For the bin-integrated version (as needed, e.g., by Xspec, see "get_bbody_spec_binned()")
 * @note The spectral values are set to zero outside the predefined energy band
 *       [EMIN_XILLVER, EMAX_XILLVER].
 * @note Memory for the 2D spectral array is allocated internally. The caller is
 *       responsible for releasing this memory using the appropriate deallocation function.
 */
static VecD_2D get_bbody_specs(const double *ener, int nener, const VecD &temperature, const int* status)
{
    CHECK_STATUS_RET(*status, VecD_2D());

    const int nrad = static_cast<int>(temperature.size());
    VecD_2D spec_array(nrad, VecD(nener));

    for (int ii = 0; ii < nrad; ii++)
    {
        bbody_photon_flux_per_keV(ener, nener, spec_array[ii].data(), temperature[ii], 1.0);
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
 * @brief Direct (primary) disk emission of the given radial zones, calculated on a fine radial grid.
 *
 * Each zone [rlo, rhi] is divided into logarithmic sub-zones with r_hi/r_lo <= RRAD_BBODY_MAX_RATIO_SUBZONES. Each
 * sub-zone emits a blackbody with the temperature of the profile tprofile at its mean radius, and the spectrum of the
 * zone is the average of these blackbodies, weighted with the (flat-space) area of the sub-zones. This gives the
 * correct spectrum also for coarse zones, for which a single temperature is not a good approximation (e.g., the outer
 * zones of the returning-radiation table, see doc/returnrad_disk_emission.md).
 *
 * @param ener Energy bin boundaries [keV] (size nener+1).
 * @param nener Number of energy bins.
 * @param rlo Inner radii of the zones [r_g].
 * @param rhi Outer radii of the zones [r_g].
 * @param Rin Inner edge of the disk [r_g].
 * @param kT Temperature of the profile [keV]: maximal temperature (TPROFILE_ALPHA, TPROFILE_NOVIKOV_THORNE), or
 *           temperature at Rin (TPROFILE_DISKBB).
 * @param spin Spin (only used for TPROFILE_NOVIKOV_THORNE).
 * @param tprofile Temperature profile (TPROFILE_ALPHA, TPROFILE_DISKBB, or TPROFILE_NOVIKOV_THORNE).
 * @param status Status flag.
 *
 * @return Bin-integrated spectra per unit area [nzones][nener], set to zero outside [EMIN_XILLVER, EMAX_XILLVER].
 */
VecD_2D get_disk_emission_zones_binned(const double *ener,
                                       int nener,
                                       const VecD &rlo,
                                       const VecD &rhi,
                                       double Rin,
                                       double kT,
                                       double spin,
                                       int tprofile,
                                       int *status) {
    CHECK_STATUS_RET(*status, VecD_2D());
    const int nzones = static_cast<int>(rlo.size());
    assert(nzones > 0 && rhi.size() == rlo.size());

    // logarithmic sub-zones of each zone
    VecD rlo_sub;
    VecD rhi_sub;
    std::vector<int> izone_sub;
    for (int ii = 0; ii < nzones; ii++) {
        const int nsub =
                std::max(1,
                         static_cast<int>(ceil(log(rhi[ii] / rlo[ii]) / log(RRAD_BBODY_MAX_RATIO_SUBZONES) - 1e-10)));
        for (int kk = 0; kk < nsub; kk++) {
            rlo_sub.push_back(rlo[ii] * pow(rhi[ii] / rlo[ii], static_cast<double>(kk) / nsub));
            rhi_sub.push_back(rlo[ii] * pow(rhi[ii] / rlo[ii], static_cast<double>(kk + 1) / nsub));
            izone_sub.push_back(ii);
        }
    }
    const int nsub_total = static_cast<int>(rlo_sub.size());

    VecD temperature_sub = get_tprofile(rlo_sub, rhi_sub, nsub_total, Rin, kT, spin, tprofile, status);
    CHECK_STATUS_RET(*status, VecD_2D());

    // sum of the blackbody spectra of the sub-zones of each zone, weighted with their area; the blackbody is only
    // calculated up to E = 60 kT (above, it is below exp(-60) of its maximum)
    VecD_2D spec_zones(nzones, VecD(nener, 0.0));
    VecD area_zones(nzones, 0.0);
    VecD spec_sub(nener);
    for (int kk = 0; kk < nsub_total; kk++) {
        const int nener_bb = static_cast<int>(
            std::upper_bound(ener, ener + nener, 60 * temperature_sub[kk]) - ener); // bins below 60 kT
        bbody_photon_flux_per_keV(ener, std::min(nener_bb, nener), spec_sub.data(), temperature_sub[kk], 1.0);

        const double area = M_PI * (rhi_sub[kk] * rhi_sub[kk] - rlo_sub[kk] * rlo_sub[kk]);
        area_zones[izone_sub[kk]] += area;
        for (int jj = 0; jj < std::min(nener_bb, nener); jj++) {
            spec_zones[izone_sub[kk]][jj] += area * spec_sub[jj];
        }
    }

    // average over the area of each zone, bin-integrated, and zero outside [EMIN_XILLVER, EMAX_XILLVER]
    for (int ii = 0; ii < nzones; ii++) {
        for (int jj = 0; jj < nener; jj++) {
            if (ener[jj] < EMIN_XILLVER || ener[jj + 1] > EMAX_XILLVER) {
                spec_zones[ii][jj] = 0.0;
            } else {
                spec_zones[ii][jj] *= (ener[jj + 1] - ener[jj]) / area_zones[ii];
            }
        }
    }
    return spec_zones;
}


/**
 * @brief Radial range of each zone as covered by the relativistic line profiles.
 *
 * The line profiles of the zones (calc_relline_profile) are integrated over N_FRAD radial points from rin to rout,
 * uniform in 1/sqrt(r) (get_fine_radial_grid). Each point represents the cell between the midpoints to its neighbors
 * and is assigned completely to the zone in which it is located. The area of a zone in its line profile is therefore
 * the union of the cells of its points, which differs from [rlo, rhi] by up to half a cell at each boundary (up to
 * 4% of the area for the outer zones of the returning-radiation table). Calculating the direct disk emission of each
 * zone on this range keeps its area consistent with the line profile, such that no area is assigned to the spectrum
 * of a neighboring zone (see doc/returnrad_disk_emission.md).
 *
 * @param rlo Inner radii of the zones [r_g].
 * @param rhi Outer radii of the zones [r_g].
 * @param rin Inner edge of the disk, as used for the line profiles [r_g].
 * @param rout Outer edge of the disk, as used for the line profiles [r_g].
 * @param range_lo Output: inner radius of the range of each zone [r_g].
 * @param range_hi Output: outer radius of the range of each zone [r_g].
 *
 * @note A zone without any radial point (not the case for the zones of the returning-radiation table) keeps its
 * nominal range.
 */
void get_zone_ranges_relline_profile(const VecD &rlo,
                                     const VecD &rhi,
                                     double rin,
                                     double rout,
                                     VecD &range_lo,
                                     VecD &range_hi) {
    const int nzones = static_cast<int>(rlo.size());
    assert(nzones > 0 && rhi.size() == rlo.size());

    VecD rgrid(rlo);
    rgrid.push_back(rhi.back());

    VecD re(N_FRAD);
    get_fine_radial_grid(rin, rout, re.data(), N_FRAD); // descending, as for the line profiles

    range_lo.assign(nzones, 0.0);
    range_hi.assign(nzones, 0.0);
    std::vector<bool> has_points(nzones, false);
    for (int ii = 0; ii < N_FRAD; ii++) {
        const int izone = binary_search(rgrid.data(), nzones + 1, re[ii]); // as in calc_relline_profile
        const double r_cell_hi = (ii == 0) ? re[ii] : 0.5 * (re[ii - 1] + re[ii]);
        const double r_cell_lo = (ii == N_FRAD - 1) ? re[ii] : 0.5 * (re[ii] + re[ii + 1]);
        if (!has_points[izone]) {
            range_lo[izone] = r_cell_lo;
            range_hi[izone] = r_cell_hi;
            has_points[izone] = true;
        } else {
            range_lo[izone] = fmin(range_lo[izone], r_cell_lo);
            range_hi[izone] = fmax(range_hi[izone], r_cell_hi);
        }
    }
    for (int ii = 0; ii < nzones; ii++) {
        if (!has_points[ii]) {
            range_lo[ii] = rlo[ii];
            range_hi[ii] = rhi[ii];
        }
    }
}


/**
 * @brief Relativistic kernel of the direct disk emission, for a maximal temperature of 1 keV (kTbb = 1).
 *
 * Per unit area, the bin-integrated blackbody photon spectrum of the temperature T is T^3 times the one of 1 keV,
 * shifted in energy by the factor T. With the temperature profile T(r) = kTbb t(r), the direct emission of the disk is
 * therefore the convolution of the blackbody of kTbb with a kernel, which is the relativistic line profile of each
 * radial point at the line energy t(r), weighted with t(r)^3 (calc_relline_profile_radial_line_energy). The kernel is
 * integrated over the same N_FRAD radial points as the line profiles of the zones, but with the temperature of each
 * point instead of a single temperature per zone (see doc/returnrad_disk_emission.md).
 *
 * The kernel does not depend on kTbb. It is cached and only re-calculated if the relativistic parameters or the
 * temperature profile change.
 *
 * @param ener Energy grid of the convolution [keV] (size n_ener+1).
 * @param n_ener Number of energy bins.
 * @param rel_param Relativistic parameters.
 * @param sysPar System parameters of rel_param (radial points, transfer functions).
 * @param tprofile Temperature profile (TPROFILE_ALPHA, TPROFILE_DISKBB, or TPROFILE_NOVIKOV_THORNE).
 * @param status Status flag.
 *
 * @return Kernel [photons/bin] (size n_ener), to be convolved with the bin-integrated blackbody photon spectrum
 *         (bbody_photon_flux_per_keV) per unit area.
 */
static const VecD &get_disk_emission_kernel(const double *ener,
                                            int n_ener,
                                            const relParam *rel_param,
                                            RelSysPar *sysPar,
                                            int tprofile,
                                            int *status) {
    static relParam *cached_rel_param = nullptr;
    static int cached_tprofile = -1;
    static VecD cached_kernel;

    if (did_rel_param_change(cached_rel_param, rel_param) == 0 && tprofile == cached_tprofile
        && cached_kernel.size() == static_cast<size_t>(n_ener)) {
        return cached_kernel;
    }
    cached_tprofile = -1; // invalid until the new kernel is calculated

    // temperature of each radial point for kTbb = 1, with the radii in ascending order (TPROFILE_DISKBB needs the inner
    // edge as first radius); the innermost radial point is Rin up to rounding errors
    const int nr = sysPar->nr;
    VecD radius(nr);
    for (int ii = 0; ii < nr; ii++) {
        radius[ii] = fmax(sysPar->re[nr - 1 - ii], rel_param->rin);
    }
    const VecD temperature = get_tprofile(radius, radius, nr, rel_param->rin, 1.0, rel_param->a, tprofile, status);
    CHECK_STATUS_RET(*status, cached_kernel);

    VecD line_ener(nr);
    VecD weight(nr);
    for (int ii = 0; ii < nr; ii++) {
        line_ener[ii] = temperature[nr - 1 - ii]; // same (descending) order as sysPar->re
        weight[ii] = pow(line_ener[ii], 3);
    }

    cached_kernel.assign(n_ener, 0.0);
    calc_relline_profile_radial_line_energy(ener,
                                            n_ener,
                                            sysPar,
                                            line_ener.data(),
                                            weight.data(),
                                            cached_kernel.data(),
                                            status);
    CHECK_STATUS_RET(*status, cached_kernel);

    set_cached_rel_param(rel_param, &cached_rel_param, status);
    cached_tprofile = tprofile;
    return cached_kernel;
}


/**
 * @brief Observed direct disk emission of the relxillBB models [photons/bin], on the energy grid of the convolution.
 *
 * Convolution of the blackbody of the maximal temperature kTbb (bin-integrated photon spectrum per unit area) with the
 * kernel of the disk emission (get_disk_emission_kernel). The FFT convolution is cyclic, and the kernel shifts the
 * energies down by up to a factor ~100 (outer disk). The blackbody is therefore set to zero in the bins that the
 * kernel can shift outside the energy grid, as they would re-appear at the other end of the grid. This only changes
 * the spectrum observed below ~0.05 keV (~0.1 keV for the diskbb profile).
 */
static VecD calc_disk_emission_observed(double *ener,
                                        int n_ener,
                                        double kTbb,
                                        const relParam *rel_param,
                                        RelSysPar *sysPar,
                                        specCache *spec_cache,
                                        int *status) {
    VecD spec_obs(n_ener, 0.0);
    const VecD &kernel =
            get_disk_emission_kernel(ener, n_ener, rel_param, sysPar, get_tprofile_type_from_env(), status);
    CHECK_STATUS_RET(*status, spec_obs);

    // maximal shift of the kernel to lower and higher energies (in bins, 1 keV is no shift as in fftw_conv_spectrum)
    const auto is_positive = [](double val) { return val > 0; };
    const int ibin_1keV = binary_search(ener, n_ener + 1, 1.0);
    const int nbins_shift_down =
            ibin_1keV - static_cast<int>(std::find_if(kernel.begin(), kernel.end(), is_positive) - kernel.begin());
    const int nbins_shift_up =
            static_cast<int>(kernel.rend() - std::find_if(kernel.rbegin(), kernel.rend(), is_positive)) - 1 - ibin_1keV;

    VecD spec_bbody(n_ener);
    bbody_photon_flux_per_keV(ener, n_ener, spec_bbody.data(), kTbb, 1.0);
    for (int jj = 0; jj < n_ener; jj++) {
        if (jj < nbins_shift_down || jj >= n_ener - nbins_shift_up) {
            spec_bbody[jj] = 0.0;
        } else {
            spec_bbody[jj] *= ener[jj + 1] - ener[jj]; // bin-integrated
        }
    }

    // the FFT arrays of zone 0 are only used as buffer (the FFTs are always re-calculated in relxillbb_kernel)
    convolveSpectrumFFTNormalized(ener,
                                  spec_bbody.data(),
                                  kernel.data(),
                                  spec_obs.data(),
                                  n_ener,
                                  1,
                                  1,
                                  0,
                                  spec_cache,
                                  status);

    // convolveSpectrumFFTNormalized normalizes the number of photons in the band [EMIN_XILLVER, EMAX_XILLVER], but the
    // cool outer disk is observed with up to ~1% of its photons below EMIN_XILLVER: normalize it on the whole energy
    // grid instead, where the convolution conserves the number of photons (no bins are wrapped around)
    const double sum_obs = std::accumulate(spec_obs.begin(), spec_obs.end(), 0.0);
    const double norm_fac = std::accumulate(kernel.begin(), kernel.end(), 0.0)
                            * std::accumulate(spec_bbody.begin(), spec_bbody.end(), 0.0) / sum_obs;
    for (double &val: spec_obs) {
        val *= norm_fac;
    }
    return spec_obs;
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
 * @note The temperature profile is selected with RELXILL_RRAD_BBODY_TPROFILE (default: alpha disk with zero torque
 * at the inner edge, see doc/returnrad_disk_emission.md)
 * @note The primary (direct) disk emission of each zone (specPri) is calculated on a fine radial grid within the
 * radial range the zone covers in its relativistic line profile (see get_disk_emission_zones_binned and
 * get_zone_ranges_relline_profile), whereas the returning radiation (specRet) is calculated from a single blackbody
 * per zone. The observed direct emission of the relxillBB models does not use specPri, but is calculated for the
 * whole disk at once (calc_disk_emission_observed).
 * @note Memory for the returned structure must be freed by the caller
 * @note The spin parameter follows the convention where positive values indicate prograde rotation
 *
 * @see returnSpec2D for the structure definition
 * @see spec_multicolor_disk for the same disk emission without relativistic effects (used in the unit tests)
 */
returnSpec2D* spec_returnrad_blackbody(double* ener, int nener,
                                       double Tin,
                                       double Rin,
                                       double Rout,
                                       double spin,
                                       int *status) {
    CHECK_STATUS_RET(*status, NULL);

    // 1 - get the fractions from the table (plus interpolation to current parameters)
    returningFractions *dat = get_rrad_fractions(spin, Rin, Rout, status);


    // 2 - temperature profile of the whole disk (temperature of each zone, used for the returning radiation)
    VecD temperature = get_temperature_profile_disk_zones(dat, Rin, Tin, status);

    // 3 - get primary and returning 2D spectrum (bin integrated spectra): the primary (direct) disk emission of each
    //     zone is calculated on a fine radial grid within the radial range the zone covers in its relativistic line
    //     profile, the returning radiation from the zone temperatures
    VecD range_lo;
    VecD range_hi;
    get_zone_ranges_relline_profile(dat->rlo, dat->rhi, Rin, Rout, range_lo, range_hi);
    VecD_2D spec_prim_zones = get_disk_emission_zones_binned(ener, nener, range_lo, range_hi, Rin, Tin, dat->a,
                                                             get_tprofile_type_from_env(), status);

    int num_rrad_order = get_from_env_rrad_bbody_order();
  //  VecD_2D spec_rr_zones_primary = get_returnrad_specs(ener, nener, dat, temperature, status, 1);
    VecD_2D spec_rr_zones = get_returnrad_specs(ener, nener, dat, temperature, status, num_rrad_order);

    if (is_debug_run())
    {
        for (int ii=0; ii<dat->nrad; ii++)
        {
            assert_double_array_isnormal(spec_rr_zones[ii].data(), nener);
            assert_double_array_isnormal(spec_prim_zones[ii].data(), nener);
        }
    }

    returnSpec2D *returnSpec = get_returnrad_output_structure(dat, spec_rr_zones, spec_prim_zones, ener, nener, status);

    if (shouldOutfilesBeWritten()) {
        std::string fname = "!undefined.fits";
        fname = fbase_rrad_bbody + "rframe-bbodyEmit.fits";
        fits_rr_write_2Dspec(fname.c_str(), returnSpec->specPri, ener, nener,
                             returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);

        fname = fbase_rrad_bbody + "rframe-bbodyRet.fits";
        fits_rr_write_2Dspec(fname.c_str(), returnSpec->specRet, ener, nener,
                             returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);


        fname = fbase_rrad_bbody_dat + "rframe-bbodyTemp.dat";
        FILE *fp = fopen(fname.c_str(), "w+");
        fprintf(fp, "#zone \t rlo \t rhi \t T_col \n");
        for (int ii = 0; ii < returnSpec->nrad; ii++) {
            fprintf(fp,
                    "%i \t %e \t %e \t %e \n",
                    ii,
                    returnSpec->rlo[ii],
                    returnSpec->rhi[
                        +ii],
                    temperature[ii]);
        }
        fclose_errormsg(fp, fname.c_str());
    }

    delete_returningFractions(&dat);

    return returnSpec;
}


/*
 * Helper Routines for Black Body returning radiation calculations
 */


static void calc_rr_bbspec_gzone(double *ener,
                                 int nener,
                                 double *spec_ret,
                                 const VecD &spec_pri,
                                 const double *gfac,
                                 int ng,
                                 const double* frac_g)
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
 */
static void calc_rr_bbspec_ring(VecD& spec,
                                double* ener,
                                int nener,
                                int irad,
                                const VecD_2D &spec_pri,
                                const returningFractions* dat)
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
 * @brief Compute returning radiation spectra from the given primary spectra for all zones.
 *
 * Loops over all radial zones and calls `calc_rr_bbspec_ring` to apply the
 * returning fractions and energy shifts. The result is not bin-integrated.
 *
 * @param      dat      Returning fractions data (g-factors, transfer fractions, f_ret).
 * @param      nener    Number of energy bins.
 * @param      ener     Energy grid boundaries (size nener+1) [keV].
 * @param      spec_pri Primary spectra per zone (ph/s/keV/cm²).
 *
 * @return 2D spectrum [nrad][nener] in ph/s/keV/cm² (not bin-integrated).
 */
static VecD_2D calc_returnrad_from_input_spectra(const returningFractions *dat,
                                                 int nener,
                                                 double *ener,
                                                 const VecD_2D &spec_pri) {
    VecD_2D spec_zones(dat->nrad, VecD(nener));
    for (int ii = 0; ii < dat->nrad; ii++) {
        calc_rr_bbspec_ring(spec_zones[ii], ener, nener, ii, spec_pri, dat); // return ph / cm²/s/keV ( not bin integ.)
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
 * @param      norder     Number of returning radiation calculations (combined orders are returned!)
 *
 * @return 2D spectrum [nrad][nener_inp] in cts/bin, or empty on error.
 */
static VecD_2D get_returnrad_specs(double *ener_inp,
                                   int nener_inp,
                                   const returningFractions* dat,
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
    VecD_2D spec_sum_out(dat->nrad, VecD(nener));
    for (int iorder = 0; iorder < norder; iorder++) {
        spec_out[iorder] = calc_returnrad_from_input_spectra(dat, nener, ener, iorder == 0 ? spec_pri : spec_out[iorder - 1]);
        if (iorder == 0) {
            spec_sum_out = spec_out[iorder];
        } else {
            add_2d(spec_sum_out, spec_out[iorder]);
        }
    }

    // match the spectra to the output grid
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
void bin_integrate_flux(int nener, const double *ener, VecD_2D &spec) {
    const int nrad_spec = static_cast<int>(spec.size());
    for (int ii = 0; ii < nrad_spec; ii++) {
        for (int jj = 0; jj < nener; jj++) {
            spec[ii][jj] *= (ener[jj + 1] - ener[jj]); // now make it cts/bin
        }
    }
}

/**
 * multicolor disk spectrum [photons/bin, arbitrary normalization] for the given radial zones [rlo, rhi]
 *
 * The spectrum of each zone is calculated as the direct disk emission in the model (get_disk_emission_zones_binned,
 * i.e., on a fine radial grid within each zone), and is weighted with its flat-space ring area (i.e., without
 * relativistic effects):
 * - TPROFILE_ALPHA: alpha disk with zero torque at the inner edge rlo[0], with the maximal temperature kT
 * - TPROFILE_DISKBB: diskbb profile, with the temperature kT at the inner edge rlo[0]
 *
 * It uses the same temperature profiles and blackbody spectra as the disk emission of the relxillBB models, and is
 * used in the unit tests to compare them with the Xspec models ezdiskbb and diskbb, and to estimate the precision of
 * the radial zones (see doc/returnrad_disk_emission.md).
 */
void spec_multicolor_disk(const double *ener,
                          double *spec,
                          int n,
                          const VecD &rlo,
                          const VecD &rhi,
                          double kT,
                          int tprofile,
                          int *status) {
    CHECK_STATUS_VOID(*status);
    const int nrad = static_cast<int>(rlo.size());
    assert(nrad > 0 && rhi.size() == rlo.size());
    assert(tprofile != TPROFILE_NOVIKOV_THORNE); // relativistic profile (needs the spin), not for this function

    VecD_2D spec_arr = get_disk_emission_zones_binned(ener, n, rlo, rhi, rlo[0], kT, 0.0, tprofile, status);
    CHECK_STATUS_VOID(*status);

    for (int jj = 0; jj < n; jj++) {
        spec[jj] = 0.0;
        for (int ii = 0; ii < nrad; ii++) {
            spec[jj] += spec_arr[ii][jj] * M_PI * (rhi[ii] * rhi[ii] - rlo[ii] * rlo[ii]);
        }
    }
}


VecD get_radial_grid_from_returntab(const returnSpec2D *spec) {
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

static double get_xillver_norm_factor_from_primary_spectrum(const double *spec_in,
                                                            const double *ener,
                                                            int n_ener,
                                                            int *status) {
    CHECK_STATUS_RET(*status, 0.0);

    EnerGrid* egrid = get_coarse_xillver_energrid();

    auto xillverInputSpec = new double[egrid->nbins];
    CHECK_MALLOC_RET_STATUS(xillverInputSpec, status, 0.0)

    _rebin_spectrum(egrid->ener, xillverInputSpec, egrid->nbins, ener, spec_in, n_ener);

    // divide by the primary normalization factor to get the scaling of the xillver reflection spectrum
    double normFactorXill =
            calc_xillver_norm_from_primary_spectrum(xillverInputSpec, egrid->ener, egrid->nbins, status);

    delete[] xillverInputSpec;

    return normFactorXill;
}


/**
 * Computes the reflected return flux for a specific zone of a multi-zone accretion disk model.
 *
 * This function calculates the zone-specific reflected return flux in the disk frame
 * by calculating the xillver reflection spectrum, which is normalized by a normalization
 * factor derived from the primary (returning) spectrum.
 * The final flux is multiplied with the absolute value of the boost parameter. The direct disk emission is not
 * included, as it is calculated for the whole disk at once (see calc_disk_emission_observed).
 *
 * Assumptions:
 * - A single blackbody temperature (kTbb) is assumed for the reflection of all zones.
 *
 * @see @ref returnrad_normalizations for a derivation of the normalization
 *      factors applied below and the rationale for splitting them this way.
 *
 * Parameters:
 * @param xill_param          Pointer to the structure containing xillver model parameters (e.g., spectral parameters and disk settings).
 * @param rel_profile         Pointer to the multi-zone diskline profile, including relativistic correction data for all zones.
 * @param returnSpec          Pointer to the 2D return spectrum structure containing zone-specific primary and reflected flux data.
 * @param xill_flux_returnrad Array where the calculated reflected return flux for the specified zone is stored (length `returnSpec->n_ener`).
 * @param izone               Index of the zone for which the return flux is calculated.
 * @param status              Pointer to the status flag, used to report success or failure.
 */
void calc_zone_xillver_reflected_flux(const xillParam &xill_param,
                                      const relline_spec_multizone *rel_profile,
                                      const returnSpec2D *returnSpec,
                                      VecD &xill_flux_returnrad,
                                      int izone,
                                      int *status) {
    get_angle_weighted_xill_spec(xill_flux_returnrad,
                                 returnSpec->ener,
                                 returnSpec->n_ener, &xill_param,
                                 rel_profile->rel_cosne->dist[izone],
                                 status);

    double xillverReflectionNormFactor
            = get_xillver_norm_factor_from_primary_spectrum(returnSpec->specRet[izone],
                                                            returnSpec->ener,
                                                            returnSpec->n_ener,
                                                            status);
    CHECK_STATUS_VOID(*status);

    for (int jj = 0; jj < returnSpec->n_ener; jj++) {
        xill_flux_returnrad[jj] *= fabs(xill_param.boost);
        xill_flux_returnrad[jj] *= xillverReflectionNormFactor;
    }
}

/**
   * @brief Compute the per-zone flux of the returning radiation ("mirrored" instead of reflected).
   *
   * The returning radiation is multiplied with the absolute value of the boost parameter. The direct disk emission is
   * not included, as it is calculated for the whole disk at once (see calc_disk_emission_observed).
   *
   * @param xill_param   Xillver parameters; only @c boost is read here.
   * @param returnSpec   Precomputed return-radiation spectra
   * @param returnFlux   Output flux for zone ii, sized in returnSpec.
   * @param ii           Radial zone index.
   */
void calc_zone_mirrored_flux(const xillParam &xill_param,
                             const returnSpec2D *returnSpec,
                             VecD &returnFlux,
                             int ii) {
    for (int jj = 0; jj < returnSpec->n_ener; jj++) {
        returnFlux[jj] = returnSpec->specRet[ii][jj] * fabs(xill_param.boost);
    }
}


static void set_low_values_to_zero(double* spec, int n)
{
    double maxVal = 0.0;
    constexpr double lowestModelValue = 1e-8;

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

static void set_values_outside_to_zero(double* spec, const double* ener, int n)
{
    for (int ii = 0; ii < n; ii++)
    {
        if (ener[ii] < EMIN_XILLVER || ener[ii + 1] > EMAX_XILLVER) {
            spec[ii] = 0.0;
        }
    }
}


void check_and_reset_value(double* ptr_val, double val_ref, const char* param_name)
{
    constexpr double prec = 1e-6;
    if (fabs(*ptr_val-val_ref) > prec)
    {
        *ptr_val = val_ref;
        printf(" re-setting parameter %s to %.3e (only this value allowed)\n", param_name, val_ref);
    }
}

void check_parameter_bounds_relxillbb(xillParam* xill_param,relParam* rel_param)
{
    if (rel_param->model_type == MOD_TYPE_RELXILLBB) {
        check_and_reset_value(&rel_param->a, 0.998, "a");
        check_and_reset_value(&rel_param->rin, kerr_rms(0.998), "Rin");
        check_and_reset_value(&rel_param->rout, 1000, "Rout");
        check_and_reset_value(&xill_param->lxi, 2, "logxi");
        check_and_reset_value(&xill_param->dens, 18, "logN");
    }
}

IonGradient calc_ion_gradient_from_returnrad(
    const returnSpec2D *returnSpec,
    const ModelDefinition &params,
    int *status) {
    const auto radialGrid = get_radial_grid_from_returntab(returnSpec);
    IonGradient ion_gradient{radialGrid, ION_GRAD_TYPE_ALPHA, 0};
    const auto primary_source_params = PrimarySourceParameters{params};
    auto emis_profile =
            calculate_emis_profile_from_spectra(returnSpec->specRet,
                                                returnSpec->rlo,
                                                returnSpec->rhi,
                                                returnSpec->nrad,
                                                returnSpec->ener,
                                                returnSpec->n_ener,
                                                status);
    ion_gradient.calculate_gradient(*emis_profile, primary_source_params);
    free_emisProfile(emis_profile);

    if (shouldOutfilesBeWritten()) {
        ion_gradient.write_to_file((fbase_rrad_bbody_dat + "emis_profile.dat").c_str());
    }

    return ion_gradient;
}


// boost=0: only the direct (primary) emission, without reflection
static bool is_boost_primary_only(double boost) {
    return fabs(boost) < 1e-8;
}

// filename suffix of the diagnose outfiles, following the boost sign convention
// (boost<0: reflected only, boost=0: primary only, boost>0: primary + reflected)
static std::string get_bbret_diagnose_suffix(double boost) {
    std::string prefix = mirrored_reflection_instead_of_xillver() ? "mirror-" : "";
    if (is_boost_primary_only(boost)) {
        return prefix + "primary";
    }
    if (boost < 0) {
        return prefix + "reflect";
    }
    return prefix + "total";
}

static void write_bbret_diagnose_outfiles(const xillParam *xill_param,
                                          const VecD_2D &specZonesConvOutput,
                                          const VecD_2D &specZonesEmittedCombined,
                                          const double *ener,
                                          int n_ener,
                                          const returnSpec2D *returnSpec,
                                          int *status) {
    const std::string suffix = get_bbret_diagnose_suffix(xill_param->boost);

    std::string fname = fbase_rrad_bbody + "obs-" + suffix + ".fits";
    std::cout << " writing BBret diagnose outfile: " << fname.c_str() << std::endl;

    fits_rr_write_2Dspec(fname.c_str(),
                         specZonesConvOutput,
                         ener,
                         n_ener,
                         returnSpec->rlo,
                         returnSpec->rhi,
                         returnSpec->nrad,
                         nullptr,
                         status);

    // the spectra of the zones only contain the reflected emission (the direct emission is calculated for the whole
    // disk, see calc_disk_emission_observed), i.e., they are zero for boost=0; the direct emission of each zone in the
    // rest frame (specPri) is written to rframe-bbodyEmit.fits by spec_returnrad_blackbody()
    if (is_boost_primary_only(xill_param->boost)) {
        return;
    }

    fname = fbase_rrad_bbody + "rframe-" + suffix + ".fits";
    fits_rr_write_2Dspec(fname.c_str(),
                         specZonesEmittedCombined,
                         ener,
                         n_ener,
                         returnSpec->rlo,
                         returnSpec->rhi,
                         returnSpec->nrad,
                         nullptr,
                         status);
}

static void set_xill_param_from_iongrad(xillParam &xill_param_zone, const IonGradient &ion_gradient, int izone) {
    xill_param_zone.dens = ion_gradient.dens[izone];
    xill_param_zone.lxi = ion_gradient.lxi[izone];
    xill_param_zone.radzone = static_cast<double>(izone);
}

/**
 * @brief Computes the relativistic reflection spectrum using a blackbody as the primary radiation source.
 *
 * This function generates the relativistic convolution of the incident spectrum using the parameters of the
 * blackbody source and the relativistic disk. It combines contributions from multiple zones to produce the
 * full output spectrum.
 *
 * @param spectrum The input/output XspecSpectrum containing the energy grid and flux array.
 * @param params Model definition holding all model parameters (rel + xill are extracted internally).
 * @param status Pointer to the status flag. Will be checked and updated during execution for error handling.
 *
 * The function performs the following major operations:
 * 1. Initialize and assert various parameters required for the computation.
 * 2. Generate a standard energy grid used for convolution.
 * 3. Load the XILLVER table data based on the model and primary source type.
 * 4. Calculate the ionization gradient from the return radiation emissivity profile.
 * 5. Compute the observed direct disk emission (boost >= 0) for the whole disk with a single convolution
 *    (calc_disk_emission_observed).
 * 6. Compute the reflected spectrum for each zone in the relativistic disk.
 * 7. Convolve the resulting spectrum using relativistic effects for each reflection zone.
 * 8. Rebin the convolved spectrum onto the input energy grid and combine contributions from all zones.
 * 9. Apply corrections to clean and format the output spectrum.
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
void relxillbb_kernel(const XspecSpectrum &spectrum,
                      const ModelDefinition &params,
                      int* status)
{
    double *ener_inp = spectrum.energy;
    double *specOutput = spectrum.flux;
    int n_ener_inp = spectrum.num_flux_bins();

    auto rel_param = get_rel_params(params);
    auto xill_param = get_xill_params(params);

    check_parameter_bounds_relxillbb(xill_param, rel_param);

    CHECK_STATUS_VOID(*status);
    assert(xill_param->model_type == MOD_TYPE_RELXILLBB_NS || xill_param->model_type == MOD_TYPE_RELXILLBB);

    // get a standard grid for the convolution (is rebinned later to the input grid)
    auto ener_grid = get_relxill_conv_energy_grid();
    int n_ener = ener_grid->nbins;
    double* ener = ener_grid->ener;

    xillTable *xill_tab = nullptr;
    get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, status);

    returnSpec2D* returnSpec = spec_returnrad_blackbody(ener, n_ener, xill_param->kTbb,
                                                        rel_param->rin, rel_param->rout, rel_param->a,
                                                        status);
    CHECK_RELXILL_ERROR("calculating return rad spec, failed", status);
    auto ion_gradient = calc_ion_gradient_from_returnrad(returnSpec, params, status);


    auto radialGrid = get_radial_grid_from_returntab(returnSpec);
    RelSysPar* systemParameters = get_system_parameters(rel_param, status);
    relline_spec_multizone
        * relProfile = relbase_profile(ener,
                                       n_ener,
                                       rel_param,
                                       systemParameters,
                                       xill_tab,
                                       radialGrid,
                                       status);


    // ========== //
    VecD specSingleZone(n_ener_inp);
    VecD_2D specZonesConvOutput(returnSpec->nrad);
    VecD_2D specZonesEmittedCombined(returnSpec->nrad);
    // ========== //

    specCache *spec_cache = init_global_specCache(status);
    CHECK_STATUS_VOID(*status);
    setArrayToZero(specOutput, n_ener_inp);

    // direct disk emission (boost >= 0): a single convolution for the whole disk, with the temperature of each radial
    // point (the zones below only contain the reflected emission)
    if (xill_param->boost >= 0) {
        VecD specDirectObs = calc_disk_emission_observed(ener,
                                                         n_ener,
                                                         xill_param->kTbb,
                                                         rel_param,
                                                         systemParameters,
                                                         spec_cache,
                                                         status);
        _rebin_spectrum(ener_inp, specOutput, n_ener_inp, ener, specDirectObs.data(), n_ener);
    }

    for (int izone = 0; izone < relProfile->n_zones; izone++)
    {
        assert(returnSpec->n_ener == n_ener);

        // copy parameters for bug safety to not make changes to the global xill_param
        xillParam xill_param_zone = *xill_param;  // flat copy works, as no arrays contained

        specZonesEmittedCombined[izone] = VecD(n_ener, 0.0);
        specZonesConvOutput[izone] = VecD(n_ener, 0.0);

        // the zones only contain the reflected emission, i.e., nothing for boost=0 (the normalization of the FFT
        // convolution is not defined for a spectrum without flux)
        if (is_boost_primary_only(xill_param->boost)) {
            continue;
        }

        if (mirrored_reflection_instead_of_xillver()) {
            // TODO: remove this once everything is working; only for testing
            calc_zone_mirrored_flux(xill_param_zone, returnSpec, specZonesEmittedCombined[izone], izone);
        } else {
            // very dangerous; need to check that it is within the table boundaries as well
            // (will only be used for the model with the xillverRTR table (MOD_TYPE_RELXILLBB))
            xill_param_zone.radzone = static_cast<double>(izone);
            if (use_iongrad_relxillbb()) {
                set_xill_param_from_iongrad(xill_param_zone, ion_gradient, izone);
            }
            calc_zone_xillver_reflected_flux(xill_param_zone, relProfile, returnSpec, specZonesEmittedCombined[izone],
                                             izone, status);
        }


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
        for (int jj = 0; jj < n_ener_inp; jj++) {
            specOutput[jj] += specSingleZone[jj];
        }
    }

    // clean spectrum
    set_low_values_to_zero(specOutput, n_ener_inp);
    set_values_outside_to_zero(specOutput, ener_inp, n_ener_inp);

    if (shouldOutfilesBeWritten()) {
        write_bbret_diagnose_outfiles(xill_param, specZonesConvOutput, specZonesEmittedCombined,
                                      ener, n_ener, returnSpec, status);
    }

    delete_returnSpec2D(&returnSpec);
    delete_rel_params(rel_param);
    delete xill_param;
}

