# Thermal Returning Radiation: Disk Emission and Temperature Profile {#returnrad_disk_emission}

This page documents how the direct thermal emission of the accretion disk is calculated for `relxillBB` and
`relxillBB_NS`, how the available temperature profiles are defined, and how precise the current implementation is
(state of 2026-10). The normalization of the reflected spectrum is described in \ref returnrad_normalizations.

## Where it is calculated

\ref spec_returnrad_blackbody (`src/Relreturn_BlackBody.cpp`) calculates the spectrum emitted by each radial zone of
the disk:

1. `get_rrad_fractions(spin, Rin, Rout)` sets up the radial zones of the returning-radiation table: 50 zones,
   starting at the inner disk radius `Rin` and ending at `Rout` (at most 1000 Rg).
2. `get_temperature_profile_disk_zones()` calls `get_tprofile()` (`src/Relphysics.cpp`) for the temperature of each
   zone. The profile is selected with the environment variable `RELXILL_RRAD_BBODY_TPROFILE` (see below).
3. `get_bbody_specs_binned()` calculates the blackbody photon spectrum of each zone, integrated over the energy bins,
   on the energy grid of the convolution (`N_ENER_CONV` = 4096 bins from 0.00035 to 2000 keV, i.e., bins of 0.4%).

The result is `specPri[izone]` (per unit area of the zone). In \ref relxillbb_kernel it is folded with the relativistic
line profile of each zone and rebinned to the energy grid of the model. It is also the input for the returning
spectrum `specRet[izone]`.

\ref spec_diskbb performs the same calculation with the diskbb temperature profile and weights the zones with their
flat-space ring area \f$\pi (r_\mathrm{hi}^2 - r_\mathrm{lo}^2)\f$. It is only used in the unit tests, for the
comparison with the diskbb model of Xspec.

## Temperature profiles

For both profiles, the temperature of each zone is evaluated at its mean radius \f$\frac{1}{2}(r_\mathrm{lo} +
r_\mathrm{hi})\f$, and each zone emits a single blackbody of this temperature.

### Zero-torque profile (default, `RELXILL_RRAD_BBODY_TPROFILE=1`)

The temperature profile of a standard alpha disk with zero torque at the inner edge (`disk_Tprofile_alpha()`):

\f[
T(r) \;=\; K \left(\frac{r}{R_\mathrm{in}}\right)^{-3/4} \left(1 - \sqrt{\frac{R_\mathrm{in}}{r}}\right)^{1/4}
\f]

- \f$R_\mathrm{in}\f$ is the inner radius of the disk (model parameter `Rin`), at which the radial zones start.
  The temperature is zero at \f$R_\mathrm{in}\f$.
- \f$K\f$ is chosen such that \f$T(R_\mathrm{max}) = kT_\mathrm{bb}\f$ at
  \f$R_\mathrm{max} = 1.5^{4/5} R_\mathrm{in} = 1.383\,R_\mathrm{in}\f$ (Poutanen et al. 2007). The exact maximum of
  the profile is at \f$\frac{49}{36} R_\mathrm{in} = 1.361\,R_\mathrm{in}\f$, so the maximal temperature of the disk
  is \f$1.0003\,kT_\mathrm{bb}\f$.
- \f$kT_\mathrm{bb}\f$ is therefore the maximal temperature of the disk, as the parameter `T_max` of the Xspec model
  `ezdiskbb`. At large radii, \f$T(r) \approx 2.05\,kT_\mathrm{bb}\,(r/R_\mathrm{in})^{-3/4}\f$.

### diskbb profile (for testing, `RELXILL_RRAD_BBODY_TPROFILE=2`)

The temperature profile of the Xspec model `diskbb` (`disk_Tprofile_diskbb()`):

\f[
T(r) \;=\; kT_\mathrm{bb} \left(\frac{r}{R_\mathrm{in}}\right)^{-3/4}
\f]

with \f$kT_\mathrm{bb}\f$ the temperature at the inner edge \f$R_\mathrm{in}\f$ of the disk (the inner edge of the
first zone), as the parameter `Tin` of `diskbb`.

## Precision of the current implementation

The precision is given as the shape deviation from the corresponding Xspec model: the maximal relative deviation of
the ratio model/reference in the given band, after normalizing the ratio to its mean in the band. All values are for
an energy grid of 50 logarithmic bins from 0.01 to 20 keV (as in the unit tests), \f$kT_\mathrm{bb}\f$ = 0.5, 1, and
2 keV, and \f$R_\mathrm{in} = r_\mathrm{isco}(a=0.998)\f$, with the spectra calculated on the fine energy grid of the
model.

### Disk emission without relativistic effects

| temperature profile | calculation                                      | reference  | 0.1-10 kTbb  | 0.3-10 kTbb  |
|---------------------|--------------------------------------------------|------------|--------------|--------------|
| zero-torque         | continuous disk, \f$R_\mathrm{out} \to \infty\f$ | `ezdiskbb` | 0.2%         | 0.2%         |
| zero-torque         | continuous disk, \f$R_\mathrm{out}\f$ = 1000 Rg  | `ezdiskbb` | 0.3%         | 0.2%         |
| zero-torque         | **50 radial zones (as in the model)**            | `ezdiskbb` | **3.8-4.1%** | **1.0%**     |
| diskbb              | continuous disk, \f$R_\mathrm{out}\f$ = 1000 Rg  | `diskbb`   | 0.01%        | 0.01%        |
| diskbb              | **50 radial zones (as in the model)**            | `diskbb`   | **1.5-1.7%** | **0.3-0.4%** |

Both temperature profiles therefore agree with the Xspec models, including the definition of
\f$kT_\mathrm{bb}\f$ and of the inner edge of the disk. The deviation in the model comes from the discretization into
50 radial zones:

- The zones of the returning-radiation table are fine close to the inner edge (\f$r_\mathrm{hi}/r_\mathrm{lo}\f$ =
  1.04), but coarse at large radii: the last zones are 227-417 Rg and 417-1000 Rg (\f$r_\mathrm{hi}/r_\mathrm{lo}\f$
  up to 2.4).
- The coarse outer zones mainly affect the lowest energies. For the zero-torque profile, the ratio to `ezdiskbb` is
  0.96 at 0.12 kTbb and 0.98 at 0.19 kTbb, and within 1% above 0.5 kTbb.
- The zero-torque profile is more affected than the diskbb profile, as its outer disk is relatively hotter (factor
  2.05 at large radii, see above), so the outer zones emit at higher energies relative to kTbb.

The blackbody of each zone is evaluated at the bin centres of the energy grid. On the fine energy grid of the model,
this is negligible. On a coarse grid (50 bins from 0.01 to 20 keV, i.e., bins of 16.5%), it would add a deviation of
about 3%. Therefore, the unit test calculates \ref spec_diskbb on the energy grid of the model and rebins it.

### Full model (with relativistic effects)

The full model (`relxillBB_NS` with `boost=0`, i.e., only the direct emission of the disk) additionally contains the
relativistic effects (energy shift, light bending), which the Xspec models do not. The comparison is therefore only
approximate (band 0.15-10 kTbb):

- **diskbb profile vs. `diskbb`**: 4.9-7.9% for \f$R_\mathrm{in} = 10\,r_\mathrm{isco}\f$ and an inclination of 50
  deg, where the gravitational redshift and the Doppler shift partly cancel. For \f$R_\mathrm{in}\f$ = 8-12
  \f$r_\mathrm{isco}\f$ and an inclination of 40-50 deg, the deviation is 7-20%, and 27-41% for 30 deg.
- **zero-torque profile vs. `diskbb`**: 33-40%, i.e., the comparison is sensitive to the temperature profile.
- **zero-torque profile vs. `ezdiskbb`**: not a precise comparison in the full model. For
  \f$R_\mathrm{in} = 10\,r_\mathrm{isco}\f$, the disk only extends over a factor of 80 in radius
  (\f$R_\mathrm{out}\f$ = 1000 Rg), which truncates the emission of the relatively hot outer disk. Even without
  relativistic effects, this gives a deviation of 56% (0.1-10 kTbb) or 8% (0.3-10 kTbb). For smaller
  \f$R_\mathrm{in}\f$, the relativistic effects dominate.

## Unit tests

In `test/unit/tests-bbody-returnrad.cpp` (run with `bin/tests "[bbret]"`):

- *" DiskBB Spectrum "*: \ref spec_diskbb for kTbb = 0.5, 1, 2 keV, compared to `diskbb`. The shape deviation must
  be below 2% in 0.1-10 kTbb.
- *" relxillBB_NS direct disk emission vs. Xspec diskbb for the diskbb and zero-torque temperature profiles "*
  (`[diskbb]`): full model for \f$R_\mathrm{in} = 10\,r_\mathrm{isco}\f$ and 50 deg. With the diskbb profile, the
  shape deviation from `diskbb` must be below 10% in 0.15-10 kTbb, and with the zero-torque profile above 10%.

The reference spectra of `diskbb` (norm = 1) are in `test/unit/refdata/`: `refvalues_diskbb.fits` (kTbb = 1 keV,
column `diskbb`; its column `spec` is an earlier output of \ref spec_diskbb and is no longer used),
`refvalues_diskbb_Tin0.5.fits`, and `refvalues_diskbb_Tin2.0.fits`, created with `create_refvalues_diskbb.sl`.

The comparison of the zero-torque profile with `ezdiskbb` (table above) is not yet a unit test.

## History

- 2026-10: `disk_temperature_diskbb()` calculated \f$(kT_\mathrm{bb}\, r/R_\mathrm{in})^{-3/4}\f$ instead of
  \f$kT_\mathrm{bb} (r/R_\mathrm{in})^{-3/4}\f$, i.e., an inner temperature of \f$kT_\mathrm{bb}^{-3/4}\f$. This had
  no effect for kTbb = 1 keV, the only value that was tested.
- 2026-10: for the diskbb profile, kTbb was the temperature at the mean radius of the first zone instead of the inner
  edge of the disk, i.e., the disk was 1.5% hotter than `diskbb` (shape deviation of 8-9%). On the coarse energy grid
  of the unit test, this was partly compensated by the evaluation at the bin centres (3.7% in total).
- Both only affected `RELXILL_RRAD_BBODY_TPROFILE=2` and \ref spec_diskbb, not the default zero-torque profile.

## Open issues

- Precision of the 50 radial zones, in particular the coarse zones at large radii: a better distribution of the
  zones, or as few additional zones as possible.
- \f$R_\mathrm{max} = 1.5^{4/5} R_\mathrm{in}\f$ instead of \f$\frac{49}{36} R_\mathrm{in}\f$ (0.03% in the maximal
  temperature).
- The parameter `shiftTmaxRRad` of `relxillBB_NS` is read, but not used in the code.
