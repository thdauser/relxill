# Thermal Returning Radiation: Disk Emission and Temperature Profile {#returnrad_disk_emission}

This page documents how the direct thermal emission of the accretion disk is calculated for `relxillBB` and
`relxillBB_NS`, how the available temperature profiles are defined, how they are tested, and how precise the current
implementation is (state of 2026-10). The normalization of the reflected spectrum is described in
\ref returnrad_normalizations.

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

\ref spec_multicolor_disk calculates the same disk emission without relativistic effects, for any given radial zones:
it uses the same temperature profiles and blackbody spectra, and weights each zone with its flat-space ring area
\f$\pi (r_\mathrm{hi}^2 - r_\mathrm{lo}^2)\f$. It is only used in the unit tests, for the comparison with the Xspec
models and for the precision of the radial zones.

## Temperature profiles

For both profiles, the temperature of each zone is evaluated at its mean radius \f$\frac{1}{2}(r_\mathrm{lo} +
r_\mathrm{hi})\f$, and each zone emits a single blackbody of this temperature.

### Alpha-disk profile (zero torque at the inner edge; default, `RELXILL_RRAD_BBODY_TPROFILE=1`)

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
  `ezdiskbb`, which uses the same temperature profile. At large radii,
  \f$T(r) \approx 2.05\,kT_\mathrm{bb}\,(r/R_\mathrm{in})^{-3/4}\f$.

In the code, this profile is called `TPROFILE_ALPHA`. The Xspec model `ezdiskbb` is only used as the reference in the
unit tests (and in the names of their reference files).

### diskbb profile (for testing, `RELXILL_RRAD_BBODY_TPROFILE=2`)

The temperature profile of the Xspec model `diskbb` (`disk_Tprofile_diskbb()`, `TPROFILE_DISKBB`):

\f[
T(r) \;=\; kT_\mathrm{bb} \left(\frac{r}{R_\mathrm{in}}\right)^{-3/4}
\f]

with \f$kT_\mathrm{bb}\f$ the temperature at the inner edge \f$R_\mathrm{in}\f$ of the disk (the inner edge of the
first zone), as the parameter `Tin` of `diskbb`.

## Tests and precision

The disk emission is tested in three steps, which correspond to the unit tests listed below:

- **(1) Temperature profile and normalization**: \ref spec_multicolor_disk on a fine radial grid compared to the Xspec
  model with the same temperature profile. This tests the implementation of the temperature profiles, including the
  definition of \f$kT_\mathrm{bb}\f$ and of the inner edge of the disk, the blackbody spectra, and how the
  normalization changes with \f$kT_\mathrm{bb}\f$.
- **(2a) Radial zones**: \ref spec_multicolor_disk on the 50 radial zones of the model compared to a fine radial grid
  over the same radii. This is a code-internal comparison (as (1) shows that the fine radial grid is correct) and
  gives the precision of the disk emission of the model due to its radial zones.
- **(2b) Full model**: the direct disk emission of the full `relxillBB_NS` model (with relativistic effects) compared
  to the Xspec models. This tests how the disk emission is implemented in the model, but is only approximate because
  of the relativistic effects.

The precision is given as the shape deviation: the maximal relative deviation of the ratio of two spectra in the
given band, after normalizing the ratio to its mean in the band. All values are for an energy grid of 50 logarithmic
bins from 0.01 to 20 keV (as in the unit tests), and \f$kT_\mathrm{bb}\f$ = 0.5, 1, and 2 keV. The spectra of
\ref spec_multicolor_disk are calculated on the fine energy grid of the model and rebinned (evaluating the blackbody at
the bin centres of the coarse grid directly would add a deviation of about 3%).

### (1) Temperature profile and normalization

Fine logarithmic radial grid with 1000 zones from \f$R_\mathrm{in}\f$ to \f$10^4\,R_\mathrm{in}\f$ (large enough
for the outer edge not to affect the band 0.1-10 kTbb):

| temperature profile | reference  | shape deviation (0.1-10 kTbb) | normalization for kTbb = 0.5, 2 keV (rel. to 1 keV) |
|---------------------|------------|-------------------------------|-----------------------------------------------------|
| alpha disk          | `ezdiskbb` | 0.16-0.18%                    | 0.99996, 0.99996                                    |
| diskbb              | `diskbb`   | 0.007%                        | 1.00000, 1.00000                                    |

- Both temperature profiles agree with the Xspec models, including the definition of \f$kT_\mathrm{bb}\f$ and of the
  inner edge, and the normalization changes with \f$kT_\mathrm{bb}\f$ in the same way.
- For the alpha disk, the remaining deviation (0.2% at 10 kTbb) is due to the radius of the maximal temperature,
  \f$1.5^{4/5} R_\mathrm{in}\f$ instead of \f$\frac{49}{36} R_\mathrm{in}\f$. With the exact radius, the alpha disk and
  `ezdiskbb` agree to 0.01%.
- The outer radius of the model (1000 Rg) only affects the lowest energies: for the alpha disk with
  \f$R_\mathrm{in} = r_\mathrm{isco}(a=0.998)\f$, the deviation is 1% at 0.1 kTbb and below 0.03% above 0.2 kTbb.

### (2a) Radial zones of the model

The 50 radial zones of the model (\f$R_\mathrm{in} = r_\mathrm{isco}(a=0.998)\f$ to 1000 Rg) compared to a fine
logarithmic grid with 500 zones over the same radii (kTbb = 1 keV, the shape only depends on \f$E/kT_\mathrm{bb}\f$):

| temperature profile | shape deviation (0.1-10 kTbb) |
|---------------------|-------------------------------|
| alpha disk          | 3.5%                          |
| diskbb              | 1.5-1.7%                      |

- The zones of the returning-radiation table are fine close to the inner edge (\f$r_\mathrm{hi}/r_\mathrm{lo}\f$ =
  1.04), but coarse at large radii: the last zones are 227-417 Rg and 417-1000 Rg (\f$r_\mathrm{hi}/r_\mathrm{lo}\f$
  up to 2.4).
- The coarse outer zones mainly affect the lowest energies. For the alpha disk, the deviation is 3.5% at 0.12 kTbb,
  1.8% at 0.19 kTbb, and within 1% above 0.3 kTbb.
- The alpha disk is more affected than the diskbb profile, as its outer disk is relatively hotter (factor 2.05 at
  large radii, see above), so the outer zones emit at higher energies relative to kTbb.

### (2b) Full model (with relativistic effects)

The full model (`relxillBB_NS` with `boost=0`, i.e., only the direct emission of the disk) additionally contains the
relativistic effects (energy shift, light bending), which the Xspec models do not. The comparison is therefore only
approximate:

- **diskbb profile vs. `diskbb`** (0.15-10 kTbb): 4.9-7.9% for \f$R_\mathrm{in} = 10\,r_\mathrm{isco}\f$ and an
  inclination of 50 deg, where the gravitational redshift and the Doppler shift partly cancel. For
  \f$R_\mathrm{in}\f$ = 8-12 \f$r_\mathrm{isco}\f$ and an inclination of 40-50 deg, the deviation is 7-20%, and
  27-41% for 30 deg. The alpha disk deviates by 33-40%, i.e., the comparison is sensitive to the temperature profile.
- **alpha disk vs. `ezdiskbb`** (0.3-10 kTbb): 7.4-7.5% for \f$R_\mathrm{in} = 10\,r_\mathrm{isco}\f$ and an
  inclination of 40 deg; the diskbb profile deviates by 30-31%. Lower energies are not compared: for
  \f$R_\mathrm{in} = 10\,r_\mathrm{isco}\f$, the disk only extends over a factor of 80 in radius
  (\f$R_\mathrm{out}\f$ = 1000 Rg), which truncates the emission of the relatively hot outer disk (even without
  relativistic effects, the deviation would be 56% in 0.1-10 kTbb).

A comparison with a relativistic disk model (e.g., `kerrbb`) instead of the non-relativistic Xspec models is still to
be investigated.

## Unit tests

In `test/unit/tests-bbody-returnrad.cpp`; all disk emission tests run with `bin/tests "[diskemission]"`:

- **(1)** *" disk emission (1): temperature profile on a fine radial grid vs. Xspec diskbb and ezdiskbb "*: shape
  deviation below 0.3% (alpha disk) and 0.05% (diskbb), and normalization independent of kTbb to 0.1%.
- **(2a)** *" disk emission (2a): radial zones of the model vs. a fine radial grid "*: shape deviation below 4%
  (alpha disk) and 2% (diskbb).
- **(2b)** *" relxillBB_NS direct disk emission vs. Xspec diskbb for the diskbb and zero-torque temperature profiles "*
  (`[diskbb]`) and *" relxillBB_NS direct disk emission vs. Xspec ezdiskbb for the zero-torque and diskbb temperature
  profiles "* (`[ezdiskbb]`): the matching profile deviates by less than 10%, the other profile by more than 10%.

The reference spectra of the Xspec models (norm = 1) are in `test/unit/refdata/`:
`refvalues_diskbb_Tin{0.5,1.0,2.0}.fits`
(column `diskbb`) and `refvalues_ezdiskbb_Tmax{0.5,1.0,2.0}.fits` (column `ezdiskbb`), created with
`create_refvalues_diskbb.sl`.

## History

- 2026-10: `disk_temperature_diskbb()` calculated \f$(kT_\mathrm{bb}\, r/R_\mathrm{in})^{-3/4}\f$ instead of
  \f$kT_\mathrm{bb} (r/R_\mathrm{in})^{-3/4}\f$, i.e., an inner temperature of \f$kT_\mathrm{bb}^{-3/4}\f$. This had
  no effect for kTbb = 1 keV, the only value that was tested.
- 2026-10: for the diskbb profile, kTbb was the temperature at the mean radius of the first zone instead of the inner
  edge of the disk, i.e., the disk was 1.5% hotter than `diskbb` (shape deviation of 8-9%). On the coarse energy grid
  of the unit test, this was partly compensated by the evaluation at the bin centres (3.7% in total).
- Both only affected `RELXILL_RRAD_BBODY_TPROFILE=2` and the test function `spec_diskbb()` (now
  \ref spec_multicolor_disk), not the default alpha-disk profile.

## Open issues

- Precision of the 50 radial zones, in particular the coarse zones at large radii: a better distribution of the
  zones, or as few additional zones as possible. The zones are those of the returning-radiation table, so changing
  them also affects how the returning fractions are tabulated or interpolated.
- \f$R_\mathrm{max} = 1.5^{4/5} R_\mathrm{in}\f$ instead of \f$\frac{49}{36} R_\mathrm{in}\f$: the maximal temperature
  is 0.03% higher than \f$kT_\mathrm{bb}\f$ (0.2% in the spectrum at 10 kTbb, see (1)).
- A test of the full model against a relativistic disk model (e.g., `kerrbb`), see (2b).
- The parameter `shiftTmaxRRad` of `relxillBB_NS` is read, but not used in the code.
