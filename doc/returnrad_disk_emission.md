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

### Novikov-Thorne profile (`RELXILL_RRAD_BBODY_TPROFILE=3`)

The temperature profile of a relativistic thin disk around a Kerr black hole with zero torque at the inner edge
(Novikov & Thorne 1973; `disk_Tprofile_novikov_thorne()`, `TPROFILE_NOVIKOV_THORNE`), as used by the Xspec model
`kerrbb` (with `eta=0`). The emitted flux is \f$F = \frac{\dot{M}}{4\pi r} f\f$, with the closed form of \f$f\f$ from
Page & Thorne (1974, ApJ 191, 499, eq. 15n; see also Li et al. 2005, ApJS 157, 335, eqs. D11 and D13):

\f[
f \;=\; \frac{3}{2}\, \frac{1}{x^2 (x^3 - 3x + 2a)} \left[ x - x_0 - \frac{3}{2} a \ln\frac{x}{x_0}

- \sum_{i=1}^{3} \frac{3 (x_i - a)^2}{x_i \prod_{j \neq i} (x_i - x_j)} \ln\frac{x - x_i}{x_0 - x_i} \right]
  \f]

with \f$x = \sqrt{r}\f$, \f$x_0 = \sqrt{R_\mathrm{in}}\f$ (radii in units of \f$r_g = GM/c^2\f$), the spin
\f$a\f$, and the roots of \f$x^3 - 3x + 2a = 0\f$: \f$x_1 = 2\cos\left(\frac{1}{3}\arccos a - \frac{\pi}{3}\right)\f$,
\f$x_2 = 2\cos\left(\frac{1}{3}\arccos a + \frac{\pi}{3}\right)\f$, \f$x_3 = -2\cos\left(\frac{1}{3}\arccos a\right)\f$.
The temperature is \f$T(r) = kT_\mathrm{bb}\, \left(F(r) / F(R_\mathrm{max})\right)^{1/4}\f$:

- As for the alpha disk, \f$kT_\mathrm{bb}\f$ is the maximal temperature of the disk. The radius of the maximum,
  \f$R_\mathrm{max}\f$, is determined numerically (golden-section search between \f$R_\mathrm{in}\f$ and
  \f$4 R_\mathrm{in}\f$). For \f$R_\mathrm{in} = r_\mathrm{isco}\f$, it is at 1.59, 1.48, and 1.28
  \f$R_\mathrm{in}\f$ for \f$a\f$ = 0, 0.9, and 0.998 (1.36 \f$R_\mathrm{in}\f$ for the alpha disk).
- The inner edge must be at or outside the ISCO (\f$R_\mathrm{in} \geq r_\mathrm{isco}\f$).
- The closed form agrees with a numerical integration of \f$f\f$ (Li et al. 2005, eq. D13, with the circular orbits of
  Bardeen et al. 1972) to \f$2.5 \times 10^{-6}\f$ for \f$a\f$ = -0.9 to 0.998.
- Compared to the alpha disk with the same maximal temperature, the outer disk is hotter for low spin (by 22% for
  \f$a = 0\f$ at large radii) and cooler for high spin (by 8-13% for \f$a = 0.998\f$).

The figure compares the Novikov-Thorne profile with the alpha disk (created with
`test/manual/disk_temperature_profile/plot_temperature_profiles.py`, which uses the same equations as the code):

\image html disk_temperature_profiles.png "Temperature profiles normalized to the maximal temperature, for an inner edge
at the ISCO (top), and the ratio of the Novikov-Thorne profile to the alpha disk (bottom)."

With the Novikov-Thorne profile, the full model agrees with the relativistic Xspec model `kerrbb` without any free
parameter, see (2b).

## Tests and precision

The disk emission is tested in three steps, which correspond to the unit tests listed below:

- **(1) Temperature profile and normalization**: \ref spec_multicolor_disk on a fine radial grid compared to the Xspec
  model with the same temperature profile. This tests the implementation of the temperature profiles, including the
  definition of \f$kT_\mathrm{bb}\f$ and of the inner edge of the disk, the blackbody spectra, and how the
  normalization changes with \f$kT_\mathrm{bb}\f$.
- **(2a) Radial zones**: \ref spec_multicolor_disk on the 50 radial zones of the model compared to a fine radial grid
  over the same radii. This is a code-internal comparison (as (1) shows that the fine radial grid is correct) and
  gives the precision of the disk emission of the model due to its radial zones.
- **(2b) Full model**: the direct disk emission of the full `relxillBB_NS` model with the Novikov-Thorne profile
  compared to the relativistic disk model `kerrbb`. This tests how the disk emission is implemented in the model
  (temperature profile, radial zones, relativistic effects), without any free parameter.

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

**Current state and possible improvement.** The same 50 zones are used for the direct disk emission and for the
returning radiation (`specRet[izone]`) and its reflection, as they are the zones of the returning-radiation table.
Optimizing them (a different distribution or more zones) would therefore require changing the table, and has to take
the returning radiation into account as well. A better approach for later is to calculate the direct disk emission on
a fine radial grid, which is not limited by the table, and to use the coarse zones of the table only for the
returning radiation and its reflection.

### (2b) Full model (with relativistic effects)

The full model (`relxillBB_NS` with `boost=0`, i.e., only the direct emission of the disk, \f$R_\mathrm{in} =
r_\mathrm{isco}\f$, \f$R_\mathrm{out}\f$ = 1000 Rg) with the Novikov-Thorne profile and \f$kT_\mathrm{bb}\f$ = 1 keV
compared to the Xspec model `kerrbb` (`eta=0`, `hd=1`, `rflag=0`, `lflag=0`; \f$M\f$ = 10 \f$M_\odot\f$ and the
accretion rate chosen such that the maximal temperature of the disk is 1 keV), i.e., without any free parameter.
Shape deviation in 0.5-10 kTbb (50 energy bins):

| spin  | Novikov-Thorne, incl = 30 deg | incl = 60 deg | alpha disk, incl = 30 deg | incl = 60 deg |
|-------|-------------------------------|---------------|---------------------------|---------------|
| 0     | 0.35%                         | 0.52%         | 7.9%                      | 5.8%          |
| 0.9   | 0.81%                         | 0.56%         | 16%                       | 3.1%          |
| 0.998 | 0.36%                         | 0.10%         | 124%                      | 45%           |

- With the Novikov-Thorne profile, the implementation of the disk emission in the model (temperature profile,
  radial zones, relativistic effects) agrees with `kerrbb` to better than 1%.
- Lower energies are not compared: the outer radius of the model (1000 Rg, i.e., only \f$167\,r_\mathrm{isco}\f$ for
  \f$a = 0\f$) truncates the disk emission there, which gives a deviation of 15% in 0.2-10 kTbb for \f$a = 0\f$
  (2.3-2.5% for \f$a = 0.9\f$, 0.8% for \f$a = 0.998\f$).
- The alpha disk deviates by 3-124% for the same maximal temperature (and still by 5-10% for the best-fitting
  \f$kT_\mathrm{bb}\f$), due to the different temperature profile, see the figure above.

Before, the full model was compared to the non-relativistic models `diskbb` and `ezdiskbb`. These comparisons were
limited by the relativistic effects (5-8% for \f$R_\mathrm{in} = 10\,r_\mathrm{isco}\f$ and 40-50 deg, where the
gravitational redshift and the Doppler shift partly cancel) and were replaced by the comparison with `kerrbb`.

## Unit tests

In `test/unit/tests-bbody-returnrad.cpp`; all disk emission tests run with `bin/tests "[diskemission]"`:

- **(1)** *" disk emission (1): temperature profile on a fine radial grid vs. Xspec diskbb and ezdiskbb "*: shape
  deviation below 0.3% (alpha disk) and 0.05% (diskbb), and normalization independent of kTbb to 0.1%.
- **(2a)** *" disk emission (2a): radial zones of the model vs. a fine radial grid "*: shape deviation below 4%
  (alpha disk) and 2% (diskbb).
- **(2b)** *" disk emission (2b): full relxillBB_NS model with the Novikov-Thorne profile vs. Xspec kerrbb "*
  (`[kerrbb]`): for \f$a\f$ = 0, 0.9, 0.998 and an inclination of 30 and 60 deg, the shape deviation of the
  Novikov-Thorne profile is below 2% in 0.5-10 kTbb, and the one of the alpha disk above 2%.

The reference spectra of the Xspec models are in `test/unit/refdata/`: `refvalues_diskbb_Tin{0.5,1.0,2.0}.fits`
(column `diskbb`) and `refvalues_ezdiskbb_Tmax{0.5,1.0,2.0}.fits` (column `ezdiskbb`), created with
`create_refvalues_diskbb.sl`, and `refvalues_kerrbb_a{0.000,0.900,0.998}_i{30,60}.fits` (column `kerrbb`, with the
parameters in the header), created with `create_refvalues_kerrbb.sl`.

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

- Precision of the 50 radial zones, in particular the coarse zones at large radii. The zones are those of the
  returning-radiation table, so changing them also changes the returning radiation. The preferred solution is to
  calculate the direct disk emission on a fine radial grid and to use the zones of the table only for the returning
  radiation and its reflection, see (2a).
- \f$R_\mathrm{max} = 1.5^{4/5} R_\mathrm{in}\f$ instead of \f$\frac{49}{36} R_\mathrm{in}\f$: the maximal temperature
  is 0.03% higher than \f$kT_\mathrm{bb}\f$ (0.2% in the spectrum at 10 kTbb, see (1)).
- Which temperature profile should be the default: the alpha disk is Newtonian, whereas the relativistic effects are
  calculated in the Kerr metric (see the Novikov-Thorne profile above).
- The parameter `shiftTmaxRRad` of `relxillBB_NS` is read, but not used in the code.
