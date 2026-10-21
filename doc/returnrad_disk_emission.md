# Thermal Returning Radiation: Disk Emission and Temperature Profile {#returnrad_disk_emission}

This page documents how the direct thermal emission of the accretion disk is calculated for `relxillBB` and
`relxillBB_NS`, how the available temperature profiles are defined, how they are tested, and how precise the current
implementation is (state of 2026-10). The normalization of the reflected spectrum is described in
\ref returnrad_normalizations.

## Where it is calculated

The direct emission of the disk and the reflection of the returning radiation are calculated separately, but from the
same emitted primary (see (3)). Both are integrated over the `N_FRAD` = 1000 radial points of the relativistic line
profiles (from \f$R_\mathrm{in}\f$ to \f$R_\mathrm{out}\f$, uniform in \f$1/\sqrt{r}\f$, `get_fine_radial_grid()`),
and each point emits a blackbody with the temperature of the temperature profile at its radius,
\f$T(r) = kT_\mathrm{bb}\, t(r)\f$ (`get_tprofile()` in `src/Relphysics.cpp`, selected with the environment variable
`RELXILL_RRAD_BBODY_TPROFILE`, see below).

**Direct emission.** In \ref relxillbb_kernel, `calc_disk_emission_observed()` calculates the observed direct emission
of the whole disk with a single convolution:

- Per unit area, the bin-integrated photon spectrum of a blackbody of temperature \f$T\f$ is \f$T^3\f$ times the one of
  1 keV, shifted in energy by the factor \f$T\f$. The direct emission of the disk is therefore the blackbody of
  \f$kT_\mathrm{bb}\f$ convolved with a kernel, which is the relativistic line profile of each radial point at the line
  energy \f$t(r)\f$, weighted with \f$t(r)^3\f$ (`get_disk_emission_kernel()`, using
  \ref calc_relline_profile_radial_line_energy).
- The kernel does not depend on \f$kT_\mathrm{bb}\f$. It is cached and only re-calculated if a relativistic parameter
  or the temperature profile changes.
- The convolution is done with an FFT on the energy grid of the convolution (`N_ENER_CONV` = 4096 bins from 0.00035
  to 2000 keV, i.e., bins of 0.4%). As the FFT convolution is cyclic and the kernel shifts the energies down by up to
  a factor of about 100 (outer disk), the blackbody is set to zero in the bins that would be shifted outside the grid.
  This only changes the spectrum below about 0.05 keV (0.1 keV for the diskbb profile).
- The number of photons is normalized on the whole energy grid. The standard normalization of the FFT convolution only
  in 0.01-1000 keV would overestimate the flux by up to 1%, as the cool outer disk emits a part of its photons below
  0.01 keV.

**Primary of the zones and returning radiation.** \ref spec_returnrad_blackbody calculates the spectrum emitted by
each radial zone and the returning radiation:

1. `get_rrad_fractions(spin, Rin, Rout)` sets up the radial zones of the returning-radiation table: 50 zones,
   starting at the inner disk radius `Rin` and ending at `Rout` (at most 1000 Rg).
2. `get_disk_emission_zones()` bins the emission of the same radial points (with the same temperatures) into the
   zones: the cell of each point (between the midpoints to its neighbors) is split at the zone boundaries in proportion
   to its area. This is the primary of each zone, `specPri[izone]` (per unit area), calculated on the energy grid of the
   returning radiation (500 bins from 0.001 to 500 keV, set to zero outside 0.01-1000 keV).
3. \ref get_returnrad_specs calculates the returning radiation `specRet[izone]` from exactly this primary, which sets
   the normalization of the reflected spectrum of each zone (see \ref returnrad_normalizations).

In \ref relxillbb_kernel, the reflected spectrum of each zone is convolved with the relativistic line profile of the
zone. `calc_relline_profile()` splits the cells of the radial points at the zone boundaries in the same way, such that
each zone covers exactly its radial range \f$[r_\mathrm{lo}, r_\mathrm{hi}]\f$, as for the returning-radiation table
(see (3)). For `boost=0`, i.e., only the direct emission, the zones are not calculated.

**Diagnostic output.** With `RELXILL_WRITE_FILES=1`, the spectra of the zones are written to
`__relxillOutput-rrad-bbody-*.fits` (used by the scripts in `test/manual/returnrad/`):

- `rframe-bbodyEmit.fits`: the primary of each zone `specPri[izone]` (rest frame).
- `obs-bbodyEmit.fits`: the primary of each zone, convolved with the line profile of the zone (observed). The sum
  over all zones agrees with the direct emission of the model to about 0.1% in the number of photons and the energy
  (see (3c)); in the Wien tail, it deviates by up to a few percent per energy bin (2.6% at 10 kTbb for
  \f$a = 0.998\f$ and an inclination of 40 deg).
- `rframe-bbodyRet.fits`: the returning radiation of each zone `specRet[izone]` (rest frame).
- `rframe-reflect.fits` and `obs-reflect.fits` (only for `boost` \f$\neq\f$ 0): the reflected spectrum of each zone,
  scaled with \f$|b|\f$, in the rest frame and observed (`rframe-mirror-reflect.fits` and `obs-mirror-reflect.fits`
  for `RELXILL_RRAD_BBODY_MIRRORED=1`).

The total spectrum of a zone is the sum of its primary and its reflected spectrum.

\ref spec_multicolor_disk calculates the disk emission (with \ref get_disk_emission_zones_binned) without relativistic
effects, for any given radial zones, and weights each zone with its flat-space ring area
\f$\pi (r_\mathrm{hi}^2 - r_\mathrm{lo}^2)\f$. It is only used in the unit tests, for the comparison with the Xspec
models.

## Temperature profiles

For all profiles, the temperature is evaluated at each radial point. A single temperature per zone (at its mean radius,
`get_temperature_profile_disk_zones()`) is only written to the diagnostic output.

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

with \f$kT_\mathrm{bb}\f$ the temperature at the inner edge \f$R_\mathrm{in}\f$ of the disk, as the parameter `Tin`
of `diskbb`.

### Novikov-Thorne profile (`RELXILL_RRAD_BBODY_TPROFILE=3`)

The temperature profile of a relativistic thin disk around a Kerr black hole with zero torque at the inner edge
(Novikov & Thorne 1973; `disk_Tprofile_novikov_thorne()`, `TPROFILE_NOVIKOV_THORNE`), as used by the Xspec model
`kerrbb` (with `eta=0`). The emitted flux is \f$F = \frac{\dot{M}}{4\pi r} f\f$, with the closed form of \f$f\f$ from
Page & Thorne (1974, ApJ 191, 499, eq. 15n; see also Li et al. 2005, ApJS 157, 335, eqs. D11 and D13):

\f[
f \;=\; \frac{3}{2}\, \frac{1}{x^2 (x^3 - 3x + 2a)} \left[ x - x_0 - \frac{3}{2} a \ln\frac{x}{x_0} -
\sum_{i=1}^{3} \frac{3 (x_i - a)^2}{x_i \prod_{j \neq i} (x_i - x_j)} \ln\frac{x - x_i}{x_0 - x_i} \right]
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
`test/manual/disk_temperature_profile/plot_temperature_profiles.sl`, which uses the same equations as the code):

\image html disk_temperature_profiles.png "Temperature profiles normalized to the maximal temperature, for an inner edge
at the ISCO (top), and the ratio of the Novikov-Thorne profile to the alpha disk (bottom)."

With the Novikov-Thorne profile, the full model agrees with the relativistic Xspec model `kerrbb` without any free
parameter, see (2b).

## Runtime

Measured with an optimized build (`-O2`, 2026-10, \f$a = 0.998\f$). The absolute values scatter by up to 50% between
runs, depending on the load of the machine:

| calculation                                                                       | alpha disk | Novikov-Thorne |
|-----------------------------------------------------------------------------------|------------|----------------|
| temperature at the 1000 radial points (incl. the search for \f$R_\mathrm{max}\f$) | 0.04 ms    | 0.17 ms        |
| primary of the 50 zones from the 1000 radial points (500 energy bins)             | 3.3 ms     | 3.4 ms         |

- The temperature profile is calculated for the primary of the zones in each model evaluation, and for the kernel if a
  relativistic parameter changes. The Novikov-Thorne profile therefore adds about 0.1-0.3 ms, i.e., less than 0.3% of
  a model evaluation of 100 ms. Pre-calculating \f$R_\mathrm{max}\f$ is not necessary.
- The kernel of the direct emission (1000 radial points, 4096 energy bins) takes 4, 9, 20, and 32 ms for an inclination
  of 10, 30, 60, and 80 deg. It integrates the same radial points as the line profiles of the zones and takes about as
  long. As it is cached, it is only calculated if a relativistic parameter changes.
- Compared to the calculation of the direct emission per zone (before 2026-10, see History), a full `relxillBB_NS`
  evaluation takes 6% (30 deg) to 11% (60 deg) longer if a relativistic parameter changes, and the same time if only
  \f$kT_\mathrm{bb}\f$ changes. For `boost=0`, it is about 12% faster, as the zones are not calculated.

## Tests and precision

The disk emission is tested in four steps, which correspond to the unit tests listed below:

- **(1) Temperature profile and normalization**: \ref spec_multicolor_disk on a fine radial grid compared to the Xspec
  model with the same temperature profile. This tests the implementation of the temperature profiles, including the
  definition of \f$kT_\mathrm{bb}\f$ and of the inner edge of the disk, the blackbody spectra, and how the
  normalization changes with \f$kT_\mathrm{bb}\f$.
- **(2a) Kernel of the direct emission**: the kernel compared to the line profiles of the zones. This is a
  code-internal test of the integration of the kernel.
- **(2b) Full model**: the direct disk emission of the full `relxillBB_NS` model with the Novikov-Thorne profile
  compared to the relativistic disk model `kerrbb`. This tests how the disk emission is implemented in the model
  (temperature profile, relativistic effects, and the kernel), without any free parameter.
- **(3) Direct and reflected radiation from the same emitted primary**: the primary of the zones, from which the
  returning radiation and therefore the reflection is calculated, compared to the emission of the radial points and to
  the direct emission of the model. This is a code-internal test of the consistency of the normalization of the direct
  and the reflected radiation.

The precision of the spectral shape is given as the shape deviation: the maximal relative deviation of the ratio of
two spectra in the given band, after normalizing the ratio to its mean in the band. Unless stated otherwise, the values
are for an energy grid of 50 logarithmic bins from 0.01 to 20 keV (as in the unit tests). The spectra of
\ref spec_multicolor_disk are calculated on the fine energy grid of the model and rebinned (evaluating the blackbody at
the bin centres of the coarse grid directly would add a deviation of about 3%).

### (1) Temperature profile and normalization

Fine logarithmic radial grid with 1000 zones from \f$R_\mathrm{in}\f$ to \f$10^4\,R_\mathrm{in}\f$ (large enough
for the outer edge not to affect the band 0.1-10 kTbb), for \f$kT_\mathrm{bb}\f$ = 0.5, 1, and 2 keV:

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

### (2a) Kernel of the direct emission

- With a line energy of 1 keV and the same weight for all radial points, the kernel is identical to the sum of the line
  profiles of the 50 zones of the model (to \f$2 \times 10^{-15}\f$), i.e., it integrates the same radial points with
  the same areas and relativistic effects.
- A constant line energy of 0.5 keV shifts the kernel by this factor (mean of \f$\ln E\f$ shifted by -0.69316 instead
  of -0.69315), and keeps the number of photons (to \f$5 \times 10^{-5}\f$).
- Compared to a converged reference (the calculation per zone with 400 zones uniform in \f$1/\sqrt{r}\f$, with the
  area of each zone assigned exactly), the direct emission of the model agrees to 0.01-0.04% in 0.2-10 kTbb, for
  \f$a\f$ = 0, 0.904277, 0.997676 and an inclination of 30 and 60 deg (Novikov-Thorne profile). The radial points of
  the kernel are fine enough also close to the ISCO, where the temperature and the relativistic effects change
  strongly with the radius.

### (2b) Full model (with relativistic effects)

The full model (`relxillBB_NS` with `boost=0`, i.e., only the direct emission of the disk, \f$R_\mathrm{in} =
r_\mathrm{isco}\f$, \f$R_\mathrm{out}\f$ = 1000 Rg) with the Novikov-Thorne profile and \f$kT_\mathrm{bb}\f$ = 1 keV
compared to the Xspec model `kerrbb` (`eta=0`, `hd=1`, `rflag=0`, `lflag=0`; \f$M\f$ = 10 \f$M_\odot\f$ and the
accretion rate chosen such that the maximal temperature of the disk is 1 keV), i.e., without any free parameter. The
spins are on the grid of the `kerrbb` table, to avoid its interpolation in spin (see below). Shape deviation in 0.5-10
kTbb:

| spin     | Novikov-Thorne, incl = 30 deg | incl = 60 deg | alpha disk, incl = 30 deg | incl = 60 deg |
|----------|-------------------------------|---------------|---------------------------|---------------|
| 0        | 0.45%                         | 0.59%         | 7.7%                      | 5.9%          |
| 0.904277 | 0.19%                         | 0.39%         | 15%                       | 2.5%          |
| 0.997676 | 1.44%                         | 0.41%         | 111%                      | 40%           |

- With the Novikov-Thorne profile, the direct emission of the model agrees with `kerrbb` to 0.2-1.4%. This is the
  precision of `kerrbb` itself: the converged reference of (2a) deviates by the same amount (0.4-1.4%), see below.
- Lower energies are not compared: the outer radius of the model (1000 Rg, i.e., only \f$167\,r_\mathrm{isco}\f$ for
  \f$a = 0\f$, whereas `kerrbb` assumes \f$r_\mathrm{out} = 10^6\,M\f$) truncates the disk emission there, which
  gives a deviation of 16% at 0.2 kTbb for \f$a = 0\f$. For \f$a \geq 0.9\f$, the outer disk is not relevant above
  0.2 kTbb.
- The alpha disk deviates by 2.5-111% for the same maximal temperature (and still by 5-10% for the best-fitting
  \f$kT_\mathrm{bb}\f$), due to the different temperature profile, see the figure above.

**Precision of `kerrbb`.** `kerrbb` interpolates tabulated spectra linearly in the spin (46 values, e.g., 0.904277 and
0.997676, but not 0.9 or 0.998) and the inclination (in steps of 5 deg), and integrates the energy bins with the
trapezoidal rule on 601 tabulated energies (spaced by 2.3%, which overestimates the flux by about 0.3% at 8 kTbb).
Compared directly with the tabulated spectra at spins of the grid (0, 0.904277, 0.997676), the converged reference of
the model (see (2a)) deviates by 0.4-1.4% in 0.5-10 kTbb, mostly at 5-10 kTbb for high spin and low inclination. This
corresponds to about 0.1% in the energy scale of the inner disk, which is about the precision of the relativistic
tables of both models.

Before, the full model was compared to the non-relativistic models `diskbb` and `ezdiskbb`. These comparisons were
limited by the relativistic effects (5-8% for \f$R_\mathrm{in} = 10\,r_\mathrm{isco}\f$ and 40-50 deg, where the
gravitational redshift and the Doppler shift partly cancel) and were replaced by the comparison with `kerrbb`.

### (3) Direct and reflected radiation from the same emitted primary

The reflection is calculated from the returning radiation `specRet[izone]`, which is calculated from the primary
of the zones `specPri[izone]`, whereas the direct emission is calculated with the kernel from the radial points (see
above). Direct and reflected radiation are therefore only consistent if `specPri` is the emission of exactly these
radial points, and if each zone covers the same area in the returning-radiation table, in the primary, and in its line
profile. The tests are done for the alpha disk and the Novikov-Thorne profile, for \f$a\f$ = 0 and 0.998
(\f$R_\mathrm{in} = r_\mathrm{isco}\f$, \f$R_\mathrm{out}\f$ = 1000 Rg, \f$kT_\mathrm{bb}\f$ = 1 keV):

- **(3a)** The primary of each zone is identical to the emission of the radial points binned into the zone (to
  \f$3 \times 10^{-15}\f$), and all zones together, weighted with their area, emit the same as the disk on a fine
  radial grid (5000 logarithmic zones): to 0.004-0.013% in the number of photons and to 0.002-0.004% in the energy.
- **(3b)** The returning radiation of each zone is identical to the returning radiation of this primary (to
  \f$4 \times 10^{-16}\f$).
- **(3c)** The primary of the zones, observed with the line profiles of the zones (incl = 30 deg), has the same number
  of photons (to 0.02-0.11%) and the same energy flux (to 0.05-0.17%) in 0.1-1000 keV as the direct emission of the
  model. The remaining difference is due to the variation of the temperature within a zone combined with the
  relativistic effects, which is not taken into account for the zones (largest for \f$a = 0.998\f$).

**Area of the zones.** The returning-radiation table assumes that each zone covers exactly its radial range
\f$[r_\mathrm{lo}, r_\mathrm{hi}]\f$. Each radial point represents the cell between the midpoints to its neighbors.
Before 2026-10, `calc_relline_profile()` assigned each cell completely to the zone the point is located in. The area of
a zone in its line profile therefore differed from the nominal zone by up to half a cell at each boundary:

| zone grid                                            | area difference per zone (rms / max) |
|------------------------------------------------------|--------------------------------------|
| relxillBB, 50 zones of the returning-radiation table | 0.9% / 4.0%                          |
| relxilllp, 10 zones (\f$h = 6\f$)                    | 0.9% / 1.9%                          |
| ionization gradient, 50 zones (\f$h = 6\f$)          | 3.2% / 7.9%                          |

As the emission per unit area changes strongly between the zones of relxillBB (by up to a factor of 7 in the number of
photons between the coarse outer zones), this is not consistent with the returning-radiation table. With the primary of
each zone binned in the
same way (complete cells), the zones together emit 0.5-0.9% fewer photons and 0.3% less energy than the disk, as
assumed for the returning radiation (fails (3a)). With the primary of each zone binned exactly, but complete cells in
the line profiles, the primary of the zones is observed with 0.6-1.2% more photons and 0.4-0.8% more energy than the
direct emission (fails (3c)), and the same applies to the reflection. Therefore, the cells are split at the zone
boundaries in proportion to their area, both for the primary and in `calc_relline_profile()`. For the other multi-zone
models, this is not relevant, as the spectra of neighboring zones (e.g., with a slightly different ionization) are very
similar: it changes them by at most \f$2 \times 10^{-4}\f$ (`relxilllpCp` with 10 zones: \f$4 \times 10^{-6}\f$).

## Unit tests

In `test/unit/tests-bbody-returnrad.cpp`; all disk emission tests run with `bin/tests "[diskemission]"`:

- **(1)** *" disk emission (1): temperature profile on a fine radial grid vs. Xspec diskbb and ezdiskbb "*: shape
  deviation below 0.3% (alpha disk) and 0.05% (diskbb), and normalization independent of kTbb to 0.1%.
- **(2a)** *" disk emission (2a): kernel, i.e., line profile with a line energy for each radial point "*: identical to
  the sum of the line profiles of the zones (to \f$10^{-10}\f$), and shifted by a constant line energy of 0.5 keV
  (number of photons and mean of \f$\ln E\f$ to \f$10^{-3}\f$).
- **(2b)** *" disk emission (2b): full relxillBB_NS model with the Novikov-Thorne profile vs. Xspec kerrbb "*
  (`[kerrbb]`): for \f$a\f$ = 0, 0.904277, 0.997676 (on the spin grid of `kerrbb`) and an inclination of 30 and 60
  deg, the shape deviation of the Novikov-Thorne profile is below 2% in 0.5-10 kTbb, and the one of the alpha disk
  above 2%.
- **(3a)** *" disk emission (3a): primary of the zones is the emission of the radial points of the line profiles "*:
  `specPri` of each zone identical to the emission of the radial points binned into the zone (to \f$10^{-10}\f$), and
  all zones together vs. the disk on a fine radial grid to 0.1% (photons and energy).
- **(3b)** *" disk emission (3b): returning radiation of the zones is calculated from this primary "*: `specRet`
  identical to the returning radiation of this primary (to \f$10^{-10}\f$).
- **(3c)** *" disk emission (3c): direct emission of the model vs. the primary of the zones observed "*: number of
  photons and energy flux in 0.1-1000 keV to 0.3%.

The tests (3) run for the alpha disk and the Novikov-Thorne profile, for \f$a\f$ = 0 and 0.998.

The reference spectra of the Xspec models are in `test/unit/refdata/`: `refvalues_diskbb_Tin{0.5,1.0,2.0}.fits`
(column `diskbb`) and `refvalues_ezdiskbb_Tmax{0.5,1.0,2.0}.fits` (column `ezdiskbb`), created with
`create_refvalues_diskbb.sl`, and `refvalues_kerrbb_a{0.000000,0.904277,0.997676}_i{30,60}.fits` (column `kerrbb`, with
the parameters in the header), created with `create_refvalues_kerrbb.sl`.

## History

- 2026-10: `disk_temperature_diskbb()` calculated \f$(kT_\mathrm{bb}\, r/R_\mathrm{in})^{-3/4}\f$ instead of
  \f$kT_\mathrm{bb} (r/R_\mathrm{in})^{-3/4}\f$, i.e., an inner temperature of \f$kT_\mathrm{bb}^{-3/4}\f$. This had
  no effect for kTbb = 1 keV, the only value that was tested.
- 2026-10: for the diskbb profile, kTbb was the temperature at the mean radius of the first zone instead of the inner
  edge of the disk, i.e., the disk was 1.5% hotter than `diskbb` (shape deviation of 8-9%). On the coarse energy grid
  of the unit test, this was partly compensated by the evaluation at the bin centres (3.7% in total).
- Both only affected `RELXILL_RRAD_BBODY_TPROFILE=2` and the test function `spec_diskbb()` (now
  \ref spec_multicolor_disk), not the default alpha-disk profile.
- 2026-10: the direct disk emission of each zone was calculated on sub-zones with
  \f$r_\mathrm{hi}/r_\mathrm{lo} \leq 1.1\f$ instead of a single blackbody per zone of the returning-radiation table
  (deviation from a fine radial grid 0.02% instead of 3.5% for the alpha disk, mainly below 0.3 kTbb, as the outer
  zones of the table are coarse, \f$r_\mathrm{hi}/r_\mathrm{lo}\f$ up to 2.4). The returning radiation was still
  calculated with a single blackbody per zone.
- 2026-10: the direct disk emission of each zone was calculated on the radial range the zone covered in its
  relativistic line profile instead of the nominal zone (up to 4% in area, see (3)).
- 2026-10: the observed direct emission is calculated for the whole disk with a single convolution of a kernel, with the
  temperature of each radial point, instead of folding the direct emission of each zone with the line profile of the
  zone. Close to the ISCO, the temperature and the relativistic effects change strongly within a zone: compared to the
  converged reference, the direct emission per zone deviated by 1.1% for \f$a = 0.998\f$ and 30 deg (mainly at 5-10
  kTbb, 0.02-0.17% otherwise), now 0.01-0.04%. The zones only contain the reflected emission.
- 2026-10: the primary of the zones `specPri`, from which the returning radiation is calculated, is the emission of the
  radial points of the line profiles binned into the zones, and the cells of the radial points are split at the zone
  boundaries, also in `calc_relline_profile()` (see (3)). Before, the returning radiation was calculated from a single
  blackbody per zone (which differs from the binned emission of the radial points by up to 7% in the number of photons
  and 18% in the energy of a zone, but only by up to 0.2% in the energy of the returning radiation), and the zones of
  the line profiles were composed of complete cells. The reflected spectrum changed by -0.03% to -0.34% in the energy
  flux. The sub-zones (\ref get_disk_emission_zones_binned) are only used by \ref spec_multicolor_disk.

## Open issues

- The reflection is calculated on the 50 zones of the returning-radiation table, and each zone is convolved with its
  line profile as a whole. The variation of the temperature within a zone combined with the relativistic effects is
  therefore not taken into account for the reflection (0.05-0.17% in the energy flux of the primary, see (3c)).
- The kernel of the direct emission takes about as long as the line profiles of the zones (up to 30 ms for high
  inclinations), see Runtime. If needed, it could be calculated on fewer radial points at large radii.
- \f$R_\mathrm{max} = 1.5^{4/5} R_\mathrm{in}\f$ instead of \f$\frac{49}{36} R_\mathrm{in}\f$: the maximal temperature
  is 0.03% higher than \f$kT_\mathrm{bb}\f$ (0.2% in the spectrum at 10 kTbb, see (1)).
- Which temperature profile should be the default: the alpha disk is Newtonian, whereas the relativistic effects are
  calculated in the Kerr metric (see the Novikov-Thorne profile above).
- The parameter `shiftTmaxRRad` of `relxillBB_NS` is read, but not used in the code.
