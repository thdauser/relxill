# Xillver Normalization — Working Notes (RTR / `relxillBB`)

State of 2026-10-23 (checked against the code).

**The authoritative descriptions are the Doxygen pages** `doc/returnrad_normalizations.md` (normalization of the
reflection, `N_xill`), `doc/xillver_angular_dependence.md` (intensity vs. flux, angle weighting, direct emission), and
`doc/returnrad_disk_emission.md` (disk emission, diagnostic output files, tests). This file only holds what these pages
do not: the history of the normalization, units and grids for reading the diagnostic output, and open issues.

---

## 1. History (read this before trusting old plots or notes)

| Change                                                                                                                                                                                  | Commit / date          |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------|
| removed `calc_angle_weighted_xillver_spec`, and with it the `xill_norm_avg / xill_norm` factor that re-anchored the angle-weighted spectrum to the angle-*averaged* table normalization | `df34569` (2026-06-29) |
| removed `normfacMatchAtHighEnergies`, which rescaled the reflected continuum to match the returning blackbody above `4·kTbb`                                                            | `6f07241` (2026-07-01) |
| direct emission of the whole disk with one convolution, primary of the zones from the same radial points, cells split at the zone boundaries (`changes_2026-10-21_disk_emission.md`)    | 2026-10-21             |
| mirror and direct emission handed to the convolution as intensities (`F / 0.25`), as the xillver reflection; all diagnostic files are fluxes (`doc/xillver_angular_dependence.md`)      | 2026-10-23             |

Both removals of mid-2026 have the same root cause: each set the reflected level by a criterion *other than* the
table's own Dauser+2016 calibration, i.e., both altered the **relative** (reflection-to-primary) normalization. The
reflected term now carries exactly two factors, `|boost|` and `N_xill`. Symbols that appear in older notes and are
**gone**: `getAnglecorrXillSpec`, `calcXillverNormFromPrimarySpectrum`, `get_avg_xillver_spec_normalization`,
`calc_angle_averaged_xill_spec`, `calculateAngleWeightedXillverSpectrum`, `xill_norm_avg`, `normfacMatchAtHighEnergies`,
`calc_zone_xillver_reflected_flux` and `calc_zone_mirrored_flux` (now `..._intensity`).

**Since 2026-10-23, the reflected-to-direct ratio is a factor 4 smaller than before** (the direct emission was a factor
4 too weak). Diagnostic files and plots made before are not comparable.

---

## 2. Units and grids (needed to read the diagnostic files)

`returnSpec->specRet` / `specPri` are **bin-integrated photons/bin** (Xspec standard, documented in
`Relreturn_Datastruct.h`), per unit area of the zone. Not photons/cm²/s/keV: dividing by the bin width is a common
mistake when plotting them directly.

| Grid                                                                                     | Bins                                   | Range                | Where                            |
|------------------------------------------------------------------------------------------|----------------------------------------|----------------------|----------------------------------|
| Convolution grid — `returnSpec->ener` in the model, and the grid of the `__*.fits` files | `N_ENER_CONV = 4096` (power of 2, FFT) | `0.00035 – 2000` keV | `get_relxill_conv_energy_grid()` |
| Coarse xillver normalization grid (`N_xill`)                                             | `N_ENER_COARSE = 500`                  | `0.1 – 1000` keV     | `get_coarse_xillver_energrid()`  |
| Internal grid of the blackbody of the zones                                              | 500                                    | `0.001 – 500` keV    | `get_std_bbody_energy_grid()`    |
| Xspec input grid                                                                         | model-dependent                        | model-dependent      | caller                           |

The convolution grid fully covers the `0.1–1000` keV normalization band, so the coverage guard inside
`calc_xillver_norm_from_primary_spectrum` always passes for spectra on `returnSpec->ener`. Noted here so it is not
re-investigated as a suspected bug.

---

## 3. Environment toggles

All environment variables are listed in `doc/environment_variables.md`. The ones specific to the diagnostics:

| Variable                      | Effect                                                                                                                                                                                                                                                                                                               |
|-------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `RELXILL_RRAD_BBODY_MIRRORED` | Replaces the xillver reflection with a mirror of `specRet` (`calc_zone_mirrored_intensity`). Test path only — applies `\|boost\|`, but **no** `N_xill`; the returning flux (which already contains the angles of incidence) is re-emitted with the same spectral shape and an isotropic intensity, `specRet / 0.25`. |
| `RELXILL_RRAD_BBODY_ORDER`    | Number of orders of the returning radiation; the model returns the **sum** of all orders up to this value.                                                                                                                                                                                                           |
| `RELXILL_RRAD_BBODY_IONGRAD`  | Enables the per-zone ionization gradient (`set_xill_param_from_iongrad` overwrites `dens`, `lxi`, `radzone` from `ion_gradient`).                                                                                                                                                                                    |
| `RELXILL_WRITE_FILES`         | Writes the diagnostic `__*.dat` / `__*.fits` files (`shouldOutfilesBeWritten()`).                                                                                                                                                                                                                                    |

The output **filenames differ** between the mirrored and the xillver path for the reflected spectra of the zones
(`obs-mirror-reflect.fits` vs `obs-reflect.fits`, `write_bbret_diagnose_outfiles`). The primary of the zones
(`rframe-bbodyEmit.fits`, `obs-bbodyEmit.fits`) and the returning radiation (`rframe-bbodyRet.fits`) are the same in
both paths. The list of all files, their content and flux conversion is in `doc/returnrad_disk_emission.md`
("Diagnostic output files").

---

## 4. Open issues

### 4.1 `radzone` is silently clamped to the table boundary

`relxillbb_kernel` sets `xill_param_zone.radzone = izone` for each zone of the returning radiation, with the in-code
warning *"very dangerous; need to check that it is within the table boundaries as well"*. The path is `radzone` →
`radindex` (`get_xilltab_param()`, `xilltable_reading.c`) → `global_param_vals[PARAM_RID]` → the generic interpolation
in `Xilltable_handling.cpp`, which calls `resetInpvalsToBoundaries()`. An out-of-range zone index does **not** crash
or error: it is **clamped to the last tabulated radial zone**, and the warning is only printed in a debug run.

The zones of `relxillBB` are the 50 zones of the returning-radiation table (`get_rrad_fractions()`; not affected by
`RELXILL_NUM_RZONES`, which only sets the zones of the lamp-post and ionization-gradient models). Their number
decreases if `Rin` is increased, and the radii of the zones depend on the spin. For the current mock table (a = 0.998,
`Rin` at the ISCO, see §4.2), the unit test *"xillverRTR table: number of radial zones equals the radial grid of the
returning radiation"* checks that the table (`radindex`: 50 values) matches. A table with a different zone grid would be
clamped or mismatched without an error: the model should check the zone grid of the table explicitly. This becomes
essential for the final table with a spin dimension, where `radindex` refers to the zone grid of each spin (see
`doc/xillver_RTR.md`).

### 4.2 Incident spectrum of the table — and why a mismatch must *not* be fixed by rescaling

Two things set the reflected flux, and they must not be mixed up (see `doc/returnrad_normalizations.md`):

- **Relative** — the reflection-to-primary ratio, calibrated into the table via the Dauser+2016 primary-spectrum
  integral. This is a property of the table and must not be changed by any other means.
- **Absolute** — set by the returning radiation that actually irradiates the zone. `N_xill(specRet)` is that returning
  flux expressed in table units, which is how the absolute scale legitimately enters.

The two models sit differently with respect to the incident *shape* (`get_xilltable_id()`, `xilltable_reading.c`):

| Model                                          | Table                  | Incident spectrum of the table                                                   |
|------------------------------------------------|------------------------|----------------------------------------------------------------------------------|
| `relxillBB` (`MOD_TYPE_RELXILLBB`, -301)       | `xillverRTR_v0.3.fits` | `specRet` of each radial zone, selected by the `radindex` dimension (50 zones)   |
| `relxillBB_NS` (`MOD_TYPE_RELXILLBB_NS`, -300) | `xillverNS-2.fits`     | a **single** blackbody at `kTbb`, the same for every zone (kept for development) |

- `relxillBB`: the current table `xillverRTR_v0.3.fits` is a **mock version**, calculated for a single parameter set
  (a = 0.998, kTbb = 2 keV, logxi = 2, logN = 18; dimensions only `radindex` and the inclination). Until the final
  table is available, the model fixes these parameters (and Rin = ISCO, Rout = 1000; `check_parameter_bounds_relxillbb`
  and the parameter limits in `lmodel_relxill_devel.dat`), so only the inclination, z, boost, and the normalization
  are free. For these parameters, the incident spectrum of the table and `specRet` are identical by design. **The final
  table will in addition depend on kTbb, logxi, logN, and the spin**, with the incident spectrum of each zone
  calculated by the model for each combination of kTbb and spin.
- `relxillBB_NS`: the returning radiation is a superposition of blue- and red-shifted blackbodies from across the
  disk, with an effective temperature that varies between the zones. A single blackbody at `kTbb` can not represent
  that; this is a limitation of the model.

**The remedy for a shape mismatch is a table calculated for the correct incident spectrum, never a correction factor.**
The `radindex` dimension of `xillverRTR` is exactly that remedy. The removed `normfacMatchAtHighEnergies` (`6f07241`)
was an attempt to patch the same mismatch by rescaling, and that is why it was wrong. The final table therefore covers
kTbb and the spin with incident spectra calculated for each value, instead of rescaling the spectra of the mock
table.

---

## 5. Suggested cross-check

For a reflection-only run (`boost < 0`), the xillver path and the mirrored path (`RELXILL_RRAD_BBODY_MIRRORED=1`) can be
compared, as both hand an **intensity** to the relativistic convolution (see `doc/xillver_angular_dependence.md`): the
xillver path the angle-mixed table intensity scaled by `N_xill`, the mirrored path the isotropic intensity
`specRet / 0.25` (no xillver renormalization).

- In the rest frame (`rframe-*`), all files are fluxes. The mirror re-emits exactly the returning flux:
  `rframe-mirror-reflect = |boost| × rframe-bbodyRet`. The xillver reflection (`rframe-reflect`, angle-averaged with
  the slab weights) has 0.36-0.46 of the returning energy flux per zone, 0.43 for the whole disk (`B_slab` of
  `[xillverRTR][angleweight]`).
- Observed (`obs-*`), both are converted to a flux by the same convolution. The ratio of the integrals of
  `obs-reflect` and `obs-mirror-reflect` is not 0.43 (`B_slab`), but the angle-mixed `A_mix/4`, and depends on the
  inclination and the energy band: 0.32-0.47 for 20-80 deg in 0.1-1000 keV (see `doc/xillver_angular_dependence.md`,
  and the unit test `[mirror]`, which checks that it is below 1).

This compares *integrals*, not shapes — the mirrored path has no reflection physics at all. If the two disagree, the
conclusion is that `N_xill`, the illuminating spectrum fed to it, or the table's relative normalization is wrong. It is
**not** grounds for introducing a correction factor (§4.2).

The direct emission is handed to the convolution as the intensity `specPri / 0.25`, as the reflected paths. Its
normalization relative to the mirror is tested by the unit test (3d) in `tests-returnrad-bbody-diskemission.cpp`
(see `doc/returnrad_disk_emission.md`).
