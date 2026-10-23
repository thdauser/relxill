# Xillver Normalization — Working Notes (RTR / `relxillBB`)

Note: may be partially outdate from the large changes done in 21.10.2026

**The derivation of the normalization lives in the Doxygen page
`doc/returnrad_normalizations.md` (`\ref returnrad_normalizations`), which is
the authoritative description.** This file holds only what that page does not:
operational details needed to run the manual tests in this directory, read the
`__*.dat` / `__*.fits` diagnostic output, and interpret older plots.

---

## 1. What changed (read this before trusting `plots/` or `old/`)

The normalization was simplified substantially in mid-2026. Plots and notes
produced before then used a scheme that no longer exists:

| Removed                                                                                                                                                                         | Commit                 |
|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------|
| `normfacMatchAtHighEnergies` — rescaled the reflected continuum to match the returning blackbody above `4·kTbb`; found to be wrong                                              | `6f07241` (2026-07-01) |
| `calc_angle_weighted_xillver_spec`, and with it the `xill_norm_avg / xill_norm` factor that re-anchored the angle-weighted spectrum to the angle-*averaged* table normalization | `df34569` (2026-06-29) |

Symbols that appear in older notes and are **gone**: `getAnglecorrXillSpec`,
`calcXillverNormFromPrimarySpectrum`, `get_avg_xillver_spec_normalization`,
`calc_angle_averaged_xill_spec`, `calculateAngleWeightedXillverSpectrum`,
`xill_norm_avg`, `normfacMatchAtHighEnergies`.

Both removals have the same root cause: each set the reflected level by a
criterion *other than* the table's own Dauser+2016 calibration, i.e. both
altered the **relative** (reflection-to-primary) normalization. Neither was
about the absolute scale — that was already correctly set by `N_xill` from the
returning radiation.

The reflected term now carries exactly two factors, `|boost|` and
`N_xill` — see the Doxygen page. **Any plot in `plots/` or `old/` predating
2026-07-01 shows a different normalization and is not comparable to current
output.**

---

## 2. Units and grids (needed to read the diagnostic files)

`returnSpec->specRet` / `specPri` are **bin-integrated photons/bin** (Xspec
standard), documented at `Relreturn_Datastruct.h:33`. Not photons/cm²/s/keV —
dividing by bin width is a common mistake when plotting these directly.

Three different energy grids are in play, which is why the diagnostic FITS
files do not all share an x-axis:

| Grid                                                                                   | Bins                                   | Range                | Where                            |
|----------------------------------------------------------------------------------------|----------------------------------------|----------------------|----------------------------------|
| Convolution grid — `returnSpec->ener`, and what `__*.fits` zone spectra are written on | `N_ENER_CONV = 4096` (power of 2, FFT) | `0.00035 – 2000` keV | `get_relxill_conv_energy_grid()` |
| Coarse xillver normalization grid                                                      | `N_ENER_COARSE = 500`                  | `0.1 – 1000` keV     | `get_coarse_xillver_energrid()`  |
| Xspec input grid                                                                       | model-dependent                        | model-dependent      | caller                           |

The convolution grid fully covers the `0.1–1000` keV normalization band, so the
coverage guard inside `calc_xillver_norm_from_primary_spectrum` always passes
for spectra on `returnSpec->ener`. Noted here so it is not re-investigated as a
suspected bug.

---

## 3. Environment toggles

Both default to `0` and are read via `envoptions.c`; `print_env_options()`
reports their state.

| Variable                      | Effect                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|-------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `RELXILL_RRAD_BBODY_MIRRORED` | Replaces the xillver reflection with a plain mirror of `specRet` (`calc_zone_mirrored_intensity`). Test path only — applies `\|boost\|`, but **no** `N_xill`; the returning flux (which already contains the angles of incidence) is re-emitted with the same spectral shape and an isotropic intensity, `specRet / 0.25` (since 2026-10-23, see `doc/xillver_angular_dependence.md`). The direct emission is calculated for the whole disk, as for the reflection. |
| `RELXILL_RRAD_BBODY_IONGRAD`  | Enables the per-zone ionization gradient (`set_xill_param_from_iongrad` overwrites `dens`, `lxi`, `radzone` from `ion_gradient`).                                                                                                                                                                                                                                                                                                                                   |

Diagnostic `__*.dat` / `__*.fits` output is gated separately by
`shouldOutfilesBeWritten()`. Note the output **filenames differ** between the
mirrored and normal paths for the reflected spectra of the zones
(`obs-mirror-reflect.fits` vs `obs-reflect.fits`, `write_bbret_diagnose_outfiles`),
so a plot script hard-coding one set will silently find nothing when the other
toggle is active. The primary of the zones (`rframe-bbodyEmit.fits`,
`obs-bbodyEmit.fits`) is the same in both paths. The list of all diagnostic files
is in `doc/returnrad_disk_emission.md` ("Diagnostic output files", with the content and the flux conversion of each
file).

---

## 4. Open issues

### 4.1 `radzone` is silently clamped to the table boundary

`relxillbb_kernel` sets `xill_param_zone.radzone = izone` with the in-code
warning *"very dangerous; need to check that it is within the table boundaries
as well"* (`Relreturn_BlackBody.cpp:992-994`). The concern is real and still open.

The path is `radzone` → `radindex` (`xilltable_reading.c:667`) →
`global_param_vals[PARAM_RID]` → the generic interpolation loop in
`Xilltable_handling.cpp:640`, which calls `resetInpvalsToBoundaries()`. So an
out-of-range zone index does **not** crash or error — it is **clamped to the
last tabulated radial zone**, and the warning only prints under
`is_debug_run()`.

The loop bound is `relProfile->n_zones`, derived from the return-radiation
radial grid, which has no enforced relationship to the `radzone` dimension of
the xillverRTR table. If the number of zones exceeds the table's, every zone
past the end collapses onto the same tabulated reflection and the model stays
silent about it. **Worth an explicit check, and worth testing with
`RELXILL_NUM_RZONES` set high.**

### 4.2 Incident-shape mismatch — and why it must *not* be fixed by rescaling

Two things set the reflected flux, and they must not be mixed up (see
"Relative and absolute normalization" in the Doxygen page):

- **Relative** — the reflection-to-primary ratio, calibrated into the table via
  the Dauser+2016 primary-spectrum integral. This is a property of the table
  and must not be changed by any other means.
- **Absolute** — set by the returning radiation that actually irradiates the
  zone. `N_xill` is that returning primary expressed in table units, which is
  how the absolute scale legitimately enters.

The shape question below concerns the **relative** side: whether the table's
reflection-to-primary ratio is the right one for the illumination at hand. It
is never a reason to adjust the absolute scaling.

This matters here because the two RTR models sit differently with respect to
the incident *shape*:

| Model                                         | Table                  | Incident spectrum assumed                                   |
|-----------------------------------------------|------------------------|-------------------------------------------------------------|
| `MOD_TYPE_RELXILLBBRET` (`relxillBB`)         | `xillverNS-2.fits`     | a **single** blackbody at `kTbb`, same for every zone       |
| `MOD_TYPE_RELXILLBBRETXILL` (`relxillBBxill`) | `xillverRTR_v0.3.fits` | per-radial-zone, selected by the `radindex` table dimension |

(`get_xilltable_id`, `xilltable_reading.c:640-647`.)

The actual returning radiation is a superposition of blue- and red-shifted
blackbodies from across the disk, with an effective temperature that varies by
zone. For `relxillBB` the table cannot represent that, so the single-`kTbb`
assumption stated in the function header is a genuine limitation of that model.

**The remedy is a table computed for the correct incident spectrum, never a
correction factor.** `xillverRTR`'s `radindex` dimension is precisely that
remedy, applied properly. The removed `normfacMatchAtHighEnergies` (`6f07241`)
was an attempt to patch the same mismatch by rescaling, and that is exactly why
it was wrong.

So: quantifying how far the effective temperature of `specRet[izone]` drifts
from `kTbb` across the disk is a **useful diagnostic** — it bounds the error of
the `relxillBB` approximation and indicates when `relxillBBxill` is required.
It is *not* an input to any correction factor, and no such factor should be
added on the basis of it. See also §4.1: with `xillverRTR`, silent clamping of
`radindex` would defeat the per-zone tabulation and reintroduce the mismatch
without warning.

## 5. Suggested cross-check

For a reflection-only run (`boost < 0`), the xillver path and the mirrored path (`RELXILL_RRAD_BBODY_MIRRORED=1`) can be
compared, as both hand an **intensity**
to the relativistic convolution (see `doc/xillver_angular_dependence.md`): the
xillver path the angle-mixed table intensity scaled by `N_xill`, the mirrored
path the isotropic intensity `specRet / 0.25` (no xillver renormalization).

- In the rest frame (`rframe-*`), all files are fluxes (since 2026-10-23,
  see "Diagnostic output" in `doc/xillver_angular_dependence.md`). The mirror
  re-emits exactly the returning flux: `rframe-mirror-reflect = |boost| ×
  rframe-bbodyRet`. The xillver reflection (`rframe-reflect`, angle-averaged
  with the slab weights) has 0.36-0.46 of the returning energy flux per zone,
  0.43 for the whole disk (`B_slab` of `[xillverRTR][angleweight]`).
- Observed (`obs-*`), both are converted to a flux by the same convolution, so
  the ratio of the integrals of `obs-reflect` and `obs-mirror-reflect` is a
  test that `N_xill` reproduces the incident flux level. The ratio is not
  0.43 (`B_slab`), but the angle-mixed `A_mix/4`, and depends on the
  inclination and the energy band: 0.32-0.47 for 20-80 deg in 0.1-1000 keV (see the table in
  `doc/xillver_angular_dependence.md`, and the unit test
  `[mirror]`, which checks that it is below 1).

Diagnostic files and plots made before 2026-10-23 are not comparable:
`obs-mirror-reflect` and `obs-bbodyEmit` were a factor 4 lower,
`rframe-reflect` was the angle-mixed intensity for the observer inclination.

Note this compares *integrals*, not shapes — the mirrored path has no
reflection physics at all.

If the two disagree, the conclusion is that `N_xill`, the illuminating spectrum
fed to it, or the table's relative normalization is wrong. It is **not** grounds for introducing a
correction factor (§4.2).

The direct term of the model is handed to the convolution as the intensity
`specPri / 0.25`, as the reflected paths (since 2026-10-23, before as the flux,
i.e., a factor 4 too weak). It agrees with the sum of `obs-bbodyEmit` over the
zones, and its normalization relative to the mirror is tested by the unit test
(3d) in `tests-returnrad-bbody-diskemission.cpp`.
