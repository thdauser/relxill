# Changes of 2026-10-21: disk emission and returning radiation of relxillBB

Summary of the 15 commits of 2026-10-21 (`8e84579` to `b212378`), of which 10 change the model code. The
documentation is in `doc/returnrad_disk_emission.md` and `doc/returnrad_normalizations.md`.

After these changes, relxillBB and relxillBB_NS

- calculate the direct disk emission for the whole disk on the 1000 radial points of the line profiles,
- base the reflection on exactly the same emitted primary,
- use line-profile zones that cover exactly their radial range.

In addition, two bugs of the diskbb temperature profile are fixed, and a Novikov-Thorne temperature profile is new.

## Model code changes

### 1. diskbb temperature profile (`8e84579`, `ad7e6c7`)

Only affects `RELXILL_RRAD_BBODY_TPROFILE=2`, not the default alpha disk.

- `disk_temperature_diskbb()` calculated (kTbb·r/Rin)^(−3/4) instead of kTbb·(r/Rin)^(−3/4). It was only correct
  for kTbb = 1 keV.
- kTbb now refers to the inner edge `rlo[0]`, not to the mean radius of the first zone. Before, the disk was 1.5% too
  hot.

### 2. Novikov-Thorne temperature profile (`da5bfdf`)

Selected with `RELXILL_RRAD_BBODY_TPROFILE=3`.

- New `TPROFILE_NOVIKOV_THORNE`: the closed form of Page & Thorne (1974), as in the Xspec model `kerrbb`.
- R_max is found by a golden-section search, and Rin must be at or outside the ISCO.
- `get_tprofile()` now takes the spin. The default stays the alpha disk.

### 3. Direct emission of the disk (`09df5b0` → `515d3e6` → `505216a`)

- **Final state:** `calc_disk_emission_observed()` is one FFT convolution of the blackbody at kTbb with a kernel. The
  kernel is the line profile of each of the 1000 radial points at the line energy t(r), weighted with t(r)³
  (`get_disk_emission_kernel()`, new `calc_relline_profile_radial_line_energy()` in `Relprofile.cpp`).
- The kernel does not depend on kTbb. It is cached by the relativistic parameters plus the temperature profile.
- Blackbody bins that the cyclic FFT would shift off the energy grid are set to zero (only matters below about
  0.05 keV).
- The number of photons is normalized over the whole energy grid, not only in 0.01-1000 keV, which would have
  overestimated the flux by up to 1%.
- The zones now only carry the reflection, and for boost = 0 they are skipped. `calc_zone_direct_plus_*` were renamed
  to `calc_zone_xillver_reflected_flux` and `calc_zone_mirrored_flux`.
- The intermediate steps of the morning (sub-zones in `09df5b0`, zone ranges of the line profiles in `515d3e6`) are
  superseded and were removed from the model in `817cc9e`. The sub-zones remain only in the test helper
  `spec_multicolor_disk`.

### 4. Primary of the zones and returning radiation (`817cc9e`)

- `get_disk_emission_zones()` bins the emission of the same 1000 radial points into the 50 zones of the
  returning-radiation table. The cell of each point is split at the zone boundaries in proportion to its area. This is
  `specPri`.
- `get_returnrad_specs()` now takes these spectra instead of the zone temperatures, so `specRet` is calculated from
  exactly this primary.
- The single-temperature `get_bbody_specs()` and `get_zone_ranges_relline_profile()` are removed.
  `get_temperature_radial_points()` is shared by the kernel and the zones.

### 5. `calc_relline_profile()` splits the cells at the zone boundaries (`817cc9e`)

Affects all multi-zone models.

- Each zone now covers exactly [r_lo, r_hi]. Before, each cell was assigned completely to one zone, which changed the
  area of a zone by up to 4% (relxillBB) and 8% (ionization gradient).
- This is needed for the direct and the reflected radiation to agree (see test (3c)). It is in `CHANGELOG.txt` under
  "next version".

### 6. Diagnostic output (`b212378` and the follow-up commit)

Only with `RELXILL_WRITE_FILES=1`; the model output is unchanged.

- The observed primary of each zone (`specPri` convolved with the line profile of the zone) is written to its own file,
  `obs-bbodyEmit.fits`, for any boost (the rest-frame primary is in `rframe-bbodyEmit.fits`, as before).
- The zone files only contain what the zones contain in the model, i.e., the reflected spectrum scaled with |boost|:
  `rframe-reflect.fits` and `obs-reflect.fits` (`rframe-mirror-reflect.fits` and `obs-mirror-reflect.fits` in the mirror
  mode), written only for boost ≠ 0. The combined (boost = 1) and primary (boost = 0) zone files are removed; the total
  spectrum of a zone is the sum of its primary and its reflected spectrum.
- In `b212378`, the direct emission was first added to the zone files instead (`add_direct_emission_zones_diagnose`);
  this was replaced by the separate file.
- The scripts in `test/manual/returnrad` read the new file (`fobs.primary`), need one model evaluation per mode
  (boost = −1) instead of three, and the mirror plot shows the mirrored reflection and the primary per zone without
  the combined curves.

### Net effect on the output and the runtime

|                                               | change                                                                                                          |
|-----------------------------------------------|-----------------------------------------------------------------------------------------------------------------|
| direct emission vs. converged reference       | 0.01-0.04% (before: up to 1.1% at a = 0.998, 30°, near the ISCO)                                                |
| direct emission vs. `kerrbb` (Novikov-Thorne) | 0.2-1.4%, which is the precision of `kerrbb` itself                                                             |
| reflection of relxillBB                       | −0.03% to −0.34% in the energy flux                                                                             |
| other multi-zone models                       | ≤ 2e-4 (`relxilllpCp`: ≤ 4e-6)                                                                                  |
| runtime                                       | +6% (30°) to +11% (60°) if a relativistic parameter changes; unchanged if only kTbb changes; −12% for boost = 0 |

## Tests that verify these changes

### Unit tests

In `test/unit/tests-returnrad-bbody-diskemission.cpp`, run with `bin/tests "[diskemission]"`:

| test                                                                                       | covers                                                                                                                                                        | how                                                                                                                                                                                                             |
|--------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **(1)** profile vs. Xspec `diskbb`/`ezdiskbb`                                              | 1 (both diskbb bugs), alpha disk, kTbb normalization                                                                                                          | fine radial grid, shape within 0.05% (diskbb) and 0.3% (alpha); normalization independent of kTbb to 0.1%                                                                                                       |
| **(2a)** kernel                                                                            | 3: `calc_relline_profile_radial_line_energy`                                                                                                                  | at 1 keV equal to the sum of the 50 zone profiles (2e-15), i.e., same points, areas, relativistic transfer, and photon conversion; a line energy of 0.5 keV shifts it correctly and keeps the number of photons |
| **(2b)** full `relxillBB_NS` (boost = 0) vs. `kerrbb`                                      | 2 (Novikov-Thorne), 3 end to end (kernel, FFT, boost = 0 path), kernel cache per temperature profile                                                          | 3 spins × 2 inclinations, shape < 2% in 0.5-10 kTbb; the alpha disk must fail, which also catches a cache that ignores the temperature profile (both profiles run with identical parameters)                    |
| **(3a)** primary of the zones                                                              | 4: `get_disk_emission_zones`, cell split                                                                                                                      | `specPri` equals an independent binning of the 1000 points (1e-10); all zones × nominal area vs. a disk of 5000 zones to 0.1% (measured ≤ 0.013%)                                                               |
| **(3b)** returning radiation                                                               | 4: `specRet` is calculated from that primary                                                                                                                  | equals `get_returnrad_specs` applied to the primary of the test (1e-10)                                                                                                                                         |
| **(3c)** direct emission vs. zones observed                                                | 3 + 4 + 5: consistent normalization of the direct and reflected radiation, the cell split in `calc_relline_profile`, the absolute normalization of the kernel | photons and energy in 0.1-1000 keV to 0.3% (measured 0.02-0.17%); without the cell split it fails by 0.4-1.2%                                                                                                   |
| *"relxillBB: boost parameter …"* (`tests-returnrad-bbody-model.cpp`)                       | 3: boost logic                                                                                                                                                | the direct emission is added exactly once                                                                                                                                                                       |
| *"xillverRTR normalization …"* (`[xillverRTR][norm]`, `tests-returnrad-bbody-xillver.cpp`) | the renamed `calc_zone_xillver_reflected_flux`                                                                                                                | normalization of the reflection per zone                                                                                                                                                                        |
| full suite (125 test cases)                                                                | 5: no regression of the other models                                                                                                                          |                                                                                                                                                                                                                 |

### Manual tests

In `test/manual/returnrad` (ISIS, `make test` and `make plot`). On 2026-10-23, the two tests were moved to the unit
tests ("disk emission (1b)" and the reference-data test of the returning spectra, see `test/manual/returnrad/README.md`),
and `make test` was removed:

- **`testEmissivityBbody`:** checks item 4 against an analytic result. The energy flux of each zone agrees with SS73
  averaged over the area of the zone to 5.8e-4.
- **`testRframeReturnradRefdata`:** regression test of the summed `specPri` and `specRet` against the newly generated
  reference data (now tracked in git)
- **Diagnose plots:** visual check of item 6. Checked by hand, the 50 zone primaries (`obs-bbodyEmit.fits`) add up to
  the direct spectrum of the model to 0.06% in photons and 0.13% in energy.

### Not covered by a dedicated test

- **Cache invalidation if a relativistic parameter changes:** the cache key on the temperature profile is covered by
  (2b), but a change of, e.g., the inclination is only covered indirectly.
- **The FFT wrap-around protection:** only checked by hand against a non-cyclic prototype (identical above 0.05 keV).
- **The effect of the cell split on the other models:** only measured by hand for `relxilllpCp`; the ionization-gradient
  models were not re-measured.
- **The diskbb profile in the full model and in the tests (3):** only the alpha disk and the Novikov-Thorne profile are
  used there.
