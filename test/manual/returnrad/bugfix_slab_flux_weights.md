# Bug fix: inclinations converted to degrees twice in `get_semi_infinite_slab_flux_weights()`

## Summary

`get_semi_infinite_slab_flux_weights()` (`src/Xillspec.cpp:125`) computed the semi-infinite slab flux weights
w_i = 0.5 cos(i) / n_incl from the xillver table inclinations as

```cpp
weights[ii] = norm_factor_semi_infinite_slab(xill_spec->incl[ii] * 180 / M_PI) / xill_spec->n_incl;
```

`xillSpec::incl` is already given in **degrees** (copied from the inclination grid of the xillver table,
`Xilltable_handling.cpp:623`), and `norm_factor_semi_infinite_slab()` expects degrees. The additional
`* 180 / M_PI` therefore evaluates the cosine at meaningless angles:

| table incl [deg] | 18.2  | 31.8  | 41.4   | 49.5  | 56.6  | 63.3  | 69.5  | 75.5  | 81.4  | 87.1  | sum   |
|------------------|-------|-------|--------|-------|-------|-------|-------|-------|-------|-------|-------|
| weight (bug)     | 0.040 | 0.047 | -0.042 | 0.035 | 0.050 | 0.046 | 0.046 | 0.050 | 0.048 | 0.034 | 0.351 |
| weight (correct) | 0.048 | 0.043 | 0.038  | 0.033 | 0.028 | 0.023 | 0.018 | 0.013 | 0.008 | 0.003 | 0.250 |

The correct weights sum to 0.25 (the table inclinations are the bin centres of a grid uniform in cos(i),
i.e. the sum is the integral of 0.5 mu over mu from 0 to 1). The buggy weights are too large by about 40% in total,
contain a negative weight, and lack the cos(i) dependence (limb weighting).

The bug is also present on `master` (it predates the refactoring in `c8ed2df`, which only moved the line).

## Affected code paths

`get_semi_infinite_slab_flux_weights()` is only used through `calc_angle_averaged_xillver_flux()`, which in the
model is only used to calculate the correction factors of the **coronal returning radiation**:

1. `src/Relxill.cpp:376-378` -- `relxill_kernel()` computes `rel_param->rrad_corr_factors = calc_rrad_corr_factors(...)`
   if `return_rad != 0` and `a > SPIN_MIN_RRAD_CALC_CORRFAC` (= 0).
2. `src/Relxill.cpp:174-193` -- `calc_rrad_corr_factors()` calls `get_xillver_fluxcorrection_factors()` per zone.
3. `src/Xillspec.cpp:545` -- `get_xillver_fluxcorrection_factors()` calls `calc_angle_averaged_xillver_flux()` (`:147`).
4. `src/Xillspec.cpp:122-128` -- `get_semi_infinite_slab_flux_weights()`, **the bug (line 125)**.
5. `src/Xillspec.cpp:385-403` -- `get_xillver_fluxcorr()`: `corrfac_flux` = energy flux (0.1-1000 keV) of the
   angle-averaged reflection / energy flux of the primary. It is linear in the weights.

Where the correction factors act:

- `src/Relreturn_Corona.cpp:162` -- `emis_return->emis[i] *= corrfac_flux[i]` (returning emissivity, directly affected).
- `src/Relreturn_Corona.cpp:151`, `:110-111` -- `corrfac_gshift` is a ratio of two fluxes of the same spectrum, so the
  normalization of the weights cancels and only the angular shape enters (weakly affected).
- `src/Rellp.cpp:585-586` -- `add_returnrad_emis()` (called at `:644`) adds the returning emissivity to the
  emissivity of the primary source.

The returning radiation is switched on by default for the LampPost, Ring and Slab irradiation geometries
(`src/ModelDefinition.cpp:150-161`, `$switch_returnrad 1` in `lmodel_relxill_public.dat`). Therefore, **all lamp-post
and extended-corona models with returning radiation and spin > 0 are affected**, e.g. `relxilllp` with default
parameters.

Not affected: models without returning radiation, the `xillver` models (they use `norm_xillver_spec()` with the
inclination parameter in degrees, not `xillSpec::incl`), and the `relxillBB` kernel (it uses the GR emission-angle
distribution, not the slab weights).

## Size of the effect

Measured by linking a test program once against the unchanged library and once with only line 125 fixed
(`relxilllp`, default parameters unless stated otherwise).

Flux correction factors (as in `test/unit/tests-returnrad.cpp`):

|                               | bug   | fixed | bug / fixed |
|-------------------------------|-------|-------|-------------|
| `corrfac_flux`, logxi = 3.0   | 0.812 | 0.529 | 1.53        |
| `corrfac_flux`, logxi = 0.0   | 0.236 | 0.149 | 1.58        |
| `corrfac_gshift`, logxi = 3.0 | 1.080 | 1.092 | 0.99        |
| `corrfac_gshift`, logxi = 0.0 | 0.699 | 0.697 | 1.00        |

Contribution of the returning radiation to the `relxilllp` reflection spectrum (`refl_frac = -1`), i.e. the change
of the energy flux in the band when switching the returning radiation on:

| h [Rg] | band                                     | bug    | fixed  |
|--------|------------------------------------------|--------|--------|
| 3      | 0.1-1 keV                                | +1.9%  | +1.4%  |
| 3      | 6-7 keV                                  | +3.6%  | +2.7%  |
| 6      | 0.1-1 keV, 1-10 keV, 10-100 keV, 6-7 keV | < 0.1% | < 0.1% |

The returning radiation contribution was overestimated by a factor of about 1.3-1.5. The band-integrated fluxes change
by less than 1% (h = 3) and less than 0.1% (h = 6), but the spectral shape changes more, in particular for low source
heights and high spin (see the reference data comparison below: up to 6% rms per energy bin for h = 1.1 r_H).

## How it was found

While checking the normalization of the reflected spectrum of the `relxillBB` returning radiation (energy
conservation: the angle-averaged reflected flux must be below the returning flux), the angle-averaged slab weighting
gave a reflected/returning energy-flux ratio of 0.65 instead of the correct 0.43. The unit test
`[xillverRTR][norm]` (`test/unit/tests-returnrad-bbody.cpp`) used `get_semi_infinite_slab_flux_weights()` and only
checked the ratio to be < 1, so it passed with the wrong weights. The test
`" angle-averaged xillver flux equals angle-weighted intensity with slab flux weights "` (`[xillverRTR][angledist]`)
re-implemented the same `* 180.0 / M_PI` conversion and therefore confirmed the buggy routine.

## Fix

- `src/Xillspec.cpp:125`: `norm_factor_semi_infinite_slab(xill_spec->incl[ii])` (inclination is already in degrees).
- `test/unit/tests-returnrad-bbody.cpp`, new test `[xillverRTR][slabweights]`: checks
  `get_semi_infinite_slab_flux_weights()` against independently computed weights 0.5 cos(incl_deg) / n_incl and
  their sum of 0.25, and `calc_angle_averaged_xillver_flux()` bin by bin against the sum of the table spectra
  weighted with these weights, for the standard xillver and the xillverRTR table. `[xillverRTR][norm]` uses these
  weights for the energy-conservation check of the `relxillBB` reflection.
- `test/unit/tests-returnrad-bbody.cpp`, `[xillverRTR][angledist]`: removed. It compared
  `calc_angle_averaged_xillver_flux()` with `calc_xillver_angdep()` / `get_angle_weighted_xill_spec()` using weights
  re-implemented with the library formula; as all of them use the same summation, it could only detect a deviation
  of the weights from that formula, not an error of the formula itself, and is superseded by `[slabweights]`.
- `src/common.h`: documented the units of all inclinations: `relParam::incl` [rad], `xillParam::incl` [deg],
  `xillTableParam::incl` [deg], `xillTable::incl` [deg], `xillSpec::incl` [deg].
- `test/unit/test-xilltab.cpp`, `[xilltab][incl-units]`: new test checking these units for the standard xillver
  table (`relxill`) and the xillverRTR table (`relxillBBxill`).

## Verification

- `bin/tests "[xillverRTR]"`: `[slabweights]` passes (it fails with the buggy weights); the angle-averaged
  reflected/returning energy-flux ratio of `relxillBB` is 0.43 per zone sum (was 0.65 with the buggy weights).
- `bin/tests "[incl-units]"`: passes for both tables.
- Full unit test suite: all tests pass after lowering the threshold of `" Changing Rin should result in a change in
  line shape "` (`test/unit/tests-returnrad.cpp:475`) from 1e-3 to 1e-4, see below.
- Reference data (`cd test/refdata && make test`): 62 / 77 pass after the fix, all 77 pass with the code before the
  fix (same commit, built in a separate worktree). The 15 failing comparisons are all lamp-post and ring-source
  models with returning radiation:

| reference spectrum                             | h       | a     | deviation | rms change per bin |
|------------------------------------------------|---------|-------|-----------|--------------------|
| relxilllp_special_refdat_0009 (refl_frac = -1) | 1.1 r_H | 0.998 | 1.38e-3   | 6.2%               |
| relxilllp_special_refdat_0007                  | 1.1 r_H | 0.998 | 1.28e-3   | 5.8%               |
| relxilllpCp_special_refdat_0004                | 2       | 0.998 | 1.16e-3   | 5.2%               |
| relxilllp_special_refdat_0001                  | 2       | 0.998 | 4.16e-4   | 1.9%               |
| relxilllp_special_refdat_0003 (refl_frac = -1) | 2       | 0.998 | 3.74e-4   | 1.7%               |
| relxilllpCp_special_refdat_0006                | 3       | 0.8   | 1.58e-4   | 0.7%               |
| relxill_ring_defparam_refdat_0001              | (ring)  | 0.998 | 1.42e-4   | 0.6%               |
| relxilllpCp_defparam_refdat_0001               | 6       | 0.998 | 1.40e-4   | 0.6%               |
| relxilllp_defparam_refdat_0001                 | 6       | 0.998 | 1.30e-4   | 0.6%               |
| relxill_ring_random_refdat_0002                | (ring)  | 0.070 | 9.02e-5   | 0.4%               |
| relxill_ring_random_refdat_0003                | (ring)  | 0.346 | 4.91e-5   | 0.2%               |
| relxilllp_random_refdat_0004                   | 17.6    | 0.724 | 3.26e-5   | 0.15%              |
| relxilllp_random_refdat_0003                   | 53.5    | 0.915 | 1.22e-5   | 0.05%              |
| relxilllp_special_refdat_0006 (refl_frac = -1) | 30      | 0.998 | 1.15e-5   | 0.05%              |
| relxilllp_special_refdat_0004                  | 30      | 0.998 | 1.12e-5   | 0.05%              |

The deviation is the test metric sqrt(sum((model/ref - 1)^2)) / N with N = 2000 bins and a limit of 1e-5; the rms
change per bin is deviation * sqrt(N). The reference data need to be regenerated for these models.

### Side effect on `tests-returnrad.cpp:475`

The test (introduced in `cce7908` for the interpolation of the returning radiation table when Rin changes) evaluates
`relxilllp` (a = 0.995, returning radiation on) for Rin = r_isco(a=0.990) and r_isco(a=0.99025) and requires the
total photon flux (sum over the 0.1-1000 keV grid with 4000 bins) to change by more than 1e-3 (absolute).

|       | returning rad. | sum1     | sum2     | abs(sum1 - sum2) | max. relative change per bin |
|-------|----------------|----------|----------|------------------|------------------------------|
| bug   | off            | 128.4858 | 128.4799 | 5.9e-3           | 8.5e-5                       |
| bug   | on             | 131.6953 | 131.6976 | 2.3e-3           | 1.3e-4                       |
| fixed | on             | 130.6400 | 130.6396 | 3.5e-4           | 1.1e-4                       |

The small change in Rin changes the reflection of the primary source and of the returning radiation in opposite
directions. With the corrected (weaker) returning radiation both nearly cancel in the total flux, while the spectrum
still changes per bin by about 1e-4 (as before the fix). The test criterion (absolute change of the total flux) does
not measure the change of the line shape and fails by this coincidence; the model itself still reacts to Rin.
