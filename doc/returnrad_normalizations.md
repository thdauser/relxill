# Return-Radiation Flux Normalizations {#returnrad_normalizations}

This page documents the normalization factors applied inside
`calc_zone_direct_plus_xillver_reflected_flux` in
`src/Relreturn_BlackBody.cpp`.

## Context

For each radial zone the model combines two contributions:

1. The **direct** thermal blackbody emission from the disk surface
   (`returnSpec->specBB[izone]`).
2. The **reflected** xillver component irradiated by the returning
   thermal radiation (`returnSpec->specRet[izone]`).

Both must be folded with the relativistic line profile of the zone
(`rel_profile`) and combined on a common energy grid before being added
into the zone flux array `xill_flux_returnrad`.

## Normalizations applied

\f[
F_\mathrm{zone}(E) \;=\;
N_\mathrm{dir}\, F_\mathrm{BB,direct}(E) \;+\;
N_\mathrm{refl}\, F_\mathrm{xillver,returned}(E)
\f]

TODO — fill in once the derivation is finalized:

- **`N_dir`** — accounts for ...
- **`N_refl`** — accounts for ...
- Energy-grid resampling factor (per-bin width vs. flux density) ...
- Solid-angle / area factor between zone surface and observer frame ...

## Why the factors are split this way

TODO — record the *reason* for the chosen split (e.g. caching of the
unnormalized xillver spectrum, decoupling ionization from geometry,
matching the convention used elsewhere in `Relreturn_Corona.cpp`).

## Pitfalls

- Forgetting either factor leads to a flux that scales correctly with
  energy but is off by a constant — easy to miss in log-log plots.
- The reflected component is already convolved with the disk emissivity
  inside the RTR table, so do *not* re-apply the emissivity here.
