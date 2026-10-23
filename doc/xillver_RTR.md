# The xillverRTR Table for relxillBB {#xillver_RTR}

This page documents the xillver table of the returning thermal radiation (RTR) used by `relxillBB`: how the model uses
it, what the table has to provide, how its incident spectra are produced, and what needs to change in the model for the
final table. The current table `xillverRTR_v0.3.fits` is a **mock version** for a single parameter set; the final table
will in addition depend on \f$kT_\mathrm{bb}\f$, logxi, logN, and the spin (later also A_Fe). The normalization of the
reflection is
described in \ref returnrad_normalizations, the angular dependence in \ref xillver_angular_dependence, and the zones of
the disk in \ref returnrad_disk_emission (state of 2026-10-23).

## How relxillBB uses the table

The reflection is calculated separately for each of the radial zones of the returning radiation (\ref relxillbb_kernel):

1. **Zone and table spectra.** For zone `izone`, the model sets `radzone = izone`, which is the table parameter
   `radindex` (`get_xilltab_param()`). The table is interpolated in all its parameters (generic interpolation in
   `Xilltable_handling.cpp`), and the spectra of all table inclinations are read.
2. **Normalization to the table convention.** As for all xillver tables, the spectra are divided by
   \f$10^{\mathrm{logxi}}\f$ and \f$10^{\mathrm{logN}-15}\f$ (`renorm_xill_spec()`), i.e., they are normalized to a
   density of \f$10^{15}\,\mathrm{cm}^{-3}\f$ and \f$\xi = 1\f$.
3. **Angle weighting.** The intensities of the table inclinations are weighted with the GR emission-angle distribution
   of the zone, giving the angle-mixed intensity that is convolved (\ref calc_zone_xillver_reflected_intensity, see
   \ref xillver_angular_dependence).
4. **Absolute normalization.** The spectrum is scaled with \f$N_\mathrm{xill}\f$ of the returning spectrum `specRet` of
   the zone, i.e., its energy flux in 0.1-1000 keV relative to the reference illumination of the table (\ref
   returnrad_normalizations). Only the *shape* of the incident spectrum of the table enters; its absolute level
   is set by the model.
5. **Convolution** with the relativistic line profile of the zone.

The incident spectrum of each zone of the table must therefore be `specRet` of the same zone (for the same parameters),
so that the reflection belongs to the illumination it is scaled with. A mismatch of the shape must be removed by
tabulating the correct incident spectrum, never by a correction factor (see `xillver_normalization_summary.md`,
§4.2, in `test/manual/returnrad/`).

## The current table (mock version v0.3)

| property           | value                                                                                                                                                                                                                                  |
|--------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| file               | `xillverRTR_v0.3.fits` (`XILLTABLE_RRAD_FILENAME` in `src/common.h`), in `RELXILL_TABLE_PATH`                                                                                                                                          |
| format             | Xspec table model (`PARAMETERS`, `ENERGIES`, `SPECTRA`), model name `xillverRTR`, photons/cm²/s per bin                                                                                                                                |
| parameters         | `radindex` (50 values, 0-49), `Incl` (10 values, the bin centres of a grid uniform in \f$\cos i\f$, 18-87 deg)                                                                                                                         |
| energy grid        | 2999 bins, 0.072-1000 keV                                                                                                                                                                                                              |
| incident spectra   | `specRet` of the 50 zones for a = 0.998, \f$R_\mathrm{in}\f$ = ISCO, \f$R_\mathrm{out}\f$ = 1000 Rg, \f$kT_\mathrm{bb}\f$ = 2 keV (default temperature profile and order), calculated in 2026-03 (`relxill-table`, version `test_0.3`) |
| xillver parameters | logxi = 2, logN = 18                                                                                                                                                                                                                   |

The incident spectra were calculated with the model before the changes of 2026-10-21 (primary of the zones from the
radial points of the line profiles, see `test/manual/returnrad/changes_2026-10-21_disk_emission.md`); the returning
radiation of the zones calculated with the current model differs from them by up to about 0.2% in the energy flux.

As the table only covers this parameter set, the model fixes a, \f$R_\mathrm{in}\f$, \f$R_\mathrm{out}\f$, logxi, and
logN (`check_parameter_bounds_relxillbb()`), and \f$kT_\mathrm{bb}\f$ by the parameter limits in
`src/modelfiles/lmodel_relxill_devel.dat`. Only the inclination, z, boost, and the normalization are free.

## The zones

The zones are the radial zones of the returning-radiation table `table_returnRad_v20220301.fits`
(`get_rrad_fractions()` in `src/Relreturn_Table.cpp`):

- The table contains 43 spins from 0 to 0.9982 (no negative spins), each with 50 zones from the ISCO of this spin to
  1000 Rg (`RETURNRAD_TABLE_NR`). The radii of the zones therefore depend on the spin (e.g., zone 0 starts at
  6.0 Rg for a = 0 and at 1.24 Rg for a = 0.998).
- For a given spin, the model uses the zones of the **next tabulated spin** \f$\geq\f$ a (no interpolation in spin,
  `select_spinIndexForTable()`), trimmed to \f$[R_\mathrm{in}, R_\mathrm{out}]\f$ (`allocate_radial_grid()`): with
  \f$R_\mathrm{in}\f$ = ISCO and \f$R_\mathrm{out}\f$ = 1000 Rg, all 50 zones are used; a larger \f$R_\mathrm{in}\f$
  removes the inner zones.
- `radindex` of the table is the index of the zone in this grid. The number of zones is not affected by
  `RELXILL_NUM_RZONES` (which only sets the zones of the lamp-post and ionization-gradient models).

### The spin defines the zones

The spin sets the inner edge of the disk (\f$R_\mathrm{in}\f$ = ISCO) and therefore the radial zones: for each spin,
the 50 zones span a different radial range, and a given `radindex` corresponds to a different radius (and a different
incident spectrum) for each spin. In the final table with a spin dimension, `radindex` therefore refers to the zone
grid of each spin, and the spins of the table should be spins of the returning-radiation table.

Interpolating the table in spin is possible, but not straightforward: at a fixed `radindex`, the spectra of two
tabulated spins belong to zones at different radii, so a plain interpolation in spin mixes the reflection of different
radii, while the zones of the model (and their line profiles) are those of a single tabulated spin. A consistent
interpolation would need to account for the change of the zones with the spin (e.g., interpolating at the same radius
instead of the same `radindex`). Using only the tabulated spins avoids this, as for the returning radiation (next
tabulated spin \f$\geq\f$ a).

## Incident spectra: the relxill-table pipeline

The incident spectra are calculated with the model itself, in the (private) repository `relxill-table`.

**Dependency:** the pipeline relies on the name, the format, and the units of `rframe-bbodyRet.fits` (and on
`specRet` being the flux incident on the zone). Any change there must be propagated to `relxill-table`. The
returning spectra of the zones are protected by the unit test *"returning radiation: returning spectra of the zones
agree with the reference data"* (`test/unit/tests-returnrad-bbody-spectra.cpp`); its reference file
`test/unit/refdata/refvalues_rrad_bbody_specRet_a0.998.fits` contains them for the parameters of the mock table
(\f$kT_\mathrm{bb}\f$ = 2 keV, calculated with the current model, on a coarser energy grid).

## The final table and required changes of the model

The final table will depend on `radindex`, `Incl`, \f$kT_\mathrm{bb}\f$ (`kTbb`), logxi (`logXi`), logN (`Dens`), and
the spin, with the incident spectra of the zones calculated for each combination of spin and \f$kT_\mathrm{bb}\f$. A
future version will also include A_Fe (`A_Fe`). The following needs to change in the model:

- **Spin as a table parameter.** The table reader only knows the parameter names in `global_param_names`
  (`xilltable_reading.c`; `NAME_*` and `N_PARAM_MAX` in `src/common.h`): `Gamma`, `A_Fe`, `logXi`, `Ecut`, `kTe`,
  `Dens`, `kTbb`, `A_CO`, `Frac`, `Incl`, `radindex`. A spin parameter does not exist yet and needs to be added and
  set from the spin of the model.
- **Number of table parameters.** The reader supports at most 6 parameters, including `Incl`
  (`convertTo6dimTableIndex()`). The final table with `radindex`, `Incl`, `kTbb`, `logXi`, `Dens`, and the spin is
  exactly at this limit; with `A_Fe` (7 parameters), the reader and its interpolation need to be extended.
- **Free parameters.** Remove the fixed values in `check_parameter_bounds_relxillbb()` and open the parameter limits of
  `relxillBB` in `lmodel_relxill_devel.dat` to the range of the table.
- **Check of the zone grid.** An out-of-range `radindex` is currently clamped to the last tabulated zone without an
  error (`resetInpvalsToBoundaries()`, see `xillver_normalization_summary.md`, §4.1). The model should check that the
  zone grid of the table matches the zones of the returning radiation for the spin of the model.
- **Spin without interpolation of the zones.** The zones of the returning radiation use the next tabulated spin
  \f$\geq\f$ a, whereas the xillver table interpolates in all its parameters. The spin of the xillver table must not be
  interpolated at a fixed `radindex` (see "The spin defines the zones"): either the table is read at the same tabulated
  spin as the returning radiation, or the interpolation accounts for the change of the zones.

## Open questions for the table design

Each xillver spectrum takes days to calculate, so the number of spectra (zones × inclinations × the grid of
\f$kT_\mathrm{bb}\f$, logxi, logN, spin, and later A_Fe) is decisive. Open questions, to be worked out separately:

- How many zones are needed, and how should they be distributed (e.g., how similar are the incident and the reflected
  spectra of neighbouring zones, and how much does each zone contribute to the observed reflection)?
- Can the model interpolate in `radindex` between tabulated zones (the table interpolation is generic), and how large
  is the resulting error of the full model?
- Which spins does the table need (all 43 of the returning-radiation table, or a subset, and how is a spin between
  them treated without interpolating at a fixed `radindex`)?

## Unit tests

In `test/unit/tests-returnrad-bbody-xillver.cpp` (`bin/tests "[xillverRTR]"`):

- *"xillverRTR table: loading spectra without a radial zone (radzone=-1) fails"*
- *"xillverRTR table: number of radial zones equals the radial grid of the returning radiation"*
- *"xillverRTR table: spectra of a radial zone are loaded and non-zero"*
- the angle weighting and normalization tests `[slabweights]`, `[cosnedist]`, `[norm]`, `[angleweight]`, see
  \ref xillver_angular_dependence

and the regression test of the incident spectra (`specRet`) in `test/unit/tests-returnrad-bbody-spectra.cpp`
(`bin/tests "[refdata]"`).
