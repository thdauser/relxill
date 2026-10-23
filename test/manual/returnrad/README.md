# General
This folder holds the diagnostic plots for the returning thermal radiation (RTR, `relxillBB`). They need ISIS and the
model compiled with `make model` (loaded from `build/`).

# Tests
The tests are unit tests (`bin/tests "[rrad-bbody]"`). The former manual tests of this folder were moved there on
2026-10-23:
- `testEmissivityBbody.sl` (energy flux of the primary of each zone vs. SS73) → `"disk emission (1b)"` in
  `test/unit/tests-returnrad-bbody-diskemission.cpp`
- `testRframeReturnradRefdata.sl` (returning spectra of the zones vs. reference data) → `"returning radiation:
  returning spectra of the zones agree with the reference data"` in `test/unit/tests-returnrad-bbody-spectra.cpp`,
  with the reference data in `test/unit/refdata/refvalues_rrad_bbody_specRet_a*.fits` (written with
  `bin/tests "[.create-refdata]"`)

# Plots

- `make plot`: the main **diagnose plots** of `relxillBB` (`plot_zones_diagnose_spectra.sl`, first order of the
  returning radiation, reflection only):
  - `plots/relxillBB_perzone_rframe.pdf`: rest frame, reflected (solid), returning (dashed), and direct (dotted)
    spectrum of the zones
  - `plots/relxillBB_perzone_observed.pdf`: observed, reflected (solid), mirrored (dashed), and direct (dotted) spectrum
    of the zones, their
    sums over the zones, and direct + reflected (the model)
- `make plot` also creates `plots/relxillBB_NS_perzone_mirror_orders.pdf` (`plot_perzone_mirror_orders.sl`): observed
  spectra of the zones of `relxillBB_NS` with the mirrored returning radiation, 1st order (solid), 2nd order (dashed),
  and direct (dotted); the 2nd order alone is the difference of the runs with `RELXILL_RRAD_BBODY_ORDER=2` and 1, as
  the model returns the sum of all orders
- further plot scripts: `plot_higher_order_rrad.sl` (higher orders of the mirrored returning radiation) and
  `plot_iongradient_return.sl` (ionization gradient), run with `isis-script <script>`
- all plots are stored in `plots/`
- the rest-frame plots show the spectra of the zones normalized per zone by N_xill of the returning spectrum of the zone,
  as for the xillver table in the model (all zones equally visible, relative normalization of primary, returning, and
  reflected spectrum kept); the importance of the zones is shown by the observed plots (see "Diagnostic output files"
  in `doc/returnrad_disk_emission.md`)

# Output files
- the model writes the spectra of the zones to `__relxillOutput-rrad-bbody-*.fits` (with `RELXILL_WRITE_FILES=1`, see
  `create_output_files.sl`); all of them are fluxes (rest frame: angle-averaged). The table of all files, their content
  and flux conversion is in `doc/returnrad_disk_emission.md` ("Diagnostic output files")
- `make clean` removes them
