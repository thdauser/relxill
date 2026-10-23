# General
This readme describes the manual tests and plots for returning thermal radiation (RTR). Once this is implemented and working for the main, the tests should be moved to "tests/e2e/" and the folder changed to be simply for plotting

# Tests
- the test scripts are in the subfolder tests/ and the functions in there are called by `master_test_driver.sl`
- single tests can be run with `ìsis-scripts tests/example.sl 1`
- reference data is given in `refdata/`

# Plots
- The main **diagnose plots** are called in `plot_zones_diagnose_spectra.sl`
- Further smaller plot scripts are ther other `plot_*.sl` scripts, which rely on calling the model; some may be broken by changes to the output files in relxillBB
- all plots are stored in `plots/`

# Output files

- the model writes the spectra of the zones to `__relxillOutput-rrad-bbody-*.fits` (with `RELXILL_WRITE_FILES=1`, see
  `create_output_files.sl`); all of them are fluxes (rest frame: angle-averaged). The table of all files, their content
  and flux conversion is in `doc/returnrad_disk_emission.md` ("Diagnostic output files")
- the `*-total` and `*-primary` files in this folder are from before 2026-10-21 and are no longer written




