# General
This readme describes the manual tests and plots for returning thermal radiation (RTR). Once this is implemented and working for the main, the tests should be moved to "tests/e2e/" and the folder changed to be simply for plotting

# Tests
- the test scripts are in the subfolder tests/ and the functions in there are called by `master_test_driver.sl`
- reference data is given in `refdata/`

# Plots
- The main **diagnose plots** are called in `plot_zones_diagnose_spectra.sl`
- Further smaller plot scripts are ther other `plot_*.sl` scripts, which rely on calling the model; some may be broken by changes to the output files in relxillBB
- all plots are stored in `plots/`




