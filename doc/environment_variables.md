# RELXILL Environment Variables

Environment variables that control runtime behavior of the RELXILL model.
All variables are optional and fall back to their default values if unset.

## Table Path

| Variable | Default | Values | Description |
|----------|---------|--------|-------------|
| `RELXILL_TABLE_PATH` | set by the script to compile-time path | directory path | Path to the reflection table FITS files |

## Debug and Output

| Variable | Default | Values | Description |
|----------|---------|--------|-------------|
| `DEBUG_RELXILL` | 0 | 0, 1 | Enable debug mode |
| `RELXILL_PRINT_DETAILS` | 0 | 0, 1 | Print auxiliary model info |
| `RELXILL_WRITE_FILES` | 0 | 0, 1 | Write (many) additional output files |

## Model Configuration

| Variable | Default | Values | Description |
|----------|---------|--------|-------------|
| `RELXILL_NUM_RZONES` | not set | 1–50 | Override the number of radial zones (ignored if 0 or unset) |
| `RELXILL_CONSTANT_DENSITY` | 0 | 0, 1 | Use constant disk density |
| `RELLINE_PHYSICAL_NORM` | 0 | 0, 1 | Return the physical normalization for relline/relconv |

## Blackbody Return Radiation

| Variable | Default | Values | Description |
|----------|---------|--------|-------------|
| `RELXILL_RRAD_BBODY_ORDER` | 1 | 0–3 | Order of the blackbody return radiation calculation |
| `RELXILL_RRAD_BBODY_TPROFILE` | 1 | 1–2 | Temperature profile for blackbody return radiation |
| `RELXILL_BBRET_NOREFL` | 0 | 0, 1 | Skip xillver reflection calculation for blackbody return radiation |
