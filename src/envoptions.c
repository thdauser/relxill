/*
   This file is part of the RELXILL model code.

   RELXILL is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   RELXILL is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.
   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

    Copyright 2025 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "envoptions.h"
#include "config.h"

double get_env_otherwise_default(const char *env, double def_value) {
  char *value_string = getenv(env);
  if (value_string != NULL) {
    return strtod(value_string, NULL);
  } else {
    return def_value;
  }
}

/** get the relxill table path (dynamically from env variable)  **/
char *get_relxill_table_path(void) {
  char *path;
  path = getenv("RELXILL_TABLE_PATH");
  if (path != NULL) {
    return path;
  } else {
    return RELXILL_TABLE_PATH;
  }
}

/** check if we are currently debugging the model **/
int is_debug_run(void) {
  char *env;
  env = getenv("DEBUG_RELXILL");
  if (env != NULL) {
    int debug = (int) atof(env);
    if (debug == 1) {
      return 1;
    }
  }
  return 0;
}


/** check if we are currently debugging the model **/
int shouldAuxInfoGetPrinted(void) {
  char *env = getenv("RELXILL_PRINT_DETAILS");
  if (env != NULL) {
    int envval = (int) strtod(env, NULL);
    if (envval == 1) {
      return 1;
    }
  }
  return 0;
}

/** check if we are currently debugging the model **/
int shouldOutfilesBeWritten(void) {
  char *env;
  env = getenv("RELXILL_WRITE_FILES");
  if (env != NULL) {
    int envval = (int) strtod(env, NULL);
    if (envval == 1) {
      return 1;
    }
  }
  return 0;
}


int get_from_env_rrad_bbody_order(void) {
  char *env = getenv("RELXILL_RRAD_BBODY_ORDER");
  if (env != NULL) {
    int envval = (int) strtod(env, NULL);
    if (envval >= 0 && envval <= 3) {
      return envval;
    }
    relxill_warning("only values between 0 and 3 are allowed for RELXILL_RRAD_BBODY_ORDER");
  }
  return 1;
}

int get_from_env_rrad_bbody_tprofile(void) {
  char *env = getenv("RELXILL_RRAD_BBODY_TPROFILE");
  if (env != NULL) {
    int envval = (int) strtod(env, NULL);
    if (envval >= 1 && envval <= 2) {
      return envval;
    }
    relxill_warning("only values between 1 and 2 are allowed for RELXILL_RRAD_BBODY_TPROFILE");
  }
  return 1;
}

/** check if we are currently debugging the model **/
int constantDiskDensity(void) {
  char *env;
  env = getenv("RELXILL_CONSTANT_DENSITY");
  if (env != NULL) {
    int envval = (int) strtod(env, NULL);
    if (envval == 1) {
      return 1;
    }
  }
  return 0;
}


/** check if we should return the relline/relconv physical norm from ENV **/
int do_not_normalize_relline(void) {
  char *env;
  env = getenv("RELLINE_PHYSICAL_NORM");
  if (env != NULL) {
    int phys_norm = (int) strtod(env, NULL);
    if (phys_norm == 1) {
      return 1;
    }
  }
  return 0;
}

int get_from_env_relline_num_zones(void) {
  char *env;
  env = getenv("RELXILL_NUM_RZONES");
  if (env != NULL) {
    return (int) atof(env);
  }
  return 0;
}