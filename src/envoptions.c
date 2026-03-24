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

int is_debug_run(void) {
  return (int) get_env_otherwise_default("DEBUG_RELXILL", 0);
}

int shouldAuxInfoGetPrinted(void) {
  return (int) get_env_otherwise_default("RELXILL_PRINT_DETAILS", 0);
}

int shouldOutfilesBeWritten(void) {
  return (int) get_env_otherwise_default("RELXILL_WRITE_FILES", 0);
}


int get_from_env_rrad_bbody_order(void) {
  int def = 1;
  int val = (int) get_env_otherwise_default("RELXILL_RRAD_BBODY_ORDER", def);
  if (val < 0 || val > 3) {
    relxill_warning("only values between 0 and 3 are allowed for RELXILL_RRAD_BBODY_ORDER");
    return def;
  }
  return val;
}

int get_from_env_rrad_bbody_tprofile(void) {
  int def = 1;
  int val = (int) get_env_otherwise_default("RELXILL_RRAD_BBODY_TPROFILE", def);
  if (val < 1 || val > 2) {
    relxill_warning("only values between 1 and 2 are allowed for RELXILL_RRAD_BBODY_TPROFILE");
    return def;
  }
  return val;
}

int constantDiskDensity(void) {
  return (int) get_env_otherwise_default("RELXILL_CONSTANT_DENSITY", 0);
}

int do_not_normalize_relline(void) {
  return (int) get_env_otherwise_default("RELLINE_PHYSICAL_NORM", 0);
}

int get_from_env_relline_num_zones(void) {
  return (int) get_env_otherwise_default("RELXILL_NUM_RZONES", 0);
}

void print_env_options(void) {
  printf("\n === RELXILL Environment Variables ===\n\n");

  char *path = getenv("RELXILL_TABLE_PATH");
  printf("  RELXILL_TABLE_PATH          = %s\n", path != NULL ? path : RELXILL_TABLE_PATH);
  printf("  DEBUG_RELXILL               = %d\n", is_debug_run());
  printf("  RELXILL_PRINT_DETAILS       = %d\n", shouldAuxInfoGetPrinted());
  printf("  RELXILL_WRITE_FILES         = %d\n", shouldOutfilesBeWritten());
  printf("  RELXILL_RRAD_BBODY_ORDER    = %d\n", get_from_env_rrad_bbody_order());
  printf("  RELXILL_RRAD_BBODY_TPROFILE = %d\n", get_from_env_rrad_bbody_tprofile());
  printf("  RELXILL_CONSTANT_DENSITY    = %d\n", constantDiskDensity());
  printf("  RELLINE_PHYSICAL_NORM       = %d\n", do_not_normalize_relline());
  printf("  RELXILL_NUM_RZONES          = %d\n", get_from_env_relline_num_zones());

  printf("\n");
}
