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

#include <stdarg.h>

#define MAX_REPORTED_ENVS 32
static const char *reported_envs[MAX_REPORTED_ENVS];
static char *reported_values[MAX_REPORTED_ENVS];
static int num_reported_envs = 0;
static int suppress_env_reporting = 0;

/* Returns 1 if the stored value for env differs from value_str (or env is
 * newly tracked), and updates the stored value in either case. Returns 0 if
 * unchanged, or if the tracking capacity is exhausted. */
static int env_value_changed(const char *env, const char *value_str) {
  for (int i = 0; i < num_reported_envs; i++) {
    if (strcmp(reported_envs[i], env) == 0) {
      if (strcmp(reported_values[i], value_str) == 0) return 0;
      free(reported_values[i]);
      reported_values[i] = strdup(value_str);
      return 1;
    }
  }
  if (num_reported_envs >= MAX_REPORTED_ENVS) return 0;
  reported_envs[num_reported_envs] = env;
  reported_values[num_reported_envs] = strdup(value_str);
  num_reported_envs++;
  return 1;
}

static void report_env_if_changed(const char *env, const char *value_str, const char *fmt, ...) {
  if (suppress_env_reporting) return;
  if (!env_value_changed(env, value_str)) return;

  va_list ap;
  va_start(ap, fmt);
  printf(" *** RELXILL: env %s ", env);
  vprintf(fmt, ap);
  printf("\n");
  va_end(ap);
}

double get_env_otherwise_default(const char *env, double def_value) {
  char *value_string = getenv(env);
  if (value_string != NULL) {
    double val = strtod(value_string, NULL);
    if (val != def_value) {
      report_env_if_changed(env, value_string, "= %g (default %g)", val, def_value);
    }
    return val;
  } else {
    return def_value;
  }
}

/** get the relxill table path (dynamically from env variable)  **/
char *get_relxill_table_path(void) {
  char *path = getenv("RELXILL_TABLE_PATH");
  if (path != NULL) {
    return path;
  }
  return RELXILL_TABLE_PATH;
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

int use_iongrad_relxillbb(void) {
  return (int) get_env_otherwise_default("RELXILL_RRAD_BBODY_IONGRAD", 0);
}

int mirrored_reflection_instead_of_xillver(void) {
  return (int) get_env_otherwise_default("RELXILL_RRAD_BBODY_MIRRORED", 0);
}


void print_env_options(void) {
  suppress_env_reporting = 1;
  printf("\n === RELXILL Environment Variables ===\n\n");

  char *path = getenv("RELXILL_TABLE_PATH");
  printf("  RELXILL_TABLE_PATH          = %s\n", path != NULL ? path : RELXILL_TABLE_PATH);
  printf("  DEBUG_RELXILL               = %d\n", is_debug_run());
  printf("  RELXILL_PRINT_DETAILS       = %d\n", shouldAuxInfoGetPrinted());
  printf("  RELXILL_WRITE_FILES         = %d\n", shouldOutfilesBeWritten());
  printf("  RELXILL_RRAD_BBODY_ORDER    = %d\n", get_from_env_rrad_bbody_order());
  printf("  RELXILL_RRAD_BBODY_TPROFILE = %d\n", get_from_env_rrad_bbody_tprofile());
  printf("  RELXILL_RRAD_BBODY_IONGRAD  = %d\n", use_iongrad_relxillbb());
  printf("  RELXILL_RRAD_BBODY_MIRRORED = %d\n", mirrored_reflection_instead_of_xillver());
  printf("  RELXILL_CONSTANT_DENSITY    = %d\n", constantDiskDensity());
  printf("  RELLINE_PHYSICAL_NORM       = %d\n", do_not_normalize_relline());
  printf("  RELXILL_NUM_RZONES          = %d\n", get_from_env_relline_num_zones());

  printf("\n");
  suppress_env_reporting = 0;
}
