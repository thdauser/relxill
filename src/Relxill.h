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

    Copyright 2021 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELXILL_H_
#define RELXILL_H_

#include "XspecSpectrum.h"

extern "C" {
#include "relmodels.h"
#include "relbase.h"
#include "writeOutfiles.h"
#include "relutility.h"
#include "xillspec.h"
}


enum cached {
  yes,
  no
};

class CachingStatus {

 public:
  CachingStatus() = default;

  [[nodiscard]] int is_all_cached() const {
    if (energy_grid == cached::yes && relat == cached::yes && xill == cached::yes) {
      return 1;
    } else {
      return 0;
    }
  }

  [[nodiscard]] int recomput_relat() const {
    if (relat == cached::no) {
      return 1;
    } else {
      return 0;
    }
  }

  cached energy_grid = cached::no;
  cached relat = cached::no;
  cached xill = cached::no;
};

void relxill_kernel(const XspecSpectrum &spectrum,
                    xillParam *xill_param,
                    relParam *rel_param,
                    int *status);

#endif
