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

    Copyright 2022 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELXILL_H_
#define RELXILL_H_

#include "XspecSpectrum.h"
#include "Relbase.h"
#include "Xillspec.h"
#include "IonGradient.h"
#include "ModelParams.h"

extern "C" {
#include "writeOutfiles.h"
#include "relutility.h"
}

#define SPIN_MIN_RRAD_CALC_CORRFAC (0.0)  // minimal value for which returning radiation is calculated (for relxill_kernel)

enum cached {
  yes,
  no
};

class SpectrumZones{

 public:
  SpectrumZones(const double* _energy, int _num_flux, int _num_zones):
      num_flux_bins(_num_flux), num_zones(_num_zones), num_ener_bins(_num_flux+1) {

    flux = new double*[num_zones];
    for (int ii=0; ii<num_zones; ii++){
      flux[ii] = new double[num_flux_bins];
      for (int jj=0; jj<num_flux_bins; jj++) {
        flux[ii][jj] = 0.0;
      }
    }

    m_energy = new double[num_flux_bins+1];
    for (int ii=0; ii<num_ener_bins; ii++) {
      m_energy[ii] = _energy[ii];
    }
  }

  ~SpectrumZones(){
    for (int ii=0; ii<num_zones; ii++) {
      delete[] flux[ii];
    }
    delete[] m_energy;
    delete[] flux;
  }

  [[nodiscard]] double* energy() const{
    return m_energy;
  }

  const int num_flux_bins;
  const int num_zones;
  const int num_ener_bins;
  double** flux;

 private:
  double* m_energy;

};



class CachingStatus {

 public:
  CachingStatus() = default;

  CachingStatus(relParam *rel_param, xillParam *xill_param, specCache *spec_cache, const XspecSpectrum &spectrum) {
    check_caching_parameters(rel_param, xill_param);
    check_caching_energy_grid(spec_cache, spectrum);
  }

  [[nodiscard]] auto any_parameter_changed() const -> int {
    if (m_cache_relat == cached::no || m_cache_xill == cached::no) {
      return 1;
    }
    return 0;
  }


  [[nodiscard]] auto is_all_cached() const -> bool {
    return m_cache_energy_grid == cached::yes && m_cache_relat == cached::yes && m_cache_xill == cached::yes;
  }

  [[nodiscard]] auto only_energy_grid_changed() const -> bool {
    return m_cache_energy_grid == cached::no && m_cache_relat == cached::yes && m_cache_xill == cached::yes;
  }


  [[nodiscard]] auto recomput_relat() const -> int {
    if (m_cache_relat == cached::no) {
      return 1;
    }
    return 0;
  }


  void check_caching_parameters(const relParam *rel_param, const xillParam *xill_param);

  void check_caching_energy_grid(specCache *cache, const XspecSpectrum &spectrum);

  [[nodiscard]] auto relat() const -> cached { return m_cache_relat; }

  [[nodiscard]] auto xill() const -> cached { return m_cache_xill; }

  [[nodiscard]] auto energy_grid() const -> cached { return m_cache_energy_grid; }


 private:
  static auto did_number_of_zones_change(const relParam *rel_param) -> int;

  cached m_cache_energy_grid = cached::no;
  cached m_cache_relat = cached::no;
  cached m_cache_xill = cached::no;
};

void relxill_kernel(const XspecSpectrum &spectrum,
                    const ModelParams &params,
                    int *status);

rradCorrFactors* calc_rrad_corr_factors(xillSpec **xill_spec, const RadialGrid &rgrid,
                                        xillTableParam *const *xill_table_param, int *status);

#endif
