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

    Copyright 2020 Thomas Dauser, Remeis Observatory & ECAP
*/

#ifndef RELXILL__CPPSPECTRUM_H_
#define RELXILL__CPPSPECTRUM_H_

#include <xsTypes.h>

typedef RealArray Array; // using the Xspec defined std::valarray type




class XspecSpectrum {

 public:
  XspecSpectrum(const double *_energy, double *_flux, size_t _nbins_xspec)
      : m_flux{_flux}, m_num_flux_bins{_nbins_xspec} {
    // need to allocate the energy grid, as we are not allowed to change the original energy grid
    // (Xspec requires it to be constant, but we need to shift it in energy)
    m_ener = new double[n_energy()];
    for (size_t ii = 0; ii < n_energy(); ii++) {
      m_ener[ii] = _energy[ii];
    }
  };

  ~XspecSpectrum() {
    delete (m_ener);
  }

  [[nodiscard]] double *energy() const {
    return m_ener;
  }

  [[nodiscard]] double *flux() const {
    return m_flux;
  }

  [[nodiscard]] size_t n_energy() const {   // array holds num_bins+1 bins, as bin_lo and bin_hi are combined
    return m_num_flux_bins + 1;
  }

  [[nodiscard]] int num_flux_bins() const {   // array holds nener+1 bins, as bin_lo and bin_hi are combined
    return m_num_flux_bins;
  }

  /**
   * shift the spectrum in redshift and for a line at line for 1 keV
   */
  void shift_energy_grid_1keV(double line_energy, double z) const {
    for (size_t ii = 0; ii < m_num_flux_bins + 1; ii++) {
      m_ener[ii] *= (1 + z) / line_energy;
    }
  }

  /**
   * shift the spectrum in redshift
   */
  void shift_energy_grid(double line_energy, double z) const {
    shift_energy_grid_1keV(1.0, z);
  }

 private:
  double *m_ener{nullptr};
  double *m_flux{nullptr};
  size_t m_num_flux_bins;

};

#endif //RELXILL__CPPSPECTRUM_H_
