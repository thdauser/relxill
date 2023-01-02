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

#ifndef RELXILL__CPPSPECTRUM_H_
#define RELXILL__CPPSPECTRUM_H_

#include <xsTypes.h>

#include "common.h"

typedef RealArray Array; // using the Xspec defined std::valarray type




class XspecSpectrum {

 public:

  XspecSpectrum(const double *_energy, double *_flux, size_t _nbins_xspec)
      : flux{_flux}, m_num_flux_bins{_nbins_xspec} {
    // need to allocate the energy grid, as we are not allowed to change the original energy grid
    // (Xspec requires it to be constant, but we need to shift it in energy)
    energy = new double[n_energy()];
    for (size_t ii = 0; ii < n_energy(); ii++) {
      energy[ii] = _energy[ii];
    }
  };

  ~XspecSpectrum() {
    delete[] energy;
  }

  [[nodiscard]] size_t n_energy() const {   // array holds num_bins+1 bins, as bin_lo and bin_hi are combined
    return m_num_flux_bins + 1;
  }

  [[nodiscard]] int num_flux_bins() const {   // array holds nener bins, as bin_lo and bin_hi are combined
    return static_cast<int>(m_num_flux_bins);
  }

  /**
   * shift the spectrum such that we can calculate if the line would be at 1keV
   * (if the line is at 5keV, the energy grid is multiplied by 1/5)
   */
  void shift_energy_grid_1keV(double line_energy) const {
    for (size_t ii = 0; ii < m_num_flux_bins + 1; ii++) {
      energy[ii] /= line_energy;
    }
  }

  /**
   * shift the spectrum in redshift
   */
  void shift_energy_grid_redshift(double z) const {
    if (z > 0) {
      for (size_t ii = 0; ii < m_num_flux_bins + 1; ii++) {
        energy[ii] *= (1 + z);
      }
    }
  }

  // multiply the spectrum (i.e., the flux!!) by a given value
  const XspecSpectrum &operator*=(double value) const {
    for (size_t ii = 0; ii < m_num_flux_bins + 1; ii++) {
      flux[ii] *= value;
    }
    return *this;
  }

  void multiply_flux_by(double value) const {
    for (size_t ii = 0; ii < m_num_flux_bins + 1; ii++) {
      flux[ii] *= value;
    }
  }

  // get the energy flux in ergs/cm^2/s
  [[nodiscard]] double get_energy_flux() const {
    double ener_flux = 0.0;
    for (size_t ii = 0; ii < m_num_flux_bins - 1; ii++) {
      ener_flux += flux[ii] * 0.5 * (energy[ii] + energy[ii + 1]);
    }
    return ener_flux * CONVERT_KEV2ERG;
  }

 public:
  double *energy{nullptr};
  double *flux{nullptr};
 private:
  size_t m_num_flux_bins{};

};


class DefaultSpec {
 public:
  DefaultSpec(double emin, double emax, size_t n_bins) : num_flux_bins{n_bins} {
    size_t n_energy = n_bins + 1;
    energy = new double[n_energy];
    flux = new double[n_bins];

    set_log_grid(energy, n_energy, emin, emax);
    set_input_flux();
  };

  DefaultSpec() : DefaultSpec(0.1, 1000.0, 3000)
  { };

  ~DefaultSpec() {
    delete[] energy;
    delete[] flux;
  }

  // delete copy and move assignment constructor
  DefaultSpec(const DefaultSpec& other) = delete;
  DefaultSpec& operator= (const DefaultSpec& other) = delete;

  /**
   * return the energy and flux as an XspecSpectrum
   * @return XspecSpectrum
   */
  [[nodiscard]] XspecSpectrum get_xspec_spectrum() const {
    return {energy, flux, num_flux_bins};
  }


  /* get a logarithmic grid from emin to emax with n_ener bins  */
  static void set_log_grid(double *ener, size_t n_ener, double emin, double emax) {
    for (size_t ii = 0; ii < n_ener - 1; ii++) {
      ener[ii] = static_cast<double>(ii) / (static_cast<double>(n_ener) - 1.0) * (log(emax) - log(emin)) + log(emin);
      ener[ii] = exp(ener[ii]);
    }
    ener[n_ener - 1] = emax; // separate case (otherwise it is only approx. emax, due to the log/exp functions)
  }

 public:
  double *energy{nullptr};
  double *flux{nullptr};
  const size_t num_flux_bins;

 private:
  void set_input_flux() const {
    for (size_t ii = 0; ii < num_flux_bins; ii++) {
      flux[ii] = 0.0;
    }
    flux[num_flux_bins / 2] = 1.0;
  }
};

class Spectrum {
 public:
  Spectrum(const double *_energy, size_t n_bins) : num_flux_bins{n_bins} {
    size_t n_energy = n_bins + 1;
    m_energy = new double[n_energy];
    flux = new double[n_bins];

    for (size_t ii = 0; ii < n_energy; ii++) {
      m_energy[ii] = _energy[ii];
    }

    set_input_flux_zero();
  };

  ~Spectrum() {
    delete[] m_energy;
    delete[] flux;
  }

  // delete copy and move assignment constructor
  Spectrum(const DefaultSpec &other) = delete;
  Spectrum &operator=(const DefaultSpec &other) = delete;

  /**
   * return the energy and flux as an XspecSpectrum
   * @return XspecSpectrum
   */
  [[nodiscard]] XspecSpectrum get_xspec_spectrum() const {
    return {m_energy, flux, num_flux_bins};
  }

  [[nodiscard]] double *energy() const {
    return m_energy;
  }

 public:
  const size_t num_flux_bins;
  double *flux{nullptr};

 private:
  void set_input_flux_zero() const {
    for (size_t ii = 0; ii < num_flux_bins; ii++) {
      flux[ii] = 0.0;
    }
  }

  double *m_energy{nullptr};

};

#endif //RELXILL__CPPSPECTRUM_H_
