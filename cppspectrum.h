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




class CppSpectrum {

 public:
  CppSpectrum(const Array &_energy, Array &_flux)
  // need to allocate the energy grid, as Xspec requires it to be constant and we may shift it in energy
      : m_ener{_energy},
        m_flux{_flux} {
    assert(_energy.size()
               == _flux.size()); // Xspec Local Model Convention, although last flux bin will be chopped (see below)
    // recommended by Xspec (https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html)
    //  m_flux.resize(_flux.size() - 1);


    // double arrays have the same size as the input Array
    m_ener_double = new double[m_ener.size()];
    m_flux_double = new double[m_flux.size()];

  };

  ~CppSpectrum() {
    delete[] m_ener_double;
    delete[] m_flux_double;
  }

  [[nodiscard]] Array energy() const {
    return m_ener;
  }

  [[nodiscard]] Array flux() const {
    return m_flux;
  }

  [[nodiscard]] size_t nener_bins() const {   // array holds nener+1 bins, as bin_lo and bin_hi are combined
    return m_ener.size() - 1;
  }

  /**
   * get the energy as double array, it will be allocated at the first call and always
   * all values copied from the energy() function
   */
  [[nodiscard]] double *energy_double() {
    convert_array2double(m_ener, m_ener_double);
    return m_ener_double;
  }

  /**
   * get the flux as double array, it will be allocated at the first call and always
   * all values copied from the flux() function
   */
  [[nodiscard]] double *flux_double() {
    convert_array2double(m_flux, m_flux_double);
    return m_flux_double;
  }

  /**
   * copy values from the double array to the Array
   * @throws std::bad_alloc if the double arrays are not allocated
   */
  void copy_doubleArrays2array() {
    if (m_flux_double && m_ener_double) {
      write_double2array(m_ener, m_ener_double);
      write_double2array(m_flux, m_flux_double);
    } else {
      puts(" *** relxill-error: double-arrays of the spectrum not allocated ");
      throw std::bad_alloc(); //" copying of double arrays in Spectrum failed, as they are not allocated ");
    }

  }

  /**
   * shift the spectrum such that we can calculate the line for 1 keV
   */
  void shift_energy_grid_1keV(double line_energy, double z) {
    m_ener *= (1 + z) / line_energy;
  }

 private:
  Array m_ener;
  Array m_flux;
  double *m_ener_double{nullptr};
  double *m_flux_double{nullptr};

  /**
 * Converts a C++ container "Array" into a double-array, which is necessary for C-code.
 * The routine needs to copy the array and all single values. If ptr_double_array==nullptr
 * it will be created, otherwise it is assumed that it is correctly allocated to array.size()
 * @param array
 * @return ptr_double_array
 *
 */
  static void convert_array2double(const Array &array, double *double_array) {

    assert(double_array);
    for (int ii = 0; ii < array.size(); ii++) {
      double_array[ii] = array[ii];
    }

  }

  /**
   * write values from the double-array into the Array
   * - potentially dangerous, if used outside this scope (i.e., if the double array was not
   *   created from the Array) and therefore should stay private
   */
  static void write_double2array(Array &array, const double *double_array) {
    for (size_t ii = 0; ii < array.size(); ii++) {
      array[ii] = double_array[ii];
    }
  }

};

#endif //RELXILL__CPPSPECTRUM_H_
