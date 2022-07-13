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
#ifndef IONGRADIENT_H_
#define IONGRADIENT_H_

#include "vector"
#include <memory>
#include <array>

extern "C" {
#include "relutility.h"
#include "relphysics.h"
}



class RadialGrid{

 public:
  RadialGrid(double rmin, double rmax, int nzones, double h) : num_zones(nzones) {
    radius = calculate_radial_grid(rmin, rmax, nzones, h) ;
  }

  const int num_zones;
  const double* radius; // has length num_zones+1

  ~RadialGrid(){
      delete[] radius;
  }

 private:
  static const double* calculate_radial_grid(double rmin, double rmax, int nzones, double h);
};



class IonGradient{

 public:
  IonGradient(const RadialGrid& _radial_grid , const int ion_grad_type)
  : radial_grid{_radial_grid},  // radius is of length nzones+1
    m_nzones{_radial_grid.num_zones},
    m_ion_grad_type{ion_grad_type}
  {
    lxi = new double[m_nzones];
    dens = new double[m_nzones];
    m_rmean = new double[m_nzones];

    for ( int ii=0; ii<m_nzones; ii++){
      m_rmean[ii] = 0.5*(radial_grid.radius[ii]+radial_grid.radius[ii+1]);
      dens[ii] = 15.0;
      lxi[ii] = 0.0;
    }

    if (m_nzones>1){
      assert(m_rmean[0] < m_rmean[1]); // make sure we have an ascending radial grid
    }
  };

  ~IonGradient(){
    delete[] lxi;
    delete[] irradiating_flux;
    delete[] del_emit;
    delete[] dens;
    delete[] m_rmean;
  }

  const RadialGrid& radial_grid;
  double* lxi{nullptr};
  double *irradiating_flux{nullptr};
  double *del_emit{nullptr};
  double *dens{nullptr};

  [[nodiscard]] int nzones() const {
    return m_nzones;
  }

  void calculate(const emisProfile &emis_profile, xillParam *xill_param);

  double get_ecut_disk_zone(const relParam *rel_param, double ecut_primary, int izone) const;

  void write_to_file(const char *fout) const;

 private:
  double *m_rmean;
  const int m_nzones;
  const int m_ion_grad_type;

  void calc_ion_grad_alpha(const emisProfile &emis_profile, double param_xlxi0, double param_density);

  void calc_ion_grad_pl(double xlxi0, double xindex, double inputval_dens);

  void calc_zones_delta_emit(const emisProfile &emis_profile);

};



#endif
