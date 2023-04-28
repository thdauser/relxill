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

#include "Relphysics.h"
#include "PrimarySource.h"

extern "C" {
#include "relutility.h"
}

class RadialGrid {

 public:
  RadialGrid(double rmin, double rmax, int nzones, double h)
      : radius{calculate_radial_grid(rmin, rmax, nzones, h)} {
    assert(nzones == static_cast<int>( num_zones() ) );
  }

  [[nodiscard]] size_t num_zones() const {
    return radius.size() - 1;
  }

  std::vector<double> radius{}; // has length num_zones+1

 private:
  static std::vector<double> calculate_radial_grid(double rmin, double rmax, int nzones, double h);
};



class IonGradient{

 public:
  IonGradient(const RadialGrid &_radial_grid, int ion_grad_type, double ion_grad_index)
      : radial_grid{_radial_grid},  // radius is of length nzones+1
        m_nzones{static_cast<int>(_radial_grid.num_zones())},
        m_ion_grad_type{ion_grad_type},
        m_ion_grad_index{ion_grad_index}
  {
    lxi = new double[m_nzones];
    dens = new double[m_nzones];
    m_rmean = new double[m_nzones];
    //  scaling of the reflection spectrum for the zone (can be necessary if kTe/Ecut is shifted by energy)
    m_energy_shift_source_disk = new double[m_nzones];

    for ( int ii=0; ii<m_nzones; ii++){
      m_rmean[ii] = 0.5*(radial_grid.radius[ii]+radial_grid.radius[ii+1]);
      dens[ii] = 15.0;
      lxi[ii] = 0.0;
      m_energy_shift_source_disk[ii] = 1.0;
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
    delete[] m_physical_irrad_flux;
    delete[] m_rmean;
    delete[] m_energy_shift_source_disk;
  }

  RadialGrid radial_grid;
  double *lxi{nullptr};
  double *irradiating_flux{nullptr};
  double *del_emit{nullptr};
  double *dens{nullptr};
  double *m_energy_shift_source_disk{nullptr};

  [[nodiscard]] int nzones() const {
    return m_nzones;
  }

  void calculate_gradient(const emisProfile &emis_profile, const PrimarySourceParameters &primary_source_params);

  [[nodiscard]] xillTableParam **get_xill_param_zone(const xillTableParam *primary_source_xilltab_params) const;

  double get_ecut_disk_zone(const relParam *rel_param, double ecut_primary, int izone) const;

  void write_to_file(const char *fout) const;

  static double calculate_lxi_max_from_distance(const emisProfile &emis_profile,
                                                const PrimarySourceParameters &primary_source_params,
                                                double log_density_rin);

 private:
  double *m_rmean;
  int m_nzones;
  int m_ion_grad_type;
  double m_ion_grad_index; // only used for the PL gradient
  double m_physical_emis_normalization_factor = 0; // is set in "calculate_lxi_max_from_distance"
  double *m_physical_irrad_flux{nullptr};

  void calc_ion_grad_alpha(const emisProfile &emis_profile, double param_xlxi0, double param_density);

  void calc_ion_grad_pl(double xlxi0, double xindex, double inputval_dens);

  void calc_energy_shift_from_source_to_disk(const relParam *rel_param) const;

  void set_del_emit_for_each_zone(const emisProfile &emis_profile);

  void calculate_physical_emissivity(const PrimarySourceParameters &primary_source, const emisProfile &emis_profile);

};



#endif
