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
#ifndef IONGRADIENT_H_
#define IONGRADIENT_H_

#include "vector"

extern "C" {
#include "relutility.h"
#include "relphysics.h"
}


class IonGradient{

 public:
  IonGradient(double* radius, const int nzones , const int ion_grad_type)
  : m_radius{radius},  // radius is of length nzones+1
  m_nzones{nzones},
    m_ion_grad_type{ion_grad_type}
  {
    lxi = new double[nzones];
    fx = new double[nzones];
    del_emit = new double[nzones];
    dens = new double[nzones];

    m_rmean = new double[nzones];
    for ( int ii=0; ii<nzones; ii++){
      m_rmean[ii] = 0.5*(m_radius[ii]+m_radius[ii+1]);
    }

    if (m_nzones>1){
      assert(m_rmean[0] < m_rmean[1]); // make sure we have an ascending radial grid
    }
  };

  ~IonGradient(){
    delete lxi;
    delete fx;
    delete del_emit;
    delete dens;
    delete m_rmean;
  }

  double* lxi{nullptr};
  double* fx{nullptr};
  double* del_emit{nullptr};
  double* dens{nullptr};

  [[nodiscard]] int nzones() const{
    return m_nzones;
  }

  void calculate(relParam* rel_param, xillParam* xill_param);

  void write_to_file( const char* fout);

 private:
  double* m_radius;
  const int m_nzones;
  const int m_ion_grad_type;

  double* m_rmean;

  void calc_ion_grad_alpha(relParam* rel_param, double param_xlxi0, double param_density);

  void calc_ion_grad_pl(double xlxi0, double xindex);

};

#endif
