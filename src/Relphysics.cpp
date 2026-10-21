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

#include "Relphysics.h"
#include "Rellp_Extended.h" // for tabulated energy shifts only

extern "C" {
#include "common.h"
#include "relutility.h"
}

// calculate the u^t component of the 4-velocity of a thin accretion disk
// (see Bardeen, Press, Teukolsky; 1972)
double ut_disk(double r, double a) {
  return ((r * sqrt(r) + a) /
          (sqrt(r) * sqrt(r * r - 3 * r + 2 * a * sqrt(r))));
}

double calc_proper_area_ring(double rlo, double rhi, double a) {
  double rmean = 0.5 * (rlo + rhi);

  assert(a < 1);
  assert(a >= -1);

  double rho2 = rmean * rmean + a * a; // this is rho^2
  double del = rmean * rmean - 2 * rmean + a * a;

  double area_gr = 2 * M_PI / sqrt(del);
  //  area_gr *= sqrt((rho2 * rho2 + 2 * a * a * rmean));
  area_gr *= sqrt((rho2 * rho2 - a * a * del));

  area_gr *= (rhi - rlo);

  return area_gr;
}

/**
 * @brief Blackbody photon spectrum per unit energy (NOT bin-integrated)
 *
 * Evaluates spec[i] = (E/gfac)^2 / (exp(E/(gfac*kT)) - 1) at the bin center E.
 *
 * @details Unit: photon flux density, proportional to [photons/s/cm^2/keV] emitted per unit
 * area. The constant prefactor 2*pi/(h^3 c^2) is omitted, but the temperature dependence is
 * kept (photon flux ~ T^3, energy flux ~ T^4), so spectra of different temperatures can be
 * compared/summed directly (as needed for the multi-zone disk in relxillBB).
 * To get counts/bin (Xspec convention), multiply by the bin width (e.g., via bin_integrate_flux).
 *
 * @param ener [input] energy bin boundaries in keV (size n+1)
 * @param n [input] number of bins
 * @param spec [output] photon flux density per keV (size n)
 * @param temperature [input] blackbody temperature kT in keV
 * @param gfac [input] energy shift factor; spectrum is evaluated at E/gfac (1.0 for no shift)
 */
void bbody_photon_flux_per_keV(const double *ener, int n, double *spec, double temperature, double gfac) {
  for (int ii = 0; ii < n; ii++) {
    double emean = 0.5 * (ener[ii] + ener[ii + 1]);
    spec[ii] = pow((emean / gfac), 2) / (exp(emean / gfac / temperature) - 1);
  }
}

// calculate the temperature for a given radius r and respective Rin following SS73
static double disk_temperature_alpha(double r, double Rin) {
  return pow((r / Rin), (-3. / 4)) * pow((1 - sqrt(Rin / r)), (1. / 4));
}

// temperature of a disk as used by diskbb (same radial dependence as the alpha disk for large radii)
static double disk_temperature_diskbb(double r, double Rin, double Tin) {
  return Tin * pow(r / Rin, (-3. / 4));
}

// temperature profile for the given radii; use Tmax at Rmax following Poutanan+2007
void disk_Tprofile_alpha(const double *rlo, const double *rhi, double *temp, int n, double Rin, double Tin) {
  // use Tin as Tmax
  double Rmax = pow((1.5), (4. / 5)) * Rin; // Poutanen+2007
  double K = Tin / disk_temperature_alpha(Rmax, Rin);

  double RminGrid = 0.5 * (rlo[0] + rhi[0]);
  assert(RminGrid >= Rin);

  for (int ii = 0; ii < n; ii++) {
    temp[ii] = K * disk_temperature_alpha(0.5 * (rlo[ii] + rhi[ii]), Rin);
  }
}

// temperature profile of diskbb for the given radii: T(r) = Tin (r/Rin)^(-3/4), with Tin at the inner edge of the
// disk (as for diskbb); the temperature of each zone is evaluated at its mean radius
void disk_Tprofile_diskbb(const double *rlo, const double *rhi, double *temp, int n, double Tin) {
  double Rin = rlo[0]; // inner edge of the disk
  for (int ii = 0; ii < n; ii++) {
    temp[ii] = disk_temperature_diskbb(0.5 * (rlo[ii] + rhi[ii]), Rin, Tin);
  }
}

// radial dependence of the flux emitted by a relativistic thin disk around a Kerr black hole with zero torque at its
// inner edge Rin (Novikov & Thorne 1973), with radii in units of r_g = GM/c^2 and Rin >= r_isco. The flux is
// F = Mdot/(4 pi r) f, with the closed form of f from Page & Thorne (1974, ApJ 191, 499, eq. 15n), see also Li et al.
// (2005, ApJS 157, 335, eqs. D11 and D13, the paper of the Xspec model kerrbb):
//   f = 3/2 / (x^2 (x^3 - 3x + 2a)) * [x - x0 - 3/2 a ln(x/x0)
//                                      - sum_i 3(x_i-a)^2/(x_i(x_i-x_j)(x_i-x_k)) ln((x-x_i)/(x0-x_i))]
// with x = sqrt(r), x0 = sqrt(Rin), and x_i the roots of x^3 - 3x + 2a = 0. Only the radial dependence F ~ f/x^2 is
// returned (without the constant factors), which approaches r^-3 for large radii.
static double disk_flux_novikov_thorne(double r, double Rin, double spin) {
  const double x = sqrt(r);
  const double x0 = sqrt(Rin);
  const double x1 = 2 * cos((acos(spin) - M_PI) / 3);
  const double x2 = 2 * cos((acos(spin) + M_PI) / 3);
  const double x3 = -2 * cos(acos(spin) / 3);

  const double integ = x - x0 - 1.5 * spin * log(x / x0)
                       - 3 * pow(x1 - spin, 2) / (x1 * (x1 - x2) * (x1 - x3)) * log((x - x1) / (x0 - x1))
                       - 3 * pow(x2 - spin, 2) / (x2 * (x2 - x1) * (x2 - x3)) * log((x - x2) / (x0 - x2))
                       - 3 * pow(x3 - spin, 2) / (x3 * (x3 - x1) * (x3 - x2)) * log((x - x3) / (x0 - x3));

  // the integral vanishes at r = Rin (zero torque), a tiny negative value can only occur from rounding
  return fmax(integ / (pow(x, 4) * (x * x * x - 3 * x + 2 * spin)), 0.0);
}

// radius of the maximal flux of the Novikov-Thorne profile (golden-section search between Rin and 4 Rin; the maximum
// is at 1.3-1.6 Rin for an inner edge at the ISCO)
static double radius_max_flux_novikov_thorne(double Rin, double spin) {
  const double golden = (sqrt(5.) - 1) / 2;
  double rlo = Rin;
  double rhi = 4 * Rin;
  for (int ii = 0; ii < 100; ii++) {
    const double r1 = rhi - golden * (rhi - rlo);
    const double r2 = rlo + golden * (rhi - rlo);
    if (disk_flux_novikov_thorne(r1, Rin, spin) > disk_flux_novikov_thorne(r2, Rin, spin)) {
      rhi = r2;
    } else {
      rlo = r1;
    }
  }
  return 0.5 * (rlo + rhi);
}

// Novikov-Thorne temperature profile T ~ F^(1/4) for the given radii, with zero torque at the inner edge Rin and the
// maximal temperature Tin (Tmax, as for the alpha disk); the temperature of each zone is evaluated at its mean radius
void disk_Tprofile_novikov_thorne(const double *rlo,
                                  const double *rhi,
                                  double *temp,
                                  int n,
                                  double Rin,
                                  double Tin,
                                  double spin) {
  assert(Rin >= kerr_rms(spin) * (1 - 1e-6));
  const double flux_max = disk_flux_novikov_thorne(radius_max_flux_novikov_thorne(Rin, spin), Rin, spin);

  for (int ii = 0; ii < n; ii++) {
    temp[ii] = Tin * pow(disk_flux_novikov_thorne(0.5 * (rlo[ii] + rhi[ii]), Rin, spin) / flux_max, 0.25);
  }
}

VecD get_tprofile(const VecD &rlo,
                  const VecD &rhi,
                  const int nrad,
                  double Rin,
                  double Tin,
                  double spin,
                  int type,
                  int *status) {
  VecD tprofile(nrad);

  if (type == TPROFILE_ALPHA) {
    disk_Tprofile_alpha(rlo.data(), rhi.data(), tprofile.data(), nrad, Rin, Tin);
    if (is_debug_run()) { printf(" using temperature profile ALPHA\n"); }
  } else if (type == TPROFILE_DISKBB) {
    disk_Tprofile_diskbb(rlo.data(), rhi.data(), tprofile.data(), nrad, Tin);
    if (is_debug_run()) { printf(" using temperature profile DISKBB\n"); }
  } else if (type == TPROFILE_NOVIKOV_THORNE) {
    disk_Tprofile_novikov_thorne(rlo.data(), rhi.data(), tprofile.data(), nrad, Rin, Tin, spin);
    if (is_debug_run()) { printf(" using temperature profile NOVIKOV-THORNE\n"); }
  } else {
    RELXILL_ERROR(" failed getting the temperature profile ", status);
    printf("    reason: type=%i is unknown\n", type);
  }

  return tprofile;
}


/*** we calculate the disk density from  Shakura & Sunyaev (1973)
 *    - for zone A as describe in their publication,  formula (Eq 2.11)
 *    - only the radial dependence is picked up here  (viscosity alpha=const.)
 *    - normalized such that dens(rms) = 1
 *    - rms is given in units of Rg
 *                                                               ***/
double density_ss73_zone_a(double radius, double rms) {
  return pow((radius / rms), (3. / 2)) * pow((1 - sqrt(rms / radius)), -2);
}

/**
 * @brief wrapper around density_ss73_zone_a, to ensure the value of the function (meaning the density) is
 * dens=1 at where it is minimal (at 25/9*rin)
 * @param radius
 * @param rms
 * @return
 */
double density_ss73_zone_a_normalized(double radius, double rms) {
  const double rad_density_min = (25. / 9.) * rms;
  const double density_min = density_ss73_zone_a(rad_density_min, rms);

  return density_ss73_zone_a(radius, rms) / density_min;
}


double relat_abberation(double del, double beta) {
  return acos((cos(del) - beta) / (1 - beta * cos(del)));
}

/**
 * @brief calculates the energy shift from a lamp post source at h to infinity
 */
double calc_g_inf(double height, double a) {
  return sqrt(1.0 - (2 * height / (height * height + a * a)));
}

double calc_cos_theta_ring_source(const double height, const double ring_radius, const double spin) {
  if (ring_radius == 0.0) {
    return 1.0; // cos_theta = 1 on the rotational axis
  }
  if (spin == 0.0) {
    // special case of zero spin
    return height / sqrt(height * height + ring_radius * ring_radius);
  } // general case, solution of equation x = sqrt(h^2 * tan(theta)^2 + a^2 * sin(theta)^2) for cos(theta)
  return sqrt((sqrt(4.0 * spin * spin * height * height +
                    (height * height + ring_radius * ring_radius - spin * spin) *
                    (height * height + ring_radius * ring_radius - spin * spin)) +
               spin * spin - height * height - ring_radius * ring_radius) / (2.0 * spin * spin));
}

double calc_spherical_radius_ring_source(const double height, const double ring_radius, const double spin) {
  double cos_theta = calc_cos_theta_ring_source(height, ring_radius, spin);
  return height / cos_theta; // spherical radius by definition
}

double calc_extent_from_sph_radius_costheta(const double r_sph, const double costheta, const double spin) {
  return sqrt(1.0 - costheta * costheta) * sqrt(r_sph * r_sph + spin * spin);
}

double calc_extent_from_sph_radius_primary_source(const double height, const double r_sph, const double spin) {
  if (abs(height - r_sph) < 1e-4) {
    return 0.0; // it is lamp post
  }
  assert(height <= r_sph && ("Height cannot be more than spherical radius"));
  return sqrt((r_sph * r_sph - height * height) * (1.0 + spin * spin / (r_sph * r_sph)));
}

double calc_height_from_rcostheta(const double r_sph, const double costheta) {
  return r_sph * costheta;
}

void convert_hx_to_rtheta(const double height,
                          const double ring_radius,
                          const double spin,
                          double *r_sph,
                          double *theta) {
  *r_sph = calc_spherical_radius_ring_source(height, ring_radius, spin);
  *theta = acos(calc_cos_theta_ring_source(height, ring_radius, spin)) / CONVERT_DEG2RAD; //return in degrees!
}


/* get RMS (ISCO) for the Kerr Case */
double kerr_rms(double a) {
  //	 accounts for negative spin
  double sign = 1.0;
  if (a < 0) {
    sign = -1.0;
  }

  double Z1 = 1.0 + pow(1.0 - a * a, 1.0 / 3.0) * (pow(1.0 + a, 1.0 / 3.0) + pow(1.0 - a, 1.0 / 3.0));
  double Z2 = sqrt((3.0 * a * a) + (Z1 * Z1));

  return 3.0 + Z2 - sign * sqrt((3.0 - Z1) * (3.0 + Z1 + (2 * Z2)));
}

/* get the rplus value (size of the black hole event horizon */
double kerr_rplus(double a) {
  return 1.0 + sqrt(1.0 - a * a);
}

/** calculate the doppler factor for a moving primary source **/
double doppler_factor(double del, double bet) {
  return sqrt(1.0 - bet * bet) / (1.0 + bet * cos(del));
}

/** calculates g = E/E_i in the lamp post geometry (see, e.g., 27 in Dauser et al., 2013, MNRAS) **/
double gi_potential_lp(double r, double a, double h, double bet, double del) {
  /** ! calculates g = E/E_i in the lamp post geometry
    ! (see, e.g., page 48, Diploma Thesis, Thomas Dauser) **/
  double ut_d = ((r * sqrt(r) + a) / (sqrt(r) * sqrt(r * r - 3 * r + 2 * a * sqrt(r))));
  double ut_h = sqrt((h * h + a * a) / (h * h - 2 * h + a * a));

  double gi = ut_d / ut_h;

  // check if we need to calculate the additional factor for the velocity
  if (fabs(bet) < 1e-6) {
    return gi;
  }

  double gam = 1.0 / sqrt(1.0 - bet * bet);

  // get the sign for the equation
  double sign = 1.0;
  if (del > M_PI / 2) {
    sign = -1.0;
  }

  // ** Adam, 28.7.2021: **
  // So given the expression for Carter’s q:
  //q^2 = \sin^2\delta (h^2+a^2)^2/\Delta_h - a^2,
  //it is clear that
  //q^2 + a^2 = \sin^2\delta (h^2+a^2)^2/\Delta_h
  //Therefore:
  //(h^2+a^2)^2 - \Delta_h (q^2+a^2) = (h^2+a^2)^2 - \sin^2\delta (h^2+a^2)^2
  //= (h^2+a^2)^2 \cos^2\delta
  //
  //Therefore:
  //\sqrt{ (h^2+a^2)^2 - \Delta_h (q^2+a^2) } / (h^2+a^2) = \cos\delta
  //
  //Therefore your equation (27) becomes:
  //glp = glp(beta=0) / { \gamma [ 1 -/+ \beta\cos\delta ] }

  double delta_eq = h * h - 2 * h + a * a;
  double q2 = (pow(sin(del), 2)) * (pow((h * h + a * a), 2) / delta_eq) - a * a;

  double beta_fac = sqrt(pow((h * h + a * a), 2) - delta_eq * (q2 + a * a));
  beta_fac = gam * (1.0 + sign * beta_fac / (h * h + a * a) * bet);

  return gi / beta_fac;
}

/** similar to gi_potential_lp below, but for ring or slab source **/
double gi_potential_ext(const relParam *param, const double r) {
  const double g_sd = get_energyboost_ext_source_disk(param, r, 1.0); // gamma == 1.0
  return g_sd;
}

/**
 * @brief calculate the doppler factor from source to the observer
 * we currently assume that delta_obs and incl are the same (so no GR light-bending
 * is accounted for).
 * @param incl [in rad]
 * @param beta [in v/c]
 * @return
 */
double doppler_factor_source_obs(const relParam *rel_param) {
  //update?
  // note that the angles incl and delta_obs are defined the opposite way
  double delta_obs = M_PI - rel_param->incl;
  assert(delta_obs > M_PI / 2);
  assert(delta_obs < M_PI);

  // glp = D*glp(beta=0) = glp(beta=0) / { \gamma [ 1 -/+ \beta\cos\delta ] }
  // Adam: we use the minus sign for δ > π/2 and the plus sign for δ < π/2 to get
  // -> as always δ > π/2 in this case, we directly can use the minus sign
  double doppler = doppler_factor(delta_obs, rel_param->beta);
  // assert(doppler >= 1-1e-8);

  return doppler;
}

/**
 * @brief calculate the energy shift from the LP to the observer, also taking
 * a potential velocity of the source into account
 *
 * note that this routine uses  the function "doppler_factor_source_obs", which
 * sets δ_obs = incl, which is not fully true in GR, however, within 1% accuracy
 */
double energy_shift_source_obs(const relParam *rel_param) {
  double g_inf_0 = 1;
  // if we do NOT have a geometry, the energy shift is just 1
  // (as without geometry we can not define an energy shift)
  if (rel_param == nullptr || !is_any_primary_source(rel_param->emis_type)) {
    return g_inf_0;
  }

  if (is_ring_primary_source(rel_param->emis_type) || is_slab_primary_source(rel_param->emis_type)) {
    g_inf_0 = get_energyboost_ext_source_obs(rel_param, 1.0);
  } else {
    // lamp post
    g_inf_0 = calc_g_inf(rel_param->height, rel_param->a);
  }

  if (rel_param->beta < 1e-4) {
    return g_inf_0;
  } else {
    return g_inf_0 * doppler_factor_source_obs(rel_param);
  }
}

/**
 * @brief calculate the energy shift from the LP to the disk zone with mean radius radius_disk of photon trajectories
 * emitted under del_emit (in the frame of the non-moving source); del_emit is only important for beta>0
 * @details AN: currently used in IonGradient.cpp
 * The connection for which del_emit which radius radius_disk is hit, is stored in the LP table structure
 */
double energy_shift_source_disk(const relParam *rel_param, double radius_disk, double del_emit) {
  // if we do NOT have a geometry, the energy shift is just 1
  // (as without geometry we can not define an energy shift)
  if (rel_param == nullptr || not is_any_primary_source(rel_param->emis_type)
  )
  {
    return 1;
  }
  if (is_ring_primary_source(rel_param->emis_type) || is_slab_primary_source(rel_param->emis_type)) {
    double gshift = get_energyboost_ext_source_disk(rel_param, radius_disk, 1.0);
    return gshift;
  } else {
    return gi_potential_lp(radius_disk, rel_param->a, rel_param->height, rel_param->beta, del_emit);
  }
}

/**
 * @brief calculate the geometric/Newtonian emissivity at a given radius
 * @details the emissivity profiles are normalized such that at large height and
 * radii the GR profiles converge towards this definition
 * @return
 */
double calc_emissivity_newton(const relParam *param, const double radius) {
  double emis;
  if (is_ring_primary_source(param->emis_type)) {
    double h = param->r_src * cos(param->theta * CONVERT_DEG2RAD);
    double x = param->r_src * sin(param->theta * CONVERT_DEG2RAD);
    emis = calc_ring_emissivity_1d_newton(h, x, radius);
  } else if (is_slab_primary_source(param->emis_type)) {
    // AN: this warning is here because this function is used only for relxillAlpha type models,
    // so I do not expect to use this function yet
    printf("calc_emissivity_newton: slab emissivity may be inaccurate! \n");
    // AN: here I use the trick that slab emissivity at a given radius has the largest contribution from the ring
    // with exactly this radius (no GR)
    emis = calc_ring_emissivity_1d_newton(param->height, radius, radius);
  } else {
    emis = calc_lp_emissivity_newton(param->height, radius);
  }
  return emis;
}

/**
 * @brief calculate the geometric/Newtonian emissivity for the LP geometry
 * @details the emissivity profiles are normalized such that at large height and
 * radii the GR profiles converge towards this definition
 * @param h
 * @param r
 * @return
 */
double calc_lp_emissivity_newton(double h, double r) {
  double emis = pow(1. / (pow((r / h), 2) + 1), (3. / 2));
  emis /= (2 * M_PI * h * h);

  return emis;
}

/**
 * @brief calculate the Newtonian 1d emissivity for the RING geometry
 * @details This is only a proxy of the flux produced by the ring in newtonian limit,
 * for proper flux we need a 3d raytracing in flat space
 * However, this proxy is good enough for the current use in the code
 * @param h Height of the ring
 * @param x Radius of the ring (cylindrical)
 * @param r
 * @return
 */
double calc_ring_emissivity_1d_newton(double h, double x, double r) {
  double emis = pow(1.0 / ((r - x) * (r - x) / h / h + 1.0), (3. / 2));
  double bracket = sqrt(1.0 + x * x / h / h) + x / h;
  emis /= (2.0 * M_PI * h * h * bracket);

  return emis;
}

/**
 * @brief calculate the fluxboost from the (moving) primary source to the accretion disk
 * at radius rad, following g^Gamma*Doppler_factor^2 (see Dauser et al., 2013)
 * @details AN: used in Rellp to normalize flux of the LP model only. For ext models we do this during interpolation
 * Which is a bit inconsistent, so TODO (maybe) unify boost calculation for lp and ext models
 * @param rad
 * @param del_emit
 * @param a
 * @param height
 * @param gamma
 * @param beta
 * @return
 */
double calc_fluxboost_source_disk(double rad, double del_emit, double a, double height, double gamma, double beta) {
  // then energy shift factor gi^gamma (see Dauser et al., 2013)
  double boost = pow(gi_potential_lp(rad, a, height, beta, del_emit), gamma);

  // take the beaming of the jet into account (see Dauser et al., 2013)
  if (beta > 1e-6) {
    boost *= pow(doppler_factor(del_emit, beta), 2);
  }
  return boost;
}

double kerr_to_flat_area_ratio(double r, double a) {
  // the value below is obtained directly from ratio A_Kerr / A_flat by all possible simplifications
  return sqrt((r * r * r + a * a * r + 2.0 * a * a) /
              (r * r * r + a * a * r - 2.0 * r * r));
}

double gamma_correction(double r, double a) {
  // the value is obtained from the "Bardeen"-like form of the Lorentz factor by a few simplifications
  return (r * sqrt(r) + abs(a)) * sqrt((r * r - 2.0 * r + a * a) / (r * r +
                                                                    2.0 * a * sqrt(r) - 3.0 * r) / (
                                         r * r * r + a * a * r + 2.0 * a * a));
}