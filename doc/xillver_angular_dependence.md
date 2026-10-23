# Xillver Angular Dependence: Intensity and Flux {#xillver_angular_dependence}

The xillver tables contain the reflected spectrum for a set of emission angles as **intensities**. Whether a spectrum
is an intensity or a flux decides which conversion factor is applied, and where. This page gives the convention, the
code paths that convert intensities into fluxes, and the rule that follows for everything that is handed to the
relativistic convolution. The absolute normalization of the reflected spectrum (\f$N_\mathrm{xill}\f$) is described in
\ref returnrad_normalizations.

## Intensities in the xillver table

Each spectrum of the table at inclination \f$i\f$ is the intrinsic intensity \f$I (\mu_i)\f$ emerging from the slab
under the emission angle \f$\mu_i = \cos i\f$, for the illumination the table is normalized to (Dauser et al. 2016).
The inclinations (`xillSpec::incl`, in degrees) are the bin centres of a grid uniform in \f$\mu\f$.

## From intensity to flux

The flux of a semi-infinite slab is the angle average of the intensities (Garcia et al. 2014, Eq. 19):

\f[
F \;=\; \int_0^1 \frac{1}{2}\, \mu\, I (\mu)\, d\mu \;\approx\; \sum_i \frac{\mu_i}{2\, n_\mathrm{incl}}\, I_i
\f]

The weights \f$\mu_i / (2 n_\mathrm{incl})\f$ of the table inclinations (\ref get_semi_infinite_slab_flux_weights)
sum to \f$\int_0^1 \frac{1}{2} \mu\, d\mu = 0.25\f$ (to \f$10^{-7}\f$ for the table grid, as the midpoint rule
in \f$\mu\f$ is exact for this integrand). For an intensity that does
not depend on the emission angle (isotropic emitter), this gives

\f[
F \;=\; 0.25\, I \qquad\Longleftrightarrow\qquad I \;=\; \frac{F}{0.25}
\f]

(`SEMI_INFINITE_SLAB_FLUX_ISOTROPIC` in `src/Xillspec.h`). The factor 0.25 instead of \f$\pi\f$ comes from the units
of the table (the illumination is normalized to \f$n\xi/4\pi\f$). Only ratios within this convention enter the model.

## Where the conversion is applied

| path                              | weights of the table inclinations                                                           | result                         | used for                                                                                                                  |
|-----------------------------------|---------------------------------------------------------------------------------------------|--------------------------------|---------------------------------------------------------------------------------------------------------------------------|
| non-relativistic, one inclination | \f$\frac{1}{2}\cos i\f$ (`norm_xillver_spec()`)                                             | flux under inclination \f$i\f$ | `xillver` models                                                                                                          |
| non-relativistic, angle-averaged  | \f$\mu_i/(2 n_\mathrm{incl})\f$ (\ref calc_angle_averaged_xillver_flux)                     | angle-averaged flux            | correction factors of the coronal returning radiation (`get_xillver_fluxcorrection_factors()`), energy-conservation tests |
| relativistic                      | GR emission-angle distribution `rel_cosne->dist[izone]` (sum = 1, \ref calc_xillver_angdep) | angle-mixed **intensity**      | `relxill`, `relxillBB`: the convolution with the line profile converts it into the observed flux                          |

**Never both.** The relativistic convolution with the line profile of a zone, \f$P_i\f$ (`relProfile->flux[izone]`,
`convolveSpectrumFFTNormalized()`), integrates the intensity over the projected area of the zone on the sky of the
observer, i.e., it carries the conversion from intensity to observed flux. Applying the slab factor before the
convolution would count it twice (see the comment in Garcia et al. 2014, Eq. 19).

This can be seen directly in the normalization of the line profiles of `relxillBB` (not renormalized, see below):
at large radii, where the relativistic effects are small, the number of photons of the line profile of a zone,
\f$\sum_E P_i\f$, is \f$\frac{1}{2}\cos i\f$ times the (flat) area of the zone, i.e., exactly the slab factor of the
observer inclination \f$i\f$ (`relxillBB`, outermost zone at \f$r \approx 700\f$):

| inclination [deg]                  | 5     | 20    | 40    | 60    | 80    |
|------------------------------------|-------|-------|-------|-------|-------|
| \f$\frac{1}{2}\cos i\f$            | 0.498 | 0.470 | 0.383 | 0.250 | 0.087 |
| \f$\sum_E P_i\f$ / area of the zone | 0.497 | 0.468 | 0.382 | 0.250 | 0.088 |

Averaged over the observer inclination (uniform in \f$\cos i\f$), this is 0.25 again. Closer to the black hole, the
relativistic effects reduce the number of observed photons (e.g., to 0.10 at \f$r \approx 3.4\f$ for 40 deg).

For the relativistic models whose line profile is renormalized to 1 (`relxill` with a broken power-law emissivity,
`renorm_relline_profile()`), the renormalization removes the absolute conversion of the profile. It is therefore
re-applied explicitly with the slab factor \f$\frac{1}{2}\cos i_\mathrm{obs}\f$ of the observer inclination. The
profiles of `relxillBB` are not renormalized (`do_renorm_model()` returns 0 for the returning-radiation models).

## Rule: the relativistic convolution takes intensities

**Every spectrum that is convolved with the line profile of a zone must be an intensity in xillver table units.**
Spectra calculated as a flux, which are emitted isotropically, are converted with \f$I = F / 0.25\f$ before the
convolution. The spectra convolved in \ref relxillbb_kernel are

- **Xillver reflection** (\ref calc_zone_xillver_reflected_intensity): \f$N_\mathrm{xill} (F_\mathrm{ret}) \sum_k
  \mathrm{dist}_k\, I_k\f$, the angle-mixed intensity of the table, scaled to the returning flux. Already an
  intensity.
- **Mirrored returning radiation** (\ref calc_zone_mirrored_intensity, test path selected with
  `RELXILL_RRAD_BBODY_MIRRORED=1`): \f$I_\mathrm{mirror} = F_\mathrm{ret} / 0.25\f$, see below. A mirror conserves
  the energy exactly: the angle-averaged flux of this intensity is \f$F_\mathrm{ret}\f$.
- **Direct disk emission** (`calc_disk_emission_observed()`): \f$I_\mathrm{pri} = F_\mathrm{pri} / 0.25\f$, as the
  blackbody surface of the disk emits isotropically, see "Direct disk emission" below.

### What "mirrored" means

The mirror is defined by three steps, of which only the last one concerns the angular dependence of the emission:

1. **Incident flux.** The returning flux \f$F_\mathrm{ret}\f$ (`specRet`) of a zone already contains the angles of
   incidence: the calculation of the returning radiation weights each ray that hits the disk with the cosine of its
   angle of incidence, which gives the flux the zone actually receives. Nothing about the incident angles is lost or
   neglected in the mirror.
2. **Spectral shape.** "Mirrored" means that the emitted spectrum has the same shape as the incident one, i.e., there
   is no reflection physics (no absorption, no lines, no Compton scattering), and the full incident flux is
   re-emitted.
3. **Re-emission.** In between, a reprocessing in the slab is assumed such that the radiation is emitted with an
   isotropic intensity, i.e., independent of the emission angle. The emitted flux then follows the same
   \f$\frac{1}{2}\cos\theta\f$ dependence on the emission angle as the xillver intensities (the definition of the
   flux above), and therefore the same conversion \f$F = 0.25\, I\f$. This is what makes the mirror directly
   comparable with the xillver reflection.

## Diagnostic output: all spectra are fluxes

The diagnostic files (`RELXILL_WRITE_FILES=1`) contain only fluxes; the files, their content, and the conversion of
each are listed in \ref returnrad_disk_emission_files. In the rest frame, the fluxes are **angle-averaged** over the
emission angles. This only matters for the xillver reflection: the primary and the mirror emit isotropically (flux =
0.25 \f$\times\f$ intensity under every emission angle), whereas `rframe-reflect` is the average of the
angle-dependent xillver intensities with the slab flux weights (\f$B_\mathrm{slab}\f$), not the angle-mixed intensity
that is convolved for the observer (\f$A_\mathrm{mix}\f$, which depends on the observer inclination and has no direct
meaning as a flux).

Checked for the default parameters (40 deg, \f$b = -1\f$):

- rest frame: `rframe-mirror-reflect` / `rframe-bbodyRet` = 1 in every bin. The energy flux (0.1-1000 keV) of
  `rframe-reflect` is 0.36-0.46 of `rframe-bbodyRet` per zone, 0.43 for the whole disk (\f$B_\mathrm{slab}\f$ of the
  test `[angleweight]`), i.e., below the one of the mirror.
- observed vs. rest frame: the ratio of the energy flux observed and in the rest frame of each zone is identical for
  `bbodyEmit` and `mirror-reflect` (to \f$4 \times 10^{-4}\f$), as both are isotropic emitters. For `reflect`, it
  differs by the angular dependence of the xillver intensity (\f$A_\mathrm{mix} / (4 B_\mathrm{slab})\f$, see below).

## Comparing the observed spectra

After the convolution, both are observed fluxes and can be compared directly (reflection only, \f$b = -1\f$). Per zone,
the observed ratio xillver / mirror is \f$A_\mathrm{mix} / 4\f$, i.e., the angle-mixed xillver intensity relative to
the isotropic intensity of the mirror, and **not** the angle-averaged \f$B_\mathrm{slab} \approx 0.43\f$. It depends on
the observer inclination (through the emission-angle distribution) and, as the spectral shapes differ and are shifted
in energy differently by the line profiles, on the energy band. Energy flux of the full `relxillBB` model (default
parameters):

| incl [deg]                                                         | 20   | 40   | 60   | 80   |
|--------------------------------------------------------------------|------|------|------|------|
| xillver / mirror, 0.1-1000 keV                                     | 0.32 | 0.35 | 0.39 | 0.47 |
| xillver / mirror, 0.001-2000 keV                                   | 0.44 | 0.46 | 0.48 | 0.53 |
| \f$A_\mathrm{mix}/4\f$ (rest frame, 0.1-1000 keV, `[angleweight]`) | 0.48 | 0.47 | 0.46 | 0.45 |

The ratio is therefore only a consistency check (below 1: the xillver slab re-emits less than the mirror), not a
precise number. With the mirror handed to the convolution as a flux (before 2026-10-23), the ratio was 1.3-1.9, i.e.,
the xillver reflection appeared stronger than a perfect mirror.

## Unit tests

In `test/unit/tests-returnrad-bbody-xillver.cpp`, run with `bin/tests "[rrad-bbody][xillver]"`:

- `[slabweights]` *"xillver angle weighting: semi-infinite slab flux weights 0.5*cos (incl)/n_incl and the
  angle-averaged flux"*: the weights against independently computed values, and
  \ref calc_angle_averaged_xillver_flux. Their sum equals the constant `SEMI_INFINITE_SLAB_FLUX_ISOTROPIC` = 0.25
  of the mirror (to \f$10^{-6}\f$; it is exact to \f$10^{-7}\f$, as the midpoint rule in \f$\cos i\f$ is exact
  for the integrand \f$\frac{1}{2}\mu\f$).
- `[cosnedist]` *"xillverRTR angle weighting: GR emission-angle distribution of each zone sums to 1"*.
- `[norm]` *"xillverRTR normalization: reflected flux is below the returning flux in each zone"*: energy conservation
  of the angle-averaged reflected flux (\ref calc_zone_xillver_angle_averaged_reflected_flux, as written to
  `rframe-reflect`).
- `[mirror]` *"mirrored returning radiation: full model, the xillver reflection is weaker than the mirror"*: the
  energy flux ratio xillver / mirror of the full model (\f$b = -1\f$, 20 and 80 deg) is below 1, see the table above.
  It fails for a factor 4 error in either path.
- (3d) in `test/unit/tests-returnrad-bbody-diskemission.cpp` (`bin/tests "[diskemission]"`): direct emission vs.
  mirror, see "Direct disk emission" below.
- `[angleweight]` *"xillverRTR angle weighting: angle-mixed intensity vs. angle-averaged flux per zone"*: writes the
  different angle weightings per zone to `__relxill_unit_xillverRTR-angleweight.dat`.

## Direct disk emission

The direct emission of the disk (`calc_disk_emission_observed()`) convolves the blackbody emitted per unit area with
the kernel of the line profiles, which is the same as convolving the primary of the zones `specPri` with the line
profiles of the zones (test (3c) in \ref returnrad_disk_emission). The blackbody surface of the disk emits
isotropically, so, following the rule above, its intensity \f$F_\mathrm{pri} / 0.25\f$ is convolved, as for the
mirror. This sets the normalization of the direct emission relative to the reflection.

The direct emission and the mirror are converted by the same function, `isotropic_intensity_from_flux()`
(`src/Xillspec.h`, \f$I = F / 0.25\f$), so both isotropic emitters are treated identically in the code.

**Units of the blackbody.** `bbody_photon_flux_per_keV()` only calculates the shape of the blackbody,
\f$E^2 / (e^{E/kT} - 1)\f$, without any constant; whether it is read as an intensity or as a flux only changes the
absolute normalization of all spectra by a constant factor. This does not matter for the model: the returning
radiation `specRet` is `specPri` times dimensionless photon fractions, i.e., in the same units, and the direct
emission (\f$\propto\f$ `specPri`), the mirror (\f$\propto\f$ `specRet`), and the xillver reflection
(\f$\propto N_\mathrm{xill}(\f$`specRet`\f$)\f$) all scale with the same factor. Only their ratios are determined by
the model (the absolute flux is set by the normalization parameter), and these only require that `specPri` and
`specRet` are treated as fluxes in the same way, which the test (3d) checks.

The decisive test is (3d) (see \ref returnrad_disk_emission). The direct emission and the mirror are both isotropic
emitters observed with the same line profiles, so the ratio of their observed photons in the full model must be
\f$\sum_i W_i\, n_{\mathrm{ret},i} / \sum_i W_i\, n_{\mathrm{pri},i}\f$ (\f$W_i\f$ the photons of the line profile of
zone \f$i\f$, \f$n\f$ the photons of `specRet` and `specPri`). This expected value does not contain the conversion
factor, so the test does not repeat the implementation; if only one of the two paths converted the flux into an
intensity, it would be off by a factor 4. The chain of the tests `[slabweights]` (conversion factor = xillver flux
definition), `[mirror]` (xillver vs. mirror), and (3d) (mirror vs. direct emission) fixes the normalization of the
direct emission relative to the xillver reflection.

Until 2026-10-23, the direct emission was handed to the convolution as the flux \f$F_\mathrm{pri}\f$, i.e., a factor 4
too weak relative to the xillver reflection: the reflected-to-direct ratio at \f$b = 1\f$ was a factor 4 too large
(energy flux 0.1-100 keV, default parameters: now 0.013, 0.014, 0.018, 0.028 for 20, 40, 60, 80 deg). The shape of
the direct emission (e.g., the comparison with `kerrbb`) is not affected.

## History

- 2026-10-23: the mirrored returning radiation is handed to the convolution as the intensity \f$F_\mathrm{ret} /
  0.25\f$ instead of the flux \f$F_\mathrm{ret}\f$ (a factor 4 larger), following the same rule as the xillver
  reflection. The functions returning the intensities that are convolved were renamed to
  `calc_zone_xillver_reflected_intensity()` and `calc_zone_mirrored_intensity()` (before
  `calc_zone_xillver_reflected_flux()` and `calc_zone_mirrored_flux()`).
- 2026-10-23: all diagnostic files are fluxes (see "Diagnostic output"): `rframe-reflect` is the angle-averaged flux,
  `rframe-mirror-reflect` the returning flux, and `obs-bbodyEmit` the convolution of the intensity of the primary.
- 2026-10-23: the direct emission of the disk is handed to the convolution as the intensity \f$F_\mathrm{pri} /
  0.25\f$ instead of the flux \f$F_\mathrm{pri}\f$, i.e., it is a factor 4 larger relative to the reflection
  (new unit test (3d)).
