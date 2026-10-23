# Thermal Returning Radiation: Flux Normalizations {#returnrad_normalizations}

This page documents the normalization factors applied inside
\ref calc_zone_xillver_reflected_intensity in
`src/Relreturn_BlackBody.cpp`, which builds the per-zone reflected spectrum for
`relxillBB` before it is handed to the relativistic convolution. The important part is that the correct normalization is
applied to the reflected spectrum, following the normalization criterion for the direct **returning** spectrum of
xillver and the fraction of returning flux (see Dauser et al. 2016).

## Context

\ref relxillbb_kernel loops over the radial zones of the return-radiation grid.
For each zone `izone` it calculates the **reflected** component (by the xillver table) in the **disk frame**, for the
given returning thermal radiation, whose illuminating spectrum is `returnSpec->specRet[izone]`. The returning
radiation is calculated from the primary thermal emission of the zones, `returnSpec->specPri[izone]`.

The result (an intensity, see \ref xillver_angular_dependence) is written into `xill_intensity` and only afterwards
folded with
the relativistic line profile of that zone
(`convolveSpectrumFFTNormalized` with `relProfile->flux[izone]`) and rebinned
onto the input grid. Nothing on this page applies the relativistic broadening. That is deliberately the convolution's
job, and it matters for the intensity-to-flux question discussed below.

The **direct** thermal emission of the disk is not part of the zones. It is calculated for the whole disk at once
(`calc_disk_emission_observed()`), from the same radial points and temperatures as the primary of the zones, see
\ref returnrad_disk_emission.

## Call tree

Only the nodes that carry or consume a normalization are shown. Leaf utilities
(table indexing, interpolation, FITS row lookup) are omitted; use the
doxygen-generated call graph on the function page for the exhaustive version.

\dot
digraph returnrad_norm {
rankdir=LR;
node [shape=box, fontname="Helvetica", fontsize=10, color="#5b6b7b", fontcolor="#1b1b1b", style="filled", fillcolor="#f2f5f8"];
edge [fontname="Helvetica", fontsize=9, color="#5b6b7b", fontcolor="#1b1b1b"];

entry  [label="calc_zone_xillver_reflected_intensity", fillcolor="#dbe7f3"];

angw   [label="get_angle_weighted_xill_spec\n(Xillspec.cpp)"];
spec   [label="get_xillver_spectra\ntable lookup + cache"];
angdep [label="calc_xillver_angdep\n-> sum_xillver_spectra_over_incl"];
rebin1 [label="_rebin_spectrum\nto returnSpec->ener"];

nfac   [label="get_xillver_norm_factor_from_primary_spectrum\n(static, Relreturn_BlackBody.cpp:566)", fillcolor="#fbeedd"];
rebin2 [label="_rebin_spectrum\nto coarse xillver grid"];
cnorm  [label="calc_xillver_norm_from_primary_spectrum\n(Xillspec.cpp:208)", fillcolor="#fbeedd"];

comb   [label="per-bin combine loop", fillcolor="#dbe7f3"];

entry -> angw   [label="  F_xill(E)"];
entry -> nfac   [label="  specRet[izone]"];
angw -> spec;
angw -> angdep [label="  rel_cosne->dist[izone]"];
angw -> rebin1;
nfac -> rebin2;
nfac -> cnorm  [label="  N_xill"];
angw -> comb;
nfac -> comb;
comb -> comb   [label="  x |boost|"];
}
\enddot

## Calculating the Normalized Reflection Spectrum per Zone

Two separate things set the reflected flux, and keeping them apart is what
makes the scheme correct.

### The table fixes the relative normalization between the primary and the reflected spectrum

The xillver tables are computed and normalized according to the integral over
the primary (illuminating) spectrum defined in the Appendix of
Dauser et al. (2016). **This relative normalization between primary and reflected spectrum
is a property of the table and must not be altered by any other means.**

The normalization is defined as the integral over the *illuminating* energy flux
over the xillver normalization band (`EMIN_XILLVER_NORMALIZATION` = 0.1 keV to `EMAX_XILLVER_NORMALIZATION` = 1000 keV,
`Xillspec.h:35`) and defined to match \f$ \int_{0.1\,\mathrm{keV}}^{1000\,\mathrm{keV}} F_\mathrm{illum}(E)\, dE
= 10^{-20}\, \frac{n \xi}{4 \pi} = 10^{-20}\, \frac{10^{15}}{4 \pi} \f$. The last part follows as internally in
relxill the xillver spectra are re-normalized to a density of \f$n=10^{15}\f$ and ionization \f$\xi=1\f$.

### The returning radiation sets the absolute normalization of the reflected spectrum

To calculate the absolute flux of the reflected spectrum, the xillver spectra (from the table)
must be scaled by the returning radiation flux (`returnSpec->specRet[izone]`).

For this purpose we define the **xillver normalization factor** as the relative normalization of the given illuminating
(returning) spectrum \f$F_\mathrm{ret}(E) \f$ to the intrinsic xillver spectrum normalization  (see above):

\f[
N_\mathrm{xill}\bigl(F_\mathrm{ret}\bigr)\; \;=\;
\frac{\sum_i F_{\mathrm{ret},i} \,\bar{E}_i \, 10^{20}\, c_\mathrm{keV \to erg}}{10^{15}/4\pi}
\f]

\f$c_\mathrm{keV \to erg} = 1.602177 \times 10^{-9}\f$ and
\f$\bar{E}_i = \frac{1}{2}(E_i + E_{i+1})\f$. \f$N_\mathrm{xill}\f$ is therefore
dimensionless: *actual returning illumination / table reference illumination*.

Therefore, the reflected spectrum per zone, before the boost \f$|b|\f$ is applied, is given by:

\f[
F_\mathrm{refl}(E) \;=\; N_\mathrm{xill}\bigl(F_\mathrm{ret}\bigr)\; F_\mathrm{xill}(E)
\f]

### Combined spectrum

The observed spectrum is the sum of the direct emission of the whole disk and of the reflected spectra of the zones,
each convolved with the relativistic line profile \f$P_i\f$ of its zone \f$i\f$:

\f[
F_\mathrm{obs}(E) \;=\; \Theta(b)\, F_\mathrm{dir,obs}(E) \;+\; |b| \sum_i P_i \otimes \left[ N_\mathrm{xill}\bigl(F_{\mathrm{ret},i}\bigr)\; F_{\mathrm{xill},i} \right](E)
\f]

where \f$F_\mathrm{dir,obs}(E)\f$ is the observed direct emission of the disk (`calc_disk_emission_observed()`) and
\f$F_{\mathrm{ret},i}(E) = \f$ `specRet[i]`, where the latter only sets the correct absolute normalization of the
reflected xillver spectrum \f$F_{\mathrm{xill},i}(E)\f$. Furthermore, \f$b\f$ = `xill_param.boost`, and \f$\Theta\f$
is the Heaviside step function (\f$\Theta(b) = 1\f$ for \f$b \ge 0\f$, else 0), i.e. the direct term is only added
for a non-negative boost, while the reflected term is always scaled by \f$|b|\f$.

**The direct term carries no normalization factor of its own.** It is an absolute flux, calculated from the same
emitted primary as the returning radiation: `specPri[i]`, from which `specRet[i]` is calculated, is the emission of
the same radial points with the same temperatures as the direct emission, binned into the zones. There is no
\f$N_\mathrm{dir}\f$; only the reflected term is scaled. The consistency of the direct and the reflected term is tested
by the unit tests (3) described in \ref returnrad_disk_emission. As the reflected spectra, the direct emission is
handed to the convolution as an intensity (\f$F_\mathrm{pri} / 0.25\f$ of the isotropically emitting disk surface),
see \ref xillver_angular_dependence; this relative normalization is tested by (3d).

### What follows from this

- \f$N_\mathrm{xill}\f$ is **not** a tunable factor and **not** a correction.
  It is the absolute flux of the returning primary, expressed in table units.
  Its value is whatever the returning radiation is; it is never chosen to make
  a spectrum look right.
- **Nothing else may change the reflection-to-primary ratio.** Matching the
  continuum in some band, re-anchoring to a differently-weighted
  normalization, or forcing agreement with an assumed primary shape all
  override the table's own calibration. These are the operations that are
  forbidden.
- \f$|\mathrm{boost}|\f$ is not an exception to that. It is an explicit
  physical model parameter -- the reflection fraction the user sets -- applied
  on top of a correctly normalized spectrum, and deliberately *does* change
  the reflected-to-direct ratio. That is its purpose, and it is a model
  parameter rather than a normalization.

If the tabulated reflection does not suit the illumination at hand -- most
often because the incident *shape* differs from what the table assumed -- the
remedy is a table computed for the correct incident spectrum, never a
rescaling of the existing one. The `radindex` dimension of
`xillverRTR_v0.3.fits` exists for exactly this reason: it tabulates the
per-radial-zone returning-radiation illumination instead of correcting for it
after the fact.

Both normalization factors that were removed from this path in 2026 changed
the relative normalization; see \ref returnrad_norm_removed "Pitfalls (2)".

### The incident spectrum of the table and the returning spectrum

The reflection spectra of the xillver table are calculated for a given incident spectrum, in the following called
the *xillver primary*. Physically, it is the returning spectrum `specRet[izone]`. The xillver primary is, however, a
property of the table and is not calculated in the model, while `specRet[izone]` is calculated by the model for the
given parameters:

- `relxillBB` (`xillverRTR` table): the xillver primary of each radial zone (`radindex`) was taken from
  `specRet[izone]` when the table was calculated, so both are identical by design. The current table
  (`xillverRTR_v0.3.fits`) is a mock version for a single parameter set (a = 0.998, \f$kT_\mathrm{bb}\f$ = 2 keV,
  logxi = 2, logN = 18), which the model fixes until the final table is available. The final table will in addition
  depend on \f$kT_\mathrm{bb}\f$, logxi, logN, the spin (which defines the inner edge and therefore the zones), and
  later A_Fe, with the xillver primary of each zone calculated for
  each combination of \f$kT_\mathrm{bb}\f$ and spin, such that both spectra stay identical.
- `relxillBB_NS` (`xillverNS` table, kept in the code for development): the xillver primary is a blackbody of
  temperature \f$kT_\mathrm{bb}\f$, the same for all zones. For computational reasons, it is only an approximation of
  the returning spectrum.

The model does not assume that both spectra are identical, as this may change in the future (e.g., for a table
calculated for a different or more general incident spectrum). Only the integrated energy flux of `specRet[izone]`
(0.1-1000 keV) enters the model, through \f$N_\mathrm{xill}\f$, while the shape of the xillver primary is contained in
the table. A difference in the spectral shape would therefore not be corrected (see above), but would need to be kept
in mind when comparing the reflected spectrum with the returning spectrum.

## The factors

### calculating N_xill: two things to get right

When calling \ref calc_xillver_norm_from_primary_spectrum two details are easy to get wrong:

- The sum runs over the **coarse global xillver grid**
  (`get_coarse_xillver_energrid()`), not the model grid (for a faster model evaluation, as the normalization does not
  need the high precision of the regular xillver grid). That is why
  \ref get_xillver_norm_factor_from_primary_spectrum rebins first.
- The argument is the returning spectrum `specRet[izone]` (what *hits* the disk), never the primary emitted spectrum
  `specPri`  (what *leaves* it).

### The angle weighting is not a flux normalization

\ref calc_xillver_angdep weights the table inclinations with
`rel_profile->rel_cosne->dist[izone]` (normalized to sum = 1) and returns the **angle-mixed intrinsic intensity**
\f$F_\mathrm{xill}\f$ (despite its name, an intensity). It deliberately does *not*
apply the semi-infinite-slab conversion \f$\frac{1}{2}\cos i\f$ from intensity to flux, because the relativistic
convolution with the line profile carries this conversion. The same is done in standard relxill (\ref relxill_kernel,
step 5).

The conversion between intensity and flux, where it is applied, and the rule that every spectrum handed to the
relativistic convolution must be an intensity (which also sets the normalization of the mirrored returning radiation)
are described in \ref xillver_angular_dependence. In particular, the reflected spectrum of a zone can therefore not be
compared directly with the returning flux \f$F_\mathrm{ret}\f$; only its angle-averaged flux with the slab weights can
(test `[xillverRTR][norm]`).

### boost

`|boost|` scales the reflected component only. A **negative** boost selects
return-radiation-only output by suppressing the direct term (identical to the
convention used for the reflection fraction / boost elsewhere in relxill). For `boost=0`, only the direct term is
calculated (the zones are skipped).

## Why the split

The angle-weighted xillver spectrum and the normalization are kept apart
because they vary on different parameter axes:

- \f$F_\mathrm{xill}\f$ depends on the xillver parameters (\f$\xi\f$, density,
  \f$kT_\mathrm{bb}\f$) and on the emission-angle distribution, and comes
  from a cached table read.
- \f$N_\mathrm{xill}\f$ depends only on the returning illumination of the zone,
  which is recomputed whenever the geometry changes.

Keeping them separate lets the table cache stay valid across changes to the
returning-radiation geometry. It also mirrors
`Relreturn_Corona.cpp`, so the two return-radiation flavours stay comparable.

When \ref use_iongrad_relxillbb is active, `set_xill_param_from_iongrad()`
overwrites the zone's ionization *before* this function runs, so the
ionization gradient enters through \f$F_\mathrm{xill}\f$, not through
\f$N_\mathrm{xill}\f$.

## Not part of the model calculation

\ref calc_zone_mirrored_intensity is not a physics path: it is
a test shortcut selected by \ref mirrored_reflection_instead_of_xillver that
mirrors `specRet` instead of reflecting it (the direct term is added for the whole disk, as for the reflection). It
applies `|boost|` but no `N_xill`, as the returning flux is used directly. Following the same rule as the xillver
reflection, the mirror keeps the spectral shape and flux of the returning radiation (which already contains the
angles of incidence), re-emits it with an isotropic intensity, and is handed to the convolution as the intensity
\f$I_\mathrm{mirror} = F_\mathrm{ret} / 0.25\f$, see \ref xillver_angular_dependence (until 2026-10-23, the flux
\f$F_\mathrm{ret}\f$ itself, i.e., a factor 4 too weak compared to the xillver reflection).

## Pitfalls

(1) Omitting \f$N_\mathrm{xill}\f$ does not produce a correct spectral shape, even not for only reflection (boost<0).
The
  reason for this is that reflect spectrum of relxillBB is a sum contribution of different radii. And the combined
  reflection spectrum is only correct, if the spectra of these single radii are correctly normalized.

(2) \anchor returnrad_norm_removed
Two factors were applied in this path and have been removed. Both changed the
**relative** normalization -- the reflection-to-primary ratio that the table
alone defines -- which is why both were wrong. Neither was about the absolute
scale; that was already correctly set by \f$N_\mathrm{xill}\f$:

- **`normfacMatchAtHighEnergies`** (removed in `6f07241`) forced the reflected
  continuum to match a blackbody built from the returning spectrum above
  \f$4\,kT_\mathrm{bb}\f$. Matching a band overrides the table's own
  reflection-to-primary calibration with an unrelated criterion. It was a fudge
  parameter trying to make xillverNS work as returning radiation model, but should
  not be used for the final model and proper xillverRTR table.
- **`xill_norm_avg / xill_norm`** (removed in `df34569`) re-anchored the angle-weighted spectrum to
  the angle-*averaged* table normalization. That is a second, different
  normalization criterion layered on top of the first.
