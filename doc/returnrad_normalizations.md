# Thermal Returning Radiation: Flux Normalizations {#returnrad_normalizations}

This page documents the normalization factors applied inside
\ref calc_zone_direct_plus_xillver_reflected_flux in
`src/Relreturn_BlackBody.cpp`, which builds the per-zone emitted and reflected spectrum for
`relxillBB` before it is handed to the relativistic convolution. The important part is that the correct normalization is
applied to the reflected spectrum, following the normalization criterion for the direct **returning** spectrum of
xillver and the fraction of returning flux (see Dauser et al. 2026).

## Context

\ref relxillbb_kernel loops over the radial zones of the return-radiation grid.
For each zone `izone` it combines two contributions in the **disk frame**:

1. The **direct** primary thermal emission of the disk surface,
   `returnSpec->specPri[izone]`.
2. The **reflected** compenent (by the xillver table) for the given returning thermal
   radiation, whose illuminating spectrum is `returnSpec->specRet[izone]`.

The result is written into `xill_flux_returnrad` and only afterwards folded with
the relativistic line profile of that zone
(`convolveSpectrumFFTNormalized` with `relProfile->flux[izone]`) and rebinned
onto the input grid. Nothing on this page applies the relativistic broadening. That is deliberately the convolution's
job, and it matters for the intensity-to-flux question discussed below.

## Call tree

Only the nodes that carry or consume a normalization are shown. Leaf utilities
(table indexing, interpolation, FITS row lookup) are omitted; use the
doxygen-generated call graph on the function page for the exhaustive version.

\dot
digraph returnrad_norm {
rankdir=LR;
node [shape=box, fontname="Helvetica", fontsize=10, color="#5b6b7b", fontcolor="#1b1b1b", style="filled", fillcolor="#f2f5f8"];
edge [fontname="Helvetica", fontsize=9, color="#5b6b7b", fontcolor="#1b1b1b"];

entry  [label="calc_zone_direct_plus_xillver_reflected_flux", fillcolor="#dbe7f3"];

angw   [label="get_angle_weighted_xill_spec\n(Xillspec.cpp)"];
spec   [label="get_xillver_spectra\ntable lookup + cache"];
angdep [label="calc_xillver_angdep\n-> sum_xillver_spectra_over_incl"];
rebin1 [label="_rebin_spectrum\nto returnSpec->ener"];

nfac   [label="get_xillver_norm_factor_from_primary_spectrum\n(static, Relreturn_BlackBody.cpp:566)", fillcolor="#fbeedd"];
rebin2 [label="_rebin_spectrum\nto coarse xillver grid"];
cnorm  [label="calc_xillver_norm_from_primary_spectrum\n(Xillspec.cpp:208)", fillcolor="#fbeedd"];

comb   [label="per-bin combine loop", fillcolor="#dbe7f3"];

entry -> angw   [label="  I_xill(E)"];
entry -> nfac   [label="  specRet[izone]"];
angw -> spec;
angw -> angdep [label="  rel_cosne->dist[izone]"];
angw -> rebin1;
nfac -> rebin2;
nfac -> cnorm  [label="  N_xill"];
angw -> comb;
nfac -> comb;
comb -> comb   [label="  x |boost| ; + specPri if boost >= 0"];
}
\enddot

## The combination

\f[
F_\mathrm{zone}(E) \;=\;
|b|\; N_\mathrm{xill}\; I_\mathrm{xill}(E)
\;+\;
\Theta(b)\, F_\mathrm{pri}(E)
\f]

with \f$b\f$ = `xill_param.boost`, \f$\Theta\f$ the Heaviside step
(\f$\Theta = 1\f$ for \f$b \ge 0\f$, else 0),
\f$I_\mathrm{xill}\f$ the angle-weighted xillver spectrum and
\f$F_\mathrm{pri}\f$ = `returnSpec->specPri[izone]`.

**The direct term carries no normalization factor of its own.** It is added
exactly as `spec_returnrad_blackbody()` produced it. There is no \f$N_\mathrm{dir}\f$;
only the reflected term is scaled.

## The factors

### `N_xill` — illumination relative to the xillver table

\ref calc_xillver_norm_from_primary_spectrum integrates the *illuminating*
spectrum over the xillver normalization band
(`EMIN_XILLVER_NORMALIZATION` = 0.1 keV to `EMAX_XILLVER_NORMALIZATION` = 1000 keV,
`Xillspec.h:35`) and divides by the reference illumination the tables were
computed for:

\f[
N_\mathrm{xill} \;=\;
\frac{\sum_i F_i \,\bar{E}_i \, 10^{20}\, c_\mathrm{keV \to erg}}{10^{15}/4\pi}
\f]

\f$c_\mathrm{keV \to erg} = 1.602177 \times 10^{-9}\f$ and
\f$\bar{E}_i = \tfrac{1}{2}(E_i + E_{i+1})\f$. The denominator is the
\f$n = 10^{15}\,\mathrm{cm^{-3}}\f$ reference normalization of the xillver
tables (Appendix of Dauser et al. 2016). \f$N_\mathrm{xill}\f$ is therefore
dimensionless: *actual returning illumination / table reference illumination*.

Two details that are easy to get wrong:

- The sum runs over the **coarse global xillver grid**
  (`get_coarse_xillver_energrid()`), not the model grid (for a faster model evaluation, as the normalization does not
  need the high precision of the regular xillver grid). That is why
  \ref get_xillver_norm_factor_from_primary_spectrum rebins first.
- The argument is the returning spectrum `specRet[izone]` (what *hits* the disk), never the primary emitted spectrum
  `specPri`  (what *leaves* it).

### The angle weighting is *not* a flux normalization

\ref calc_xillver_angdep weights the table inclinations with
`rel_profile->rel_cosne->dist[izone]` (normalized to sum = 1) and returns the
**angle-mixed intrinsic intensity**. It deliberately does *not* apply the
semi-infinite-slab conversion \f$F(i) = \tfrac{1}{2}\cos i \; I(i)\f$
(`norm_factor_semi_infinite_slab()`), because for relxill-type models the
intensity-to-flux conversion is carried by the GR convolution in the kernel.

Applying the slab factor here as well would double-count it. The separate path
\ref get_semi_infinite_slab_flux_weights exists only for the angle-averaging the
non-relativistic case and must not be mixed into this one (see comment in Garcia+14, Eq. 19).

### `boost`

`|boost|` scales the reflected component only. A **negative** boost selects
return-radiation-only output by suppressing the direct term (identical to the
convention used for reflection fraction / boost elsewhere in relxill.

## Why the split

The angle-weighted xillver spectrum and the normalization are kept apart
because they vary on different parameter axes:

- \f$I_\mathrm{xill}\f$ depends on the xillver parameters (\f$\xi\f$, density,
  \f$kT_\mathrm{bb}\f$) and on the emission-angle distribution, and comes
  from a cached table read.
- \f$N_\mathrm{xill}\f$ depends only on the returning illumination of the zone,
  which is recomputed whenever the geometry changes.

Keeping them separate lets the table cache stay valid across changes to the
returning-radiation geometry. It also mirrors
`Relreturn_Corona.cpp`, so the two return-radiation flavours stay comparable.

When \ref use_iongrad_relxillbb is active, `set_xill_param_from_iongrad()`
overwrites the zone's ionization *before* this function runs, so the
ionization gradient enters through \f$I_\mathrm{xill}\f$, not through
\f$N_\mathrm{xill}\f$.

## Not part of the model calculation

`Relreturn_BlackBody.cpp` contains more normalizations than are actually
applied to the model spectra. These run **only** under `shouldOutfilesBeWritten()`
and feed diagnostic output:

- \ref get_xillver_primary_bbody and
  \ref get_xillver_primary_bbody_normalized_at_highener
- \ref normalization_factor_bbody_at_highener and
  \ref calc_normfac_bbody_at_highenergy

A high-energy matching factor of this kind *was* once applied in the real path
and was removed in commit 6f07241 as incorrect. Do not reintroduce it here.

\ref calc_zone_direct_plus_mirrored_flux is likewise not a physics path: it is
a test shortcut selected by \ref mirrored_reflection_instead_of_xillver that
mirrors `specRet` instead of reflecting it, and it applies `|boost|` and the
direct term but no `N_xill`.

## Pitfalls

- Omitting \f$N_\mathrm{xill}\f$ does not produce a correct spectral shape, even not for only reflection (boost<0). The
  reason for this is that reflect spectrum of relxillBB is a sum contribution of different radii. And the combined
  reflection spectrum is only correct, if the spectra of these single radii are correctly normalized.
