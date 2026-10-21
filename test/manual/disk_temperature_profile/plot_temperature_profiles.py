#!/usr/bin/env python3
"""
Plot of the temperature profiles of the disk emission of the relxillBB models (see doc/returnrad_disk_emission.md):
the alpha disk with zero torque at the inner edge (default, RELXILL_RRAD_BBODY_TPROFILE=1) and the Novikov-Thorne
profile (RELXILL_RRAD_BBODY_TPROFILE=3) for several spins, both normalized to their maximal temperature (= kTbb), with
the inner edge of the disk at the ISCO.

The profiles are calculated with the same equations as in src/Relphysics.cpp (disk_Tprofile_alpha() and
disk_Tprofile_novikov_thorne()); when this script was written, both agreed to 2e-8.

Usage:  python3 plot_temperature_profiles.py [output.png]
        (default output: doc/images/disk_temperature_profiles.png in the repository)
"""
import matplotlib
import numpy as np
import os
import sys

matplotlib.use('Agg')
import matplotlib.pyplot as plt

SPINS = [0.0, 0.9, 0.998]


def r_isco(a):
    """radius of the ISCO in units of r_g (Bardeen et al. 1972)"""
    z1 = 1 + (1 - a * a) ** (1 / 3) * ((1 + a) ** (1 / 3) + (1 - a) ** (1 / 3))
    z2 = np.sqrt(3 * a * a + z1 * z1)
    return 3 + z2 - np.sign(a) * np.sqrt((3 - z1) * (3 + z1 + 2 * z2))


def tprofile_alpha(r, rin):
    """alpha disk with zero torque at rin, normalized to T = 1 at Rmax = 1.5^(4/5) rin (as disk_Tprofile_alpha())"""
    t_alpha = lambda rad: (rad / rin) ** -0.75 * (1 - np.sqrt(rin / rad)) ** 0.25
    return t_alpha(r) / t_alpha(1.5 ** 0.8 * rin)


def flux_novikov_thorne(r, rin, a):
    """radial dependence of the Novikov-Thorne flux, F ~ f/x^2 (Page & Thorne 1974, eq. 15n; Li et al. 2005, D13)"""
    x, x0 = np.sqrt(r), np.sqrt(rin)
    x1 = 2 * np.cos((np.arccos(a) - np.pi) / 3)
    x2 = 2 * np.cos((np.arccos(a) + np.pi) / 3)
    x3 = -2 * np.cos(np.arccos(a) / 3)
    integ = (x - x0 - 1.5 * a * np.log(x / x0)
             - 3 * (x1 - a) ** 2 / (x1 * (x1 - x2) * (x1 - x3)) * np.log((x - x1) / (x0 - x1))
             - 3 * (x2 - a) ** 2 / (x2 * (x2 - x1) * (x2 - x3)) * np.log((x - x2) / (x0 - x2))
             - 3 * (x3 - a) ** 2 / (x3 * (x3 - x1) * (x3 - x2)) * np.log((x - x3) / (x0 - x3)))
    return np.clip(integ / (x ** 4 * (x ** 3 - 3 * x + 2 * a)), 0, None)


def tprofile_novikov_thorne(r, rin, a):
    """Novikov-Thorne profile with zero torque at rin, normalized to its maximal temperature"""
    r_fine = np.geomspace(rin * (1 + 1e-6), 4 * rin, 200001)
    flux_max = np.max(flux_novikov_thorne(r_fine, rin, a))
    return (flux_novikov_thorne(r, rin, a) / flux_max) ** 0.25


def plot(fname):
    surface, ink, ink2, grid = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
    colors = {0.0: '#2a78d6', 0.9: '#eb6834', 0.998: '#1baf7a'}
    plt.rcParams.update({'font.size': 10, 'axes.edgecolor': ink2, 'axes.labelcolor': ink, 'xtick.color': ink2,
                         'ytick.color': ink2, 'text.color': ink, 'axes.linewidth': 0.8})

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 7.2), sharex=True,
                                   gridspec_kw={'height_ratios': [3, 2], 'hspace': 0.08})
    fig.patch.set_facecolor(surface)
    for ax in (ax1, ax2):
        ax.set_facecolor(surface)
        ax.grid(True, which='major', color=grid, lw=0.6)
        ax.set_axisbelow(True)
        for side in ['top', 'right']:
            ax.spines[side].set_visible(False)

    x = np.geomspace(1 + 1e-6, 1000, 3000)  # r / Rin
    t_alpha = tprofile_alpha(x, 1.0)  # independent of the spin in units of r / Rin
    ax1.plot(x, t_alpha, color=ink2, lw=1.8, ls=(0, (5, 3)), label='alpha disk (Newtonian, zero torque), any spin')
    ax1.plot(x[np.argmax(t_alpha)], 1.0, 'o', ms=6, color=ink2, mec=surface, mew=1.5, zorder=5)

    for a in SPINS:
        rin = r_isco(a)
        t_nt = tprofile_novikov_thorne(x * rin, rin, a)
        ax1.plot(x, t_nt, color=colors[a], lw=1.8, label=f'Novikov-Thorne, a = {a:g}')
        ax1.plot(x[np.argmax(t_nt)], 1.0, 'o', ms=6, color=colors[a], mec=surface, mew=1.5, zorder=5)
        mask = x > 1.02
        ratio = t_nt[mask] / t_alpha[mask]
        ax2.plot(x[mask], ratio, color=colors[a], lw=1.8)
        ax2.annotate(f'a = {a:g}', xy=(x[mask][-1], ratio[-1]), xytext=(6, 0), textcoords='offset points',
                     va='center', fontsize=9, color=ink)

    ax1.annotate('maximum: 1.36 Rin (alpha)\n1.59 / 1.48 / 1.28 Rin (NT, a = 0 / 0.9 / 0.998)', xy=(1.45, 1.0),
                 xytext=(4.0, 0.93), fontsize=9, color=ink2, arrowprops=dict(arrowstyle='-', color=ink2, lw=0.8))
    ax1.set_xscale('log')
    ax1.set_xlim(1, 1000)
    ax1.set_ylim(0, 1.08)
    ax1.set_ylabel('T / T$_\\mathrm{max}$  (T$_\\mathrm{max}$ = kT$_\\mathrm{bb}$)')
    ax1.legend(frameon=False, loc='upper right', bbox_to_anchor=(1.0, 0.82), fontsize=9)
    ax1.set_title('Disk temperature profiles, inner edge at the ISCO', fontsize=10, color=ink, loc='left')

    ax2.axhline(1, color=ink2, lw=0.8)
    ax2.set_ylim(0.82, 1.28)
    ax2.set_ylabel('T$_\\mathrm{NT}$ / T$_\\mathrm{alpha}$')
    ax2.set_xlabel('r / R$_\\mathrm{in}$')
    ax2.text(1.12, 1.235, 'Novikov-Thorne relative to the alpha disk (same T$_\\mathrm{max}$)', fontsize=9, color=ink2)

    fig.subplots_adjust(left=0.11, right=0.89, top=0.94, bottom=0.08)
    fig.savefig(fname, dpi=150, facecolor=surface)
    print(f" plot written to {fname}")


if __name__ == '__main__':
    repo_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
    default = os.path.join(repo_dir, 'doc', 'images', 'disk_temperature_profiles.png')
    output = sys.argv[1] if len(sys.argv) > 1 else default
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    plot(output)
