from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm


def parse_header(path: Path) -> dict:
    meta = {}
    with path.open() as f:
        for line in f:
            if not line.startswith("#"):
                break
            for tok in line.lstrip("#").strip().split():
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    meta[k] = v
    return meta


def load_grid(path: Path):
    meta = parse_header(path)
    na = int(meta["na"])
    ne = int(meta["ne"])
    amin = float(meta["amin"])
    amax = float(meta["amax"])
    emin = float(meta["emin"])
    emax = float(meta["emax"])
    dt_days = float(meta["dt_days"])
    t_years = float(meta["t_years"])
    integrator = meta.get("integrator", "?")

    # Skip both '#'-prefixed metadata lines AND the column-name row, then
    # parse the rest as raw floats (NaN tokens are honored by float()).
    rows = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("a_init"):
                continue
            rows.append([float(x) for x in line.split(",")])
    arr = np.asarray(rows, dtype=float)
    if arr.shape[0] != na * ne or arr.shape[1] != 3:
        raise ValueError(f"unexpected CSV shape {arr.shape} for na*ne={na*ne}")
    da = arr[:, 2].reshape(na, ne)

    a_axis = np.linspace(amin, amax, na)
    e_axis = np.linspace(emin, emax, ne)
    return {
        "a": a_axis,
        "e": e_axis,
        "da": da,
        "dt_days": dt_days,
        "t_years": t_years,
        "integrator": integrator,
    }


def tf_contour_eccentricities(a_axis, dt_years):
    """For each a in a_axis return the e at which Tf(a, e) = dt_years.

    Tf = T_orbit * (1-e)^2 * sqrt(1-e^2), with T_orbit = 2*pi*a^(3/2)
    (so T_orbit in years = a^(3/2) when G=Msun=1, year=2*pi).

    Solve  a^{3/2} * (1-e)^2 * sqrt(1-e^2) = dt_years  for e in (0, 1).
    """
    e = np.empty_like(a_axis)
    ee = np.linspace(0.0, 0.999999, 10000)
    fac = (1 - ee) ** 2 * np.sqrt(1 - ee ** 2)
    for i, a in enumerate(a_axis):
        target = dt_years / (a ** 1.5)
        if target >= 1.0:
            e[i] = 0.0
            continue
        if target <= 0.0:
            e[i] = np.nan
            continue
        idx = np.searchsorted(-fac, -target)
        if idx == 0 or idx >= len(ee):
            e[i] = np.nan
            continue
        # linear interp between ee[idx-1] and ee[idx]
        f0, f1 = fac[idx - 1], fac[idx]
        e[i] = ee[idx - 1] + (target - f0) / (f1 - f0) * (ee[idx] - ee[idx - 1])
    return e


def make_figure(grid, save_path: Path, *, vmin=1e-16, vmax=1e-6,
                title_suffix=None):
    a_axis = grid["a"]
    e_axis = grid["e"]
    da = grid["da"]
    dt_days = grid["dt_days"]
    t_years = grid["t_years"]

    # pcolormesh expects coords of size (ne+1, na+1) for edges.
    da_amin = grid["a"][0] - 0.5 * (grid["a"][1] - grid["a"][0])
    da_amax = grid["a"][-1] + 0.5 * (grid["a"][1] - grid["a"][0])
    da_emin = grid["e"][0] - 0.5 * (grid["e"][1] - grid["e"][0])
    da_emax = grid["e"][-1] + 0.5 * (grid["e"][1] - grid["e"][0])

    A, E = np.meshgrid(a_axis, e_axis, indexing="ij")

    # Tf contour curves
    dt_years = dt_days / 365.25
    a_fine = np.linspace(a_axis.min(), a_axis.max(), 800)
    e_tf1 = tf_contour_eccentricities(a_fine, dt_years)
    e_tf2 = tf_contour_eccentricities(a_fine, 2.0 * dt_years)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    extent = [da_amin, da_amax, da_emin, da_emax]

    # Left panel
    mag = np.abs(da)
    # Floor zeros so they don't break LogNorm.
    mag_disp = np.where(mag <= 0, vmin * 0.1, mag)
    mag_disp = np.clip(mag_disp, vmin, vmax)
    norm = LogNorm(vmin=vmin, vmax=vmax)
    cmap = mpl.colormaps.get_cmap("inferno").copy()
    cmap.set_bad("white")
    nan_mask = ~np.isfinite(da)
    mag_masked = np.ma.masked_where(nan_mask, mag_disp)
    im0 = axes[0].pcolormesh(A.T, E.T, mag_masked.T, cmap=cmap, norm=norm,
                              shading="nearest", rasterized=True)
    cbar0 = fig.colorbar(im0, ax=axes[0], pad=0.02)
    cbar0.set_label(r"$|\Delta a / a|$")
    axes[0].set_xlabel("semi-major axis [AU]")
    axes[0].set_ylabel("eccentricity")

    # Right panel
    sign = np.zeros_like(da)
    sign[da > 0] =  1
    sign[da < 0] = -1
    sign_masked = np.ma.masked_where(nan_mask, sign)
    sign_cmap = ListedColormap(["#2754d4", "#ffffff", "#e6262e"])
    sign_cmap.set_bad("white")
    bnorm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], sign_cmap.N)
    im1 = axes[1].pcolormesh(A.T, E.T, sign_masked.T, cmap=sign_cmap,
                              norm=bnorm, shading="nearest", rasterized=True)
    cbar1 = fig.colorbar(im1, ax=axes[1], ticks=[-1, 0, 1], pad=0.02)
    cbar1.set_label(r"$\mathrm{sign}(\Delta a)$")
    axes[1].set_xlabel("semi-major axis [AU]")
    axes[1].set_ylabel("eccentricity")

    # NaN hatching overlay
    if nan_mask.any():
        from matplotlib.patches import PathPatch
        from matplotlib.path import Path as MplPath
        mpl.rcParams['hatch.linewidth'] = 0.5
        for ax in axes:
            ax.contourf(A.T, E.T, nan_mask.T.astype(float), levels=[0.5, 1.5],
                        colors="none", hatches=["xx"], zorder=4)

    # Tf contours + Mercury/Venus markers on both panels.
    for ax in axes:
        ax.plot(a_fine, e_tf1, ls="--", color="#1f5be0", lw=1.6, zorder=5,
                 label=r"$T_f = \Delta t$")
        ax.plot(a_fine, e_tf2, ls="--", color="#d62a2a", lw=1.6, zorder=5,
                 label=r"$T_f = 2\,\Delta t$")
        # Mercury (a=0.387 AU, e=0.206) and Venus (a=0.723 AU, e=0.007),
        # drawn as a small white tile labelled
        for (ax_xy, sym) in (((0.387, 0.206), "\u263F"),
                             ((0.723, 0.007), "\u2640")):
            ax.scatter([ax_xy[0]], [ax_xy[1]], s=260, marker="s",
                        facecolor="white", edgecolor="black", lw=0.7,
                        zorder=6)
            ax.annotate(sym, ax_xy, ha="center", va="center", fontsize=11,
                         zorder=7)

        ax.set_xlim(a_axis.min(), a_axis.max())
        ax.set_ylim(e_axis.min(), e_axis.max())

    title = (title_suffix if title_suffix
             else f"WHFast Kepler-step stability ({grid['integrator']})")
    fig.suptitle(
        f"{title}  —  $\\Delta t = {dt_days:g}$ d, $t = {t_years:g}$ yr",
        fontsize=11)

    fig.savefig(save_path, dpi=180, bbox_inches="tight")
    print(f"  wrote {save_path}", file=sys.stderr)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv", type=Path, required=True,
                    help="Sweep CSV produced by ./rebound")
    p.add_argument("--out", type=Path, required=True,
                    help="Output PNG path")
    p.add_argument("--vmin", type=float, default=1e-16)
    p.add_argument("--vmax", type=float, default=1e-6)
    p.add_argument("--title-suffix", type=str, default=None)
    args = p.parse_args()

    grid = load_grid(args.csv)
    print(f"  loaded {args.csv} : {grid['da'].shape}, "
          f"integrator={grid['integrator']}, "
          f"dt={grid['dt_days']}d, t={grid['t_years']}yr", file=sys.stderr)
    make_figure(grid, args.out, vmin=args.vmin, vmax=args.vmax,
                title_suffix=args.title_suffix)


if __name__ == "__main__":
    main()
