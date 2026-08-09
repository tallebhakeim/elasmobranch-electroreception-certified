# -*- coding: utf-8 -*-
"""Figures for the sensitivity and robustness section (JTB resubmission).

Reads sensitivity_results.json (written by sensitivity_lorenzini.py) and writes four
figures into figures/. Every panel shows a GUARANTEED interval, never a point value:
the vertical bars are the brackets themselves, and the dashed line at 1 is the claim
under test (the gel canal raises the access conductance).

Run:  python3 make_sensitivity_figures.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
RES = json.load(open(os.path.join(HERE, "sensitivity_results.json")))

OK, NO = "#1b6ca8", "#c1121f"            # disjoint / not disjoint
plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "figure.dpi": 200,
                     "savefig.bbox": "tight", "axes.spines.top": False,
                     "axes.spines.right": False})


def plain_ticks(ax, x, fmt="{:g}"):
    """Log axes default to 4x10^-1 style labels, which is unreadable for millimetres."""
    from matplotlib.ticker import NullFormatter
    ax.set_xticks(x)
    ax.set_xticklabels([fmt.format(v) for v in x], fontsize=7.5)
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="x", which="minor", length=0)


def pad(ax, x, logx, frac=0.16):
    x = np.asarray(x, float)
    if logx:
        r = (x.max() / x.min()) ** frac
        ax.set_xlim(x.min() / r, x.max() * r)
    else:
        d = (x.max() - x.min()) * frac or 0.1 * abs(x.max())
        ax.set_xlim(x.min() - d, x.max() + d)


def bars(ax, x, rows, logx=False, width=None):
    """Draw each guaranteed gain bracket as a vertical bar at its abscissa."""
    x = np.asarray(x, float)
    if width is None:
        width = 0.035 * (np.log10(x.max() / x.min()) if logx else (x.max() - x.min()) or 1)
    for xi, r in zip(x, rows):
        lo, hi = r["gain"]
        c = OK if r["disjoint"] else NO
        w = xi * (10 ** width - 1) if logx else width
        ax.add_patch(plt.Rectangle((xi - w / 2, lo), w, hi - lo, facecolor=c,
                                   edgecolor=c, alpha=0.55, lw=0.8))
        ax.plot([xi - w / 2, xi + w / 2], [lo, lo], color=c, lw=1.2)
        ax.plot([xi - w / 2, xi + w / 2], [hi, hi], color=c, lw=1.2)
    ax.axhline(1.0, ls="--", lw=1.0, color="0.35")
    if logx:
        ax.set_xscale("log")
    ax.set_ylabel("guaranteed gain $G_{canal}/G_{no\\ canal}$")


def fig_mesh():
    rows = RES["mesh"]
    h = [r["h"] * 1e3 for r in rows]
    ntri = [r["ntri"] for r in rows]
    fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.0), layout="constrained")
    bars(ax[0], h, rows, logx=True)
    pad(ax[0], h, True)
    plain_ticks(ax[0], sorted(h))
    ax[0].set_xlabel("element size $h$ (mm)")
    ax[0].set_title("(a) the claim, versus mesh size", fontsize=9, loc="left")
    ax[0].axvline(0.45, color="0.55", lw=0.9, ls=":")
    ax[0].annotate("canal resolved\n($h \\leq$ 0.4 mm)", xy=(0.42, 1.05), fontsize=7,
                   color="0.35", ha="right", va="bottom")

    ax[1].loglog(ntri, [r["halfwidth_pct"] for r in rows], "o-", color=OK, ms=4, lw=1.2)
    ax[1].set_xlabel("number of triangles")
    ax[1].set_ylabel("bracket half-width (%)")
    ax[1].set_title("(b) width of the guaranteed interval", fontsize=9, loc="left")
    ax[1].grid(True, which="both", lw=0.3, alpha=0.4)
    ax[1].yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
    ax[1].yaxis.set_minor_formatter(plt.NullFormatter())
    fig.savefig(os.path.join(FIGDIR, "Fig_S1_mesh.png"))
    plt.close(fig)
    print("Fig_S1_mesh.png")


def fig_gel():
    """The decisive panel: a failure to separate at low sigma_gel is NUMERICAL, not physical."""
    rows = RES["sig_gel_x_mesh"]
    hs = sorted({r["h"] for r in rows}, reverse=True)   # coarse -> fine, left to right
    gels = sorted({r["sig_gel"] for r in rows})
    fig, axes = plt.subplots(1, len(hs), figsize=(2.0 * len(hs) + 0.6, 3.0), sharey=True,
                             layout="constrained")
    for ax, hh in zip(np.atleast_1d(axes), hs):
        sub = [r for r in rows if r["h"] == hh]
        sub.sort(key=lambda r: r["sig_gel"])
        bars(ax, [r["sig_gel"] for r in sub], sub, logx=True)
        pad(ax, [r["sig_gel"] for r in sub], True)
        ax.set_title(f"$h$ = {hh*1e3:.1f} mm", fontsize=9)
        ax.set_xlabel("$\\sigma_{gel}$ (S/m)")
        plain_ticks(ax, gels)
        if ax is not np.atleast_1d(axes)[0]:
            ax.set_ylabel("")
    fig.suptitle("Refining the mesh separates the brackets: the limit is numerical, "
                 "not physical", fontsize=9, y=1.04)
    fig.savefig(os.path.join(FIGDIR, "Fig_S2_gel_x_mesh.png"))
    plt.close(fig)
    print("Fig_S2_gel_x_mesh.png")


def fig_materials():
    spec = [("sig_tis", "$\\sigma_{tissue}$ (S/m)", "(a) body tissue"),
            ("sig_sw", "$\\sigma_{seawater}$ (S/m)", "(b) seawater"),
            ("canal_halfwidth", "canal half-width (mm)", "(c) canal geometry")]
    fig, axes = plt.subplots(1, 3, figsize=(8.2, 3.0), sharey=True, layout="constrained")
    for ax, (key, xlab, title) in zip(axes, spec):
        rows = RES[key]
        x = [r[key] * (1e3 if key == "canal_halfwidth" else 1) for r in rows]
        lg = key != "canal_halfwidth"
        bars(ax, x, rows, logx=lg)
        pad(ax, x, lg)
        plain_ticks(ax, x)
        ax.set_xlabel(xlab)
        ax.set_title(title, fontsize=9, loc="left")
        if ax is not axes[0]:
            ax.set_ylabel("")
    fig.suptitle("The conclusion survives every material and geometric parameter swept",
                 fontsize=9, y=1.03)
    fig.savefig(os.path.join(FIGDIR, "Fig_S3_materials.png"))
    plt.close(fig)
    print("Fig_S3_materials.png")


def fig_tuning():
    """Redo the Monte Carlo here so the figure is self-contained and seeded."""
    rng = np.random.default_rng(0)
    n = 20000
    R_M, C_M, R_W, C_W = 6.6e6, 3.0e-9, 2.0e6, 0.4e-6

    def jitter(x):
        return x * np.exp(rng.uniform(np.log(0.5), np.log(1.5), n))

    tau_m, tau_w = jitter(R_M) * jitter(C_M), jitter(R_W) * jitter(C_W)
    f0 = 1 / (2 * np.pi * np.sqrt(tau_w * tau_m))
    f1, f2 = 1 / (2 * np.pi * tau_w), 1 / (2 * np.pi * tau_m)

    fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.0), layout="constrained")
    ax[0].hist(np.log10(f0), bins=70, color=OK, alpha=0.75)
    for q, ls in ((5, ":"), (50, "-"), (95, ":")):
        ax[0].axvline(np.log10(np.percentile(f0, q)), color="0.2", ls=ls, lw=1.0)
    ax[0].set_xlabel("$\\log_{10}$ peak frequency $f_0$ (Hz)")
    ax[0].set_ylabel("draws")
    ax[0].set_title("(a) peak of the passband, $\\pm$50% on every $RC$", fontsize=9, loc="left")
    ax[0].text(0.02, 0.95, f"median {np.median(f0):.2f} Hz\n"
                           f"[p05 {np.percentile(f0,5):.2f}, p95 {np.percentile(f0,95):.2f}]\n"
                           f"{100*((f0>0.1)&(f0<10)).mean():.1f}% within 0.1-10 Hz",
               transform=ax[0].transAxes, va="top", fontsize=7.5)

    ax[1].loglog(f1, f2, ".", ms=0.7, alpha=0.15, color=OK)
    lim = [min(f1.min(), f2.min()), max(f1.max(), f2.max())]
    ax[1].plot(lim, lim, "--", color="0.35", lw=1.0)
    ax[1].set_xlabel("high-pass corner $f_1$ (Hz)")
    ax[1].set_ylabel("low-pass corner $f_2$ (Hz)")
    ax[1].set_title("(b) a band-pass exists in 100% of draws", fontsize=9, loc="left")
    fig.savefig(os.path.join(FIGDIR, "Fig_S4_tuning.png"))
    plt.close(fig)
    print("Fig_S4_tuning.png")


if __name__ == "__main__":
    os.makedirs(FIGDIR, exist_ok=True)
    fig_mesh()
    fig_gel()
    fig_materials()
    fig_tuning()
