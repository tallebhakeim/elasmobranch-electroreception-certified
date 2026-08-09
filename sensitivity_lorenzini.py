# -*- coding: utf-8 -*-
"""Sensitivity and robustness campaign for the certified ampulla access conductance.

Written in response to the JTB editorial decision (JTB-D-26-01233), which asked for
"the necessary sensitivity & robustness analyses to ensure that the authors and reviewers
have a strong handle on claims".

The claim under test is NOT the value of the access conductance. It is the comparison:
does the gel canal raise it? Because the discrete geometric method returns a guaranteed
two-sided bracket [G_lo, G_hi] rather than a point value, that comparison can be settled
rigorously: if the bracket with the canal and the bracket without it are disjoint, the
conclusion holds for the exact solution.

We therefore report, for every parameter set, the GUARANTEED bracket on the gain

    gain = G_canal / G_nocanal  in  [ G_canal_lo / G_nocanal_hi ,  G_canal_hi / G_nocanal_lo ]

and whether its lower end stays above 1. That is the quantity a reviewer should look at.

Sweeps: mesh refinement, gel conductivity, tissue conductivity, seawater conductivity,
canal width, and (analytically) the membrane and wall time constants behind the tuning.

Run:  python3 sensitivity_lorenzini.py            (writes sensitivity_results.json)
      python3 sensitivity_lorenzini.py --quick    (coarse sweeps, for a smoke test)
"""
import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                          # dgm/ lives next to this file
sys.path.insert(0, os.path.join(_HERE, "scripts"))  # the study scripts

import lorenzini_certified as LC          # noqa: E402
from dgm import assemble_primal, solve_dirichlet, energy    # noqa: E402
from dgm.mixed2d import energy_lower_2d                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sensitivity_results.json")

# Reference parameter set (the one used in the manuscript).
# Reference mesh is h = 0.4 mm, i.e. one element across the canal half-width. The manuscript
# submitted to JTB used h = 0.6 mm, which does NOT resolve the canal: it gives a 13% bracket
# and, at low gel conductivity, brackets that fail to separate. See the sig_gel x h cross sweep.
REF = dict(sig_sw=4.0, sig_tis=0.3, sig_gel=4.0, h=0.0004, canal_halfwidth=0.0004)


def bracket(with_canal=True, h=None, sig_sw=None, sig_tis=None, sig_gel=None,
            canal_halfwidth=None):
    """Guaranteed [lower, upper] on the 2D access conductance, per unit depth (S/m x m).

    The lower bound comes from an RT0 divergence-conforming flux (complementary energy),
    the upper bound from the primal nodal potential. Both are computed on the same mesh;
    neither is an estimate.
    """
    old = (LC.SIG_SW, LC.SIG_TIS, LC.SIG_GEL, LC.CANAL)
    try:
        if sig_sw is not None:
            LC.SIG_SW = sig_sw
        if sig_tis is not None:
            LC.SIG_TIS = sig_tis
        if sig_gel is not None:
            LC.SIG_GEL = sig_gel
        if canal_halfwidth is not None:
            xc = 0.5 * (LC.CANAL[0] + LC.CANAL[1])
            LC.CANAL = (xc - canal_halfwidth, xc + canal_halfwidth)
        m, sig, pore, recv = LC.build(h=h or REF["h"], with_canal=with_canal)
        K = assemble_primal(m, sig)
        bc = {int(i): 1.0 for i in pore}
        bc.update({int(i): 0.0 for i in recv})
        v = solve_dirichlet(K, bc)
        return 2 * energy_lower_2d(m, sig, bc), 2 * energy(K, v), len(m.tris)
    finally:
        LC.SIG_SW, LC.SIG_TIS, LC.SIG_GEL, LC.CANAL = old


def gain_bracket(**kw):
    """Guaranteed bracket on G_canal / G_nocanal, plus the two underlying brackets."""
    lo_c, hi_c, ntri = bracket(with_canal=True, **kw)
    lo_n, hi_n, _ = bracket(with_canal=False, **kw)
    return dict(
        G_canal=[lo_c, hi_c], G_nocanal=[lo_n, hi_n],
        halfwidth_pct=100 * (hi_c - lo_c) / (hi_c + lo_c),
        gain=[lo_c / hi_n, hi_c / lo_n],
        disjoint=bool(lo_c > hi_n),          # the certified statement: canal strictly helps
        ntri=ntri,
    )


def sweep(name, key, values, quick=False):
    rows = []
    print(f"\n=== {name} ===")
    print(f"{key:>12} | {'G_canal (S/m depth)':>22} | {'G_nocanal (S/m depth)':>22} | "
          f"{'gain bracket':>18} | {'width%':>7} | disjoint")
    for val in values:
        t = time.time()
        r = gain_bracket(**{key: val})
        r[key] = val
        r["seconds"] = round(time.time() - t, 2)
        rows.append(r)
        gc, gn, g = r["G_canal"], r["G_nocanal"], r["gain"]
        print(f"{val:>12.5g} | [{gc[0]:9.4f},{gc[1]:9.4f}] | "
              f"[{gn[0]:9.4f},{gn[1]:9.4f}] | "
              f"[{g[0]:6.3f},{g[1]:6.3f}] | {r['halfwidth_pct']:6.2f} | "
              f"{'YES' if r['disjoint'] else 'no'}")
    return rows


def cross_sweep(name, k1, v1, k2, v2):
    """Two-way sweep, to separate a physical limit from a resolution limit."""
    rows = []
    print(f"\n=== {name} ===")
    print(f"{k1:>10} {k2:>10} {'ntri':>7} {'gain bracket':>18} {'width%':>7}  disjoint")
    for a in v1:
        for b in v2:
            r = gain_bracket(**{k1: a, k2: b})
            r[k1], r[k2] = a, b
            rows.append(r)
            g = r["gain"]
            print(f"{a:>10.5g} {b:>10.5g} {r['ntri']:>7d} [{g[0]:7.3f},{g[1]:7.3f}] "
                  f"{r['halfwidth_pct']:>7.2f}  {'YES' if r['disjoint'] else 'no'}")
    return rows


def tuning_robustness(n=20000, seed=0):
    """The band-pass corners and peak under +/-50% on each RC, log-uniform.

    Analytic, so this is a distribution rather than a bracket: it answers "how much of the
    tuning conclusion survives the fact that the RC values are representative, not fitted".
    """
    rng = np.random.default_rng(seed)
    R_M, C_M, R_W, C_W = 6.6e6, 3.0e-9, 2.0e6, 0.4e-6

    def jitter(x):
        return x * np.exp(rng.uniform(np.log(0.5), np.log(1.5), n))

    tau_m = jitter(R_M) * jitter(C_M)
    tau_w = jitter(R_W) * jitter(C_W)
    f1 = 1 / (2 * np.pi * tau_w)
    f2 = 1 / (2 * np.pi * tau_m)
    f0 = 1 / (2 * np.pi * np.sqrt(tau_w * tau_m))
    band = f2 > f1                       # a band-pass exists at all
    out = {}
    for lbl, arr in (("f1_highpass_Hz", f1), ("f2_lowpass_Hz", f2), ("f0_peak_Hz", f0)):
        out[lbl] = dict(median=float(np.median(arr)),
                        p05=float(np.percentile(arr, 5)),
                        p95=float(np.percentile(arr, 95)))
    out["fraction_bandpass"] = float(band.mean())
    out["fraction_peak_below_10Hz"] = float((f0 < 10).mean())
    out["fraction_peak_in_0p1_10Hz"] = float(((f0 > 0.1) & (f0 < 10)).mean())
    print("\n=== frequency tuning under +/-50% on every RC (log-uniform, "
          f"n={n}) ===")
    for k, v in out.items():
        if isinstance(v, dict):
            print(f"  {k:>18}: median {v['median']:8.3f}  [p05 {v['p05']:.3f}, "
                  f"p95 {v['p95']:.3f}]")
        else:
            print(f"  {k:>18}: {v:.4f}")
    return out


def main():
    quick = "--quick" in sys.argv
    res = {"reference": REF, "note": "2D plane problem: G is a conductance PER UNIT DEPTH, units S/m. It is not the 3D access conductance of an ampulla and must not be quoted in siemens."}

    t0 = time.time()
    ref = gain_bracket()
    res["reference_result"] = ref
    print("=== reference parameter set ===")
    print(f"  G with canal    : [{ref['G_canal'][0]:.4f}, {ref['G_canal'][1]:.4f}] S per metre "
          f"of depth  (half-width {ref['halfwidth_pct']:.2f}%)")
    print(f"  G without canal : [{ref['G_nocanal'][0]:.4f}, {ref['G_nocanal'][1]:.4f}] S per metre of depth")
    print(f"  GUARANTEED gain : [{ref['gain'][0]:.3f}, {ref['gain'][1]:.3f}]  "
          f"disjoint={ref['disjoint']}")

    hs = [0.0012, 0.0008, 0.0006] if quick else [0.0016, 0.0012, 0.0010, 0.0008, 0.0006, 0.0005, 0.0004]
    res["mesh"] = sweep("mesh refinement", "h", hs)

    gels = [0.3, 1.0, 4.0] if quick else [0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    res["sig_gel"] = sweep("gel conductivity (S/m)", "sig_gel", gels)

    tis = [0.1, 0.3, 1.0] if quick else [0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2]
    res["sig_tis"] = sweep("tissue conductivity (S/m)", "sig_tis", tis)

    sws = [3.0, 4.0, 5.0] if quick else [3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
    res["sig_sw"] = sweep("seawater conductivity (S/m)", "sig_sw", sws)

    cws = [0.0002, 0.0004] if quick else [0.00015, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006]
    res["canal_halfwidth"] = sweep("canal half-width (m)", "canal_halfwidth", cws)

    res["sig_gel_x_mesh"] = cross_sweep(
        "gel conductivity x mesh size: is the failure physical or numerical?",
        "sig_gel", [0.5, 1.0, 2.0, 4.0] if not quick else [1.0],
        "h", [0.0008, 0.0006, 0.0005, 0.0004, 0.0003] if not quick else [0.0006])

    res["tuning"] = tuning_robustness(n=2000 if quick else 20000)

    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    print(f"\nwrote {os.path.basename(OUT)}  ({res['elapsed_s']} s)")

    # The headline the manuscript should carry.
    allrows = [r for k in ("mesh", "sig_gel", "sig_tis", "sig_sw", "canal_halfwidth")
               for r in res[k]]
    bad = [r for r in allrows if not r["disjoint"]]
    print(f"\nSUMMARY: {len(allrows)-len(bad)}/{len(allrows)} parameter sets give DISJOINT "
          f"brackets (the canal strictly helps, guaranteed).")
    if bad:
        print("  not disjoint at:")
        for r in bad:
            k = [x for x in ("h", "sig_gel", "sig_tis", "sig_sw", "canal_halfwidth") if x in r]
            print(f"    {k[0]}={r[k[0]]:.5g}  gain in [{r['gain'][0]:.3f}, {r['gain'][1]:.3f}]")


if __name__ == "__main__":
    main()
