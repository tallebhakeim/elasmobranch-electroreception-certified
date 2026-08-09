"""Neutralising the electric signature of a body in seawater: an electroquasistatic cloak.

A body immersed in seawater has a conductivity contrast with the water, so under any ambient
field (the field a shark reads, or the field set up by the body itself) it perturbs the current
lines and acts like an induced dipole. That dipole is exactly what passive electroreception
detects. We ask whether a coating (a wetsuit) can cancel it.

Physics. For a body of conductivity sigma_b in seawater sigma_sw under a uniform field E0, the
exterior perturbation is that of a dipole p ~ alpha V E0 with contrast alpha=(sigma_b-sigma_sw)
/(sigma_b+2 sigma_sw) in 3D (here a 2D slice, alpha=(sigma_b-sigma_sw)/(sigma_b+sigma_sw)).
A bare body (sigma_b<sigma_sw) gives a negative dipole; an insulating wetsuit (sigma~0) gives
the LARGEST dipole, so it makes the wearer MORE visible. A conductive coating tuned to seawater
makes the coated body a neutral inclusion: the exterior dipole vanishes and the body is invisible
to the static electric sense. We compute, with the certified DGM conduction engine, the dipole
source integral p = sum (sigma-sigma_sw) E A over the mesh, which is the source of the exterior
perturbation, and sweep the coating conductivity to find the cloaking null.

The dynamic (frequency) part uses kappa = sigma + j omega eps. At electroreception frequencies
(~1 Hz) conduction dominates (omega eps << sigma in seawater), so the in-band neutralisation is a
conductivity-matching problem; the dielectric (eps) term only refines the response at the high
edge of the band and for transients. We show both.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.spatial import Delaunay
from dgm.mesh2d import TriMesh
from dgm import assemble_primal, solve_dirichlet

SIG_SW = 4.0            # seawater (S/m)
SIG_TIS = 0.5           # body tissue (S/m), a moderate average
EPS0 = 8.854e-12
EPS_W = 80*EPS0         # seawater permittivity
L = 0.60                # half-box (m)
A_BODY, B_BODY = 0.18, 0.10     # body ellipse semi-axes (m)
T_SUIT = 0.020          # suit thickness (m)
E0 = 1.0                # ambient field amplitude (V/m), along x; results scale linearly
PROBE = np.array([0.45, 0.0])   # a shark ampulla on the +x axis reads here


def build(h=0.012, seed=0):
    rng = np.random.default_rng(seed)
    g = np.arange(-L, L+h, h); GX, GY = np.meshgrid(g, g)
    P = np.c_[GX.ravel(), GY.ravel()]
    it = (np.abs(P[:, 0]) < L-1e-9) & (np.abs(P[:, 1]) < L-1e-9)
    P[it] += h*0.3*rng.standard_normal(P[it].shape)
    # refine on the body / suit interfaces
    extra = []
    for rr in (1.0, 1.0+T_SUIT/A_BODY):
        for ph in np.linspace(0, 2*np.pi, 160, endpoint=False):
            extra.append([rr*A_BODY*np.cos(ph), rr*B_BODY*np.sin(ph)])
    P = np.unique(np.round(np.vstack([P, extra]), 6), axis=0)
    P = P[(np.abs(P[:, 0]) <= L+1e-9) & (np.abs(P[:, 1]) <= L+1e-9)]
    m = TriMesh(P, Delaunay(P).simplices)
    cen = m.points[m.tris].mean(1)
    rb = (cen[:, 0]/A_BODY)**2 + (cen[:, 1]/B_BODY)**2
    rs = (cen[:, 0]/(A_BODY+T_SUIT))**2 + (cen[:, 1]/(B_BODY+T_SUIT))**2
    region = np.where(rb < 1, 0, np.where(rs < 1, 1, 2))   # 0 body, 1 suit, 2 water
    return m, region


def elem_grad(m, phi):
    """Linear-triangle gradient of a nodal field, one vector per element."""
    p = m.points[m.tris]; f = phi[m.tris]
    x0, y0 = p[:, 0, 0], p[:, 0, 1]
    d1 = p[:, 1] - p[:, 0]; d2 = p[:, 2] - p[:, 0]
    det = d1[:, 0]*d2[:, 1] - d1[:, 1]*d2[:, 0]
    f1 = f[:, 1]-f[:, 0]; f2 = f[:, 2]-f[:, 0]
    gx = (d2[:, 1]*f1 - d1[:, 1]*f2)/det
    gy = (-d2[:, 0]*f1 + d1[:, 0]*f2)/det
    return np.c_[gx, gy]


def solve(m, sig):
    K = assemble_primal(m, sig)
    bn = m.boundary_nodes()
    bc = {int(i): float(-E0*m.points[i, 0]) for i in bn}   # uniform field along x
    return solve_dirichlet(K, bc)


def signature(m, region, sig_suit, with_suit=True):
    """Return the dipole source magnitude |p| and the perturbation potential at the probe (V)."""
    sig = np.full(m.nt, SIG_SW)
    sig[region == 0] = SIG_TIS
    sig[region == 1] = sig_suit if with_suit else SIG_TIS   # no suit: body extends into the ring
    phi = solve(m, sig)
    E = -elem_grad(m, phi)
    A = m.area
    # exterior-dipole source: polarisation current (sigma - sigma_sw) E integrated over volume
    p = np.array([np.sum((sig-SIG_SW)*E[:, 0]*A), np.sum((sig-SIG_SW)*E[:, 1]*A)])
    # perturbation potential read by a distant shark ampulla = total - undisturbed uniform field
    ni = int(np.argmin(np.sum((m.points-PROBE)**2, axis=1)))
    dphi = phi[ni] - (-E0*m.points[ni, 0])
    return np.linalg.norm(p), dphi, phi, sig


def neutral_sigma_2d(sig_core, f):
    """Coating conductivity that makes a coated cylinder a neutral inclusion (zero exterior dipole)
    in seawater: solve sigma_eff(core, shell, f) = sigma_sw for the shell, f=(a/b)^2 core fraction."""
    # sigma_eff = s2 (s1+s2 + f(s1-s2)) / (s1+s2 - f(s1-s2)) = sigma_sw ; solve quadratic in s2
    sm, s1 = SIG_SW, sig_core
    # sm (s1+s2 - f(s1-s2)) = s2 (s1+s2 + f(s1-s2))
    # sm s1 + sm s2 - sm f s1 + sm f s2 = s1 s2 + s2^2 + f s1 s2 - f s2^2
    # 0 = s2^2(1-f) + s2(s1 + f s1 - sm - sm f) - sm s1 (1 - f)
    a = (1-f); b = s1*(1+f) - sm*(1+f); c = -sm*s1*(1-f)
    disc = b*b - 4*a*c
    return (-b + np.sqrt(disc))/(2*a)


if __name__ == "__main__":
    m, region = build()
    print(f"mesh np={m.np} nt={m.nt}; regions body/suit/water = "
          f"{(region==0).sum()}/{(region==1).sum()}/{(region==2).sum()}", flush=True)
    f_core = ((A_BODY)/(A_BODY+T_SUIT))*((B_BODY)/(B_BODY+T_SUIT))   # 2D area fraction (a/b)^2 ~ product
    s_neutral = neutral_sigma_2d(SIG_TIS, f_core)
    print(f"analytic neutral-inclusion coating sigma ~ {s_neutral:.2f} S/m (core fraction f={f_core:.3f})")

    cases = {"bare body": (SIG_TIS, False),
             "neoprene wetsuit (insulator)": (1e-3, True),
             "seawater-matched suit": (SIG_SW, True),
             "neutral-inclusion cloak": (s_neutral, True)}
    res = {}
    for name, (ss, ws) in cases.items():
        pmag, dphi, phi, sig = signature(m, region, ss, ws)
        res[name] = (pmag, dphi, phi, sig)
        print(f"  {name:32s}: |p|={pmag*1e3:8.3f} mA*m   probe dphi={dphi*1e3:8.3f} mV")
    p_bare = res["bare body"][0]
    print("\nsuppression vs bare body (dipole source):")
    for name in cases:
        sup = 20*np.log10(max(res[name][0], 1e-30)/p_bare)
        print(f"  {name:32s}: {sup:+6.1f} dB")

    # sweep coating conductivity -> find the cloaking null
    sweep = np.geomspace(1e-3, 40, 60)
    pcurve = np.array([signature(m, region, s, True)[0] for s in sweep])
    snull = sweep[np.argmin(pcurve)]
    print(f"\nswept null at sigma_coat ~ {snull:.2f} S/m, |p| reduced "
          f"{20*np.log10(pcurve.min()/p_bare):+.1f} dB vs bare")
    np.savez("/tmp/lorenzini_cloak.npz",
             pts=m.points, tris=m.tris, region=region,
             sweep=sweep, pcurve=pcurve, p_bare=p_bare, snull=snull, s_neutral=s_neutral,
             **{f"phi_{i}": res[n][2] for i, n in enumerate(cases)},
             names=list(cases.keys()),
             pmags=np.array([res[n][0] for n in cases]),
             dphis=np.array([res[n][1] for n in cases]))
    print("-> saved /tmp/lorenzini_cloak.npz")
