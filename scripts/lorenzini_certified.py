"""Certified bracket on the ampulla access conductance, with the pore and receptor meshed
as conductor boundaries (voids) so the complementary lower bound is rigorous (cf. Article 5).

Pore (surface electrode) at unit potential, receptor (deep electrode) at zero, gel canal
between, tissue and seawater around. The access conductance G = 2W (W the dissipation) is
bracketed by the primal energy (upper) and the RT0 flux (lower): the certified funnelling.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.spatial import Delaunay
from dgm.mesh2d import TriMesh
from dgm import assemble_primal, solve_dirichlet, energy
from dgm.mixed2d import energy_lower_2d

SIG_SW, SIG_TIS, SIG_GEL = 4.0, 0.3, 4.0
LX, LY = 0.040, 0.040
SURF = 0.030                      # tissue below, seawater above
PORE = (0.0195, 0.0205, 0.0298, 0.0304)     # surface electrode (x0,x1,y0,y1)
RECV = (0.0195, 0.0205, 0.0150, 0.0156)     # deep receptor electrode
CANAL = (0.0196, 0.0204)          # canal x-range (gel), between RECV top and PORE bottom


def build(h=0.0006, with_canal=True):
    rng = np.random.default_rng(0)
    gx = np.arange(0, LX+h, h); gy = np.arange(0, LY+h, h)
    GX, GY = np.meshgrid(gx, gy); P = np.c_[GX.ravel(), GY.ravel()]
    it = (P[:, 0] > 1e-9) & (P[:, 0] < LX-1e-9) & (P[:, 1] > 1e-9) & (P[:, 1] < LY-1e-9)
    P[it] += h*0.3*rng.standard_normal(P[it].shape)
    pts = [P]; nb = 14
    for (x0, x1, y0, y1) in (PORE, RECV):
        pts.append(np.c_[np.linspace(x0, x1, nb), np.full(nb, y0)])
        pts.append(np.c_[np.linspace(x0, x1, nb), np.full(nb, y1)])
        pts.append(np.c_[np.full(nb, x0), np.linspace(y0, y1, nb)])
        pts.append(np.c_[np.full(nb, x1), np.linspace(y0, y1, nb)])
    P = np.unique(np.round(np.vstack(pts), 7), axis=0)
    P = P[(P[:, 0] >= -1e-9) & (P[:, 0] <= LX+1e-9) & (P[:, 1] >= -1e-9) & (P[:, 1] <= LY+1e-9)]

    def inr(pt, r, tol=0.0):
        return (pt[:, 0] >= r[0]-tol) & (pt[:, 0] <= r[1]+tol) & (pt[:, 1] >= r[2]-tol) & (pt[:, 1] <= r[3]+tol)
    keep = np.ones(len(P), bool)
    for r in (PORE, RECV):
        keep &= ~inr(P, r, tol=-1e-4)
    P = P[keep]
    tri = Delaunay(P); cen = P[tri.simplices].mean(1)
    good = np.ones(len(tri.simplices), bool)
    for r in (PORE, RECV):
        good &= ~inr(cen, r)
    # Carving out the electrodes can orphan nodes that belonged only to removed triangles.
    # They carry no equation, so the stiffness matrix comes out exactly singular and the
    # solve silently returns NaN. Prune them and reindex. (No effect on the remaining
    # discretisation: an orphan node is a disconnected degree of freedom, not a constraint.)
    tris = tri.simplices[good]
    used = np.zeros(len(P), bool); used[tris.ravel()] = True
    if not used.all():
        remap = np.full(len(P), -1, int); remap[used] = np.arange(int(used.sum()))
        P = P[used]; tris = remap[tris]
    m = TriMesh(P, tris); cen = m.points[m.tris].mean(1)
    sig = np.where(cen[:, 1] < SURF, SIG_TIS, SIG_SW)
    if with_canal:
        incanal = (cen[:, 0] > CANAL[0]) & (cen[:, 0] < CANAL[1]) & (cen[:, 1] > RECV[3]) & (cen[:, 1] < PORE[2])
        sig = np.where(incanal, SIG_GEL, sig)

    def bnodes(r):
        on = (np.abs(P[:, 0]-r[0]) < 1e-7) | (np.abs(P[:, 0]-r[1]) < 1e-7) | \
             (np.abs(P[:, 1]-r[2]) < 1e-7) | (np.abs(P[:, 1]-r[3]) < 1e-7)
        return np.where(on & inr(P, r, tol=1e-7))[0]
    return m, sig, bnodes(PORE), bnodes(RECV)


def bracket(with_canal=True):
    m, sig, pore, recv = build(with_canal=with_canal)
    K = assemble_primal(m, sig)
    bc = {int(i): 1.0 for i in pore}; bc.update({int(i): 0.0 for i in recv})
    v = solve_dirichlet(K, bc)
    return 2*energy_lower_2d(m, sig, bc), 2*energy(K, v)


if __name__ == "__main__":
    lo_c, hi_c = bracket(True); lo_n, hi_n = bracket(False)
    print("certified ampulla access conductance (RT0 flux lower / primal energy upper):")
    print(f"  with gel canal : G in [{lo_c*1e6:.3f}, {hi_c*1e6:.3f}] uS  half-width {100*(hi_c-lo_c)/(hi_c+lo_c):.2f}%")
    print(f"  no canal (tissue only): G in [{lo_n*1e6:.3f}, {hi_n*1e6:.3f}] uS")
    print(f"  -> the gel canal raises the certified access conductance x{0.5*(lo_c+hi_c)/(0.5*(lo_n+hi_n)):.1f}")
