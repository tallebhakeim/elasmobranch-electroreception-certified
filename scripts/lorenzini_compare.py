"""Ray vs shark electroreception: morphology and ampulla placement (Article 7).

Two parametric body plans, same conduction engine:
  SHARK: fusiform (prolate ellipsoid), ampulla pores fanned over the SNOUT, canals radiating
         inward -> forward-biased sensitivity (hunts prey ahead).
  RAY:   flattened disc (oblate ellipsoid), ampulla pores spread over the VENTRAL surface,
         canals going up -> downward-biased sensitivity (hunts benthic prey below).
For each we apply a prey current dipole (placed where the animal hunts) and the swimming
motional field (v x B, navigation), and read the ampulla array. The flat ventral ray array and
the forward snout shark array give different detection geometries, quantified here.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.spatial import Delaunay
import scipy.sparse.linalg as spla
from dgm.mesh3d import TetMesh
from dgm.primal3d import assemble_primal_3d

SIG_SW, SIG_TIS, SIG_GEL = 4.0, 0.3, 4.0
R_CAN = 0.005
B_EARTH, V_SWIM = 50e-6, 1.0
THRESH = 5e-9/0.01
PREY_DIPOLE = 3e-7

BODIES = {
    "shark": dict(abc=(0.16, 0.045, 0.050), namp=7),   # fusiform, long along x
    "ray":   dict(abc=(0.11, 0.095, 0.022), namp=7),   # flat disc
}


def ampullae(kind, abc):
    a, b, c = abc; out = []
    n = BODIES[kind]["namp"]
    if kind == "shark":                       # pores fanned on the snout (front, +x), canals inward
        for k in range(n):
            ph = 2*np.pi*k/n; r = 0.6
            pore = np.array([0.92*a, r*b*np.cos(ph), r*c*np.sin(ph)])
            recv = np.array([0.55*a, 0.25*b*np.cos(ph), 0.25*c*np.sin(ph)])
            out.append((pore, recv))
    else:                                     # ray: pores over the ventral surface (-z), canals up
        pts = [(0.0, 0.0)] + [(0.55*np.cos(2*np.pi*k/(n-1)), 0.55*np.sin(2*np.pi*k/(n-1))) for k in range(n-1)]
        for (ux, uy) in pts:
            x, y = ux*a, uy*b
            zb = -c*np.sqrt(max(0.0, 1-(x/a)**2-(y/b)**2))
            pore = np.array([x, y, zb+0.002]); recv = np.array([x*0.6, y*0.6, zb+0.018])
            out.append((pore, recv))
    return out


def build(kind, h=0.015, seed=0):
    a, b, c = BODIES[kind]["abc"]; amp = ampullae(kind, (a, b, c))
    rng = np.random.default_rng(seed)
    box = np.array([a, b, c])*2.0 + 0.04
    g = [np.arange(-box[i], box[i]+h, h) for i in range(3)]
    GX, GY, GZ = np.meshgrid(*g); P = np.c_[GX.ravel(), GY.ravel(), GZ.ravel()]
    P += h*0.25*rng.standard_normal(P.shape)
    extra = []
    for pore, recv in amp:
        for s in np.linspace(0, 1, 12):
            extra.append(pore + s*(recv-pore) + rng.normal(0, R_CAN*0.5, 3))
    P = np.unique(np.round(np.vstack([P, np.array(extra)]), 5), axis=0)
    m = TetMesh(P, Delaunay(P).simplices)
    cen = m.points[m.tets].mean(1); cen[~np.isfinite(cen)] = 1e9
    ell = (cen[:, 0]/a)**2 + (cen[:, 1]/b)**2 + (cen[:, 2]/c)**2
    sig = np.where(ell < 1.0, SIG_TIS, SIG_SW)
    for pore, recv in amp:
        d = recv-pore; L2 = d@d
        t = np.clip(((cen-pore)@d)/L2, 0, 1); proj = pore + t[:, None]*d
        sig = np.where(np.linalg.norm(cen-proj, axis=1) < R_CAN, SIG_GEL, sig)
    return m, sig, amp, box


def nidx(m, p): return int(np.argmin(np.sum((m.points-p)**2, axis=1)))


def run_animal(kind):
    m, sig, amp, box = build(kind)
    a, b, c = BODIES[kind]["abc"]
    K = assemble_primal_3d(m, sig)
    bnd = np.where(np.abs(np.abs(m.points)-box).min(1) < 1e-6)[0]
    bnd = np.where((m.points[:, 0] <= -box[0]+1e-6) | (m.points[:, 0] >= box[0]-1e-6) |
                   (m.points[:, 1] <= -box[1]+1e-6) | (m.points[:, 1] >= box[1]-1e-6) |
                   (m.points[:, 2] <= -box[2]+1e-6) | (m.points[:, 2] >= box[2]-1e-6))[0]
    free = np.setdiff1d(np.arange(m.np), bnd)
    lu = spla.splu(K[free][:, free].tocsc())
    ridx = [nidx(m, r) for _, r in amp]; refi = nidx(m, np.array([0.0, 0, 0]))
    def ampV(phi): return phi[ridx]-phi[refi]
    # (B) motional field (uniform, along y)
    phiB = np.zeros(m.np); phiB[bnd] = -(np.nan_to_num(m.points[bnd])[:, 1]*V_SWIM*B_EARTH)
    phiB[free] = lu.solve(-(K[free][:, bnd] @ phiB[bnd])); phiB[~np.isfinite(phiB)] = 0
    VrB = ampV(phiB)
    # (E) prey dipole where the animal hunts: shark ahead (+x), ray below (-z)
    if kind == "shark": prey = np.array([a+0.05, 0, 0])
    else: prey = np.array([0, 0, -c-0.05])
    dvec = np.array([0, 0.04, 0]); I = PREY_DIPOLE/0.04
    load = np.zeros(m.np); load[nidx(m, prey+dvec/2)] += I; load[nidx(m, prey-dvec/2)] -= I
    phiE = np.zeros(m.np); phiE[free] = lu.solve(load[free]); phiE[~np.isfinite(phiE)] = 0
    VrE = ampV(phiE)
    return dict(kind=kind, amp=amp, VrB=VrB, VrE=VrE, np=m.np, abc=(a, b, c))


if __name__ == "__main__":
    res = {k: run_animal(k) for k in ["shark", "ray"]}
    np.savez("/tmp/lorenzini_compare.npz", **{f"{k}_{q}": np.array(res[k][q], dtype=object)
             for k in res for q in ("VrB", "VrE")})
    for k, r in res.items():
        print(f"\n{k.upper()} (np={r['np']}, abc={tuple(round(x,3) for x in r['abc'])})")
        print(f"  prey  |Vr|: {np.round(np.abs(r['VrE'])*1e9,1)} nV  (max {np.abs(r['VrE']).max()*1e9:.0f}, contrast {np.abs(r['VrE']).max()/max(np.abs(r['VrE']).min(),1e-18):.1f}x)")
        print(f"  navig |Vr|: {np.round(np.abs(r['VrB'])*1e9,1)} nV  (max {np.abs(r['VrB']).max()*1e9:.0f}, contrast {np.abs(r['VrB']).max()/max(np.abs(r['VrB']).min(),1e-18):.1f}x)")
