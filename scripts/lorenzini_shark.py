"""Unified electric + magnetic electroreception on a real shark geometry (STL).

The shark body (tissue) is immersed in seawater. An array of Lorenzini ampullae (gel canals)
sits on the snout. Two stimuli are applied with the same certified DGM conduction engine:
  (E) a prey: a bioelectric current dipole in the water near the snout (a local 1/r^3 field);
  (B) navigation: swimming at speed v through the Earth field B sets up a motional field
      E_mot = v x B, uniform over the body.
We read the receptor voltage at each ampulla for both, and report the levels against the known
behavioural threshold (~5 nV/cm). The key point is that the two stimuli leave different spatial
signatures across the array (dipolar vs uniform), which the organ can use to tell prey from
heading. A real anatomical head can replace this STL; the pipeline is unchanged.
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.spatial import Delaunay
import scipy.sparse.linalg as spla
import trimesh
from dgm.mesh3d import TetMesh
from dgm.primal3d import assemble_primal_3d

# Figure output directory. Defaults to the article folder, overridable so the repository
# runs anywhere (a hard-coded absolute path is not reproducible for a third party).
FIGDIR = os.environ.get("LORENZINI_FIGDIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures"))
os.makedirs(FIGDIR, exist_ok=True)


SIG_SW, SIG_TIS, SIG_GEL = 4.0, 0.3, 4.0
L_SHARK = 0.30                 # scale the model to a 0.30 m forebody (m)
NAMP = 6
R_CAN = 0.006                 # canal tube radius (m)
THRESH = 5e-9 / 0.01          # 5 nV/cm -> V/m behavioural threshold
B_EARTH = 50e-6               # Earth field (T)
V_SWIM = 1.0                  # swim speed (m/s)
PREY_DIPOLE = 3e-7            # prey current-dipole moment I*d (A*m), ~ small fish
PREY_RANGE = 0.05            # prey distance in front of the snout (m)


def load_shark():
    f = glob.glob("/Users/talleb/Downloads/shark*.stl")[0]
    m = trimesh.load(f, force="mesh")
    s = L_SHARK / (m.bounds[1, 0] - m.bounds[0, 0])     # scale so x-extent = L_SHARK
    m.apply_scale(s)
    m.apply_translation(-m.bounds.mean(0))               # centre at origin
    return m


def ampullae(mesh):
    """Pores fanned on the snout (the narrow -x end); canals run inward."""
    snout_x = mesh.bounds[0, 0] + 0.04 * (mesh.bounds[1, 0]-mesh.bounds[0, 0])
    out = []
    for k in range(NAMP):
        ph = 2*np.pi*k/NAMP
        r = 0.02
        pore = np.array([snout_x, r*np.cos(ph), r*np.sin(ph)])
        recv = np.array([snout_x + 0.05, 0.01*np.cos(ph), 0.01*np.sin(ph)])
        out.append((ph, pore, recv))
    return out


def build(mesh, h=0.017, seed=0):
    rng = np.random.default_rng(seed)
    b = mesh.bounds; pad = 0.4*(b[1]-b[0])
    lo, hi = b[0]-pad, b[1]+pad
    grid = [np.arange(lo[i], hi[i]+h, h) for i in range(3)]
    GX, GY, GZ = np.meshgrid(*grid); P = np.c_[GX.ravel(), GY.ravel(), GZ.ravel()]
    P += h*0.25*rng.standard_normal(P.shape)
    amp = ampullae(mesh); extra = []
    for _, pore, recv in amp:
        for s in np.linspace(0, 1, 14):
            extra.append(pore + s*(recv-pore) + rng.normal(0, R_CAN*0.5, 3))
    P = np.vstack([P, np.array(extra)])
    P = np.unique(np.round(P, 5), axis=0)
    in_node = mesh.contains(P)                       # inside test on NODES (~6x fewer than tets)
    m = TetMesh(P, Delaunay(P).simplices)
    inside = in_node[m.tets].mean(1) >= 0.5          # tet is tissue if most vertices are inside
    cen = m.points[m.tets].mean(1)
    sig = np.where(inside, SIG_TIS, SIG_SW)
    for _, pore, recv in amp:
        d = recv-pore; L2 = d@d
        t = np.clip(((cen-pore)@d)/L2, 0, 1); proj = pore + t[:, None]*d
        sig = np.where(np.linalg.norm(cen-proj, axis=1) < R_CAN, SIG_GEL, sig)
    return m, sig, amp, (lo, hi)


def factor(m, sig, box):
    K = assemble_primal_3d(m, sig)
    lo, hi = box
    bnd = np.where((np.abs(m.points-lo).min(1) < 1e-6) | (np.abs(m.points-hi).min(1) < 1e-6) |
                   (m.points[:, 0] < lo[0]+1e-6) | (m.points[:, 0] > hi[0]-1e-6) |
                   (m.points[:, 1] < lo[1]+1e-6) | (m.points[:, 1] > hi[1]-1e-6) |
                   (m.points[:, 2] < lo[2]+1e-6) | (m.points[:, 2] > hi[2]-1e-6))[0]
    free = np.setdiff1d(np.arange(m.np), bnd)
    return K, bnd, free, spla.splu(K[free][:, free].tocsc())


def nidx(m, pt): return int(np.argmin(np.sum((m.points-pt)**2, axis=1)))


def amp_voltages(m, amp, phi):
    ref = nidx(m, np.array([0.05, 0, 0]))   # body interior reference
    return np.array([phi[nidx(m, r)] - phi[ref] for _, _, r in amp])


def solve_uniform(m, K, bnd, free, lu, Evec):
    phi = np.zeros(m.np)
    pts = np.nan_to_num(m.points[bnd], posinf=0.0, neginf=0.0)
    phi[bnd] = -(pts @ Evec)
    phi[free] = lu.solve(-(K[free][:, bnd] @ phi[bnd]))
    phi[~np.isfinite(phi)] = 0.0
    return phi


def solve_dipole(m, K, bnd, free, lu, pplus, pminus, I):
    load = np.zeros(m.np); load[nidx(m, pplus)] += I; load[nidx(m, pminus)] -= I
    phi = np.zeros(m.np)
    phi[free] = lu.solve(load[free])
    return phi


if __name__ == "__main__":
    mesh = load_shark()
    print(f"shark scaled to {L_SHARK} m; building mesh ...", flush=True)
    m, sig, amp, box = build(mesh)
    print(f"  np={m.np} nt={m.nt}", flush=True)
    K, bnd, free, lu = factor(m, sig, box)
    # (B) geomagnetic motional field E = v x B (take B vertical, v along +x -> E along y)
    Emot = V_SWIM * B_EARTH                      # V/m
    phiB = solve_uniform(m, K, bnd, free, lu, np.array([0.0, Emot, 0.0]))
    VrB = amp_voltages(m, amp, phiB)
    # (E) prey current dipole near the snout
    snout = np.array([mesh.bounds[0, 0]-PREY_RANGE, 0, 0])
    d = 0.04; I = PREY_DIPOLE / d        # pole separation > mesh size so the two poles are distinct nodes
    phiE = solve_dipole(m, K, bnd, free, lu, snout+[0, d/2, 0], snout-[0, d/2, 0], I)
    VrE = amp_voltages(m, amp, phiE)
    np.savez("/tmp/lorenzini_shark.npz", VrB=VrB, VrE=VrE, phis=[a[0] for a in amp],
             pts=m.points, tets=m.tets, sig=sig, phiB=phiB, phiE=phiE)
    print(f"\n--- LEVELS ---")
    print(f"  (B) navigation: v={V_SWIM} m/s in B={B_EARTH*1e6:.0f} uT -> motional field {Emot*1e6:.1f} uV/m "
          f"= {Emot/THRESH:.0f}x threshold")
    print(f"      ampulla receptor voltages |Vr|: {np.round(np.abs(VrB)*1e9,1)} nV   (uniform signature)")
    print(f"  (E) prey dipole {PREY_DIPOLE*1e9:.1f} nA*m at {PREY_RANGE*100:.0f} cm:")
    print(f"      ampulla receptor voltages |Vr|: {np.round(np.abs(VrE)*1e9,1)} nV   (dipolar signature)")
    cB = np.abs(VrB).max()/max(np.abs(VrB).min(), 1e-18)
    cE = np.abs(VrE).max()/max(np.abs(VrE).min(), 1e-18)
    print(f"  spatial contrast across array (max/min): B={cB:.1f} (uniform)  E={cE:.1f} (localized)")

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    phis = np.array([a[0] for a in amp])
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    x = np.arange(NAMP); w = 0.38
    ax[0].bar(x-w/2, np.abs(VrB)*1e9, w, color="#185FA5", label="B: navigation (v x B, uniform)")
    ax[0].bar(x+w/2, np.abs(VrE)*1e9, w, color="#A32D2D", label="E: prey dipole (localized)")
    ax[0].axhline(THRESH*0.05*1e9, color="k", ls="--", lw=1, label="~receptor floor")
    ax[0].set_yscale("log"); ax[0].set_xlabel("ampulla"); ax[0].set_ylabel("|receptor voltage| (nV)")
    ax[0].set_title(f"(a) Levels per ampulla\nnavigation uniform (contrast {cB:.0f}x) vs prey localized (contrast {cE:.0f}x)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, axis="y")
    axp = fig.add_subplot(1, 2, 2, projection="polar")
    axp.plot(np.r_[phis, phis[0]], np.r_[np.abs(VrB), np.abs(VrB)[0]]*1e9, "o-", color="#185FA5", label="navigation")
    axp.plot(np.r_[phis, phis[0]], np.r_[np.abs(VrE), np.abs(VrE)[0]]*1e9, "s-", color="#A32D2D", label="prey")
    axp.set_title("(b) Spatial signature across the array (nV)\nuniform vs directional", pad=18)
    axp.legend(fontsize=8, loc="upper right", bbox_to_anchor=(1.2, 1.1))
    fig.suptitle("Unified electric + magnetic electroreception on a shark (certified DGM conduction): levels and signatures", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "Fig_shark_EB.png"), dpi=125)
    print("-> saved Fig_shark_EB.png")
