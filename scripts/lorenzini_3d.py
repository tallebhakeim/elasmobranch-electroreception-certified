"""3D electroreception on a simplified ray head: an array of Lorenzini ampullae.

A flattened ellipsoid (a ray is dorso-ventrally flattened) of body tissue sits in seawater.
An array of gel-filled ampulla canals fans out from a central ventral cluster, each canal
axis pointing a different azimuth. A uniform far field at azimuth theta drives the conduction
problem; the high-conductivity canals funnel the field to their receptors, and each ampulla is
tuned to the field direction aligned with its canal. We sweep the field azimuth and read each
receptor, giving a directional sensitivity rosette. A real ray/shark head STL can replace the
ellipsoid (gmsh meshes it the same way); the ellipsoid is the placeholder geometry.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.spatial import Delaunay
import scipy.sparse.linalg as spla
from dgm.mesh3d import TetMesh
from dgm.primal3d import assemble_primal_3d

# Figure output directory. Defaults to the article folder, overridable so the repository
# runs anywhere (a hard-coded absolute path is not reproducible for a third party).
FIGDIR = os.environ.get("LORENZINI_FIGDIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures"))
os.makedirs(FIGDIR, exist_ok=True)


SIG_SW, SIG_TIS, SIG_GEL = 4.0, 0.3, 4.0
A, B, C = 0.090, 0.065, 0.018          # head ellipsoid semi-axes (m): wide, flat
BOX = (0.13, 0.10, 0.055)              # seawater half-box (m)
NAMP = 6                               # ampullae (canals) fanning in azimuth
R_PORE = 0.055                         # pore ring radius on the ventral surface
R_CAN = 0.004                          # canal "tube" radius (gel)


def ampullae():
    """Each ampulla: (pore on ventral surface, receptor near centre), canal axis = its azimuth."""
    out = []
    for k in range(NAMP):
        ph = 2 * np.pi * k / NAMP
        zp = -C * np.sqrt(max(0.0, 1 - (R_PORE * np.cos(ph) / A) ** 2 - (R_PORE * np.sin(ph) / B) ** 2))
        pore = np.array([R_PORE * np.cos(ph), R_PORE * np.sin(ph), zp + 0.002])
        recv = np.array([0.012 * np.cos(ph), 0.012 * np.sin(ph), -0.004])
        out.append((ph, pore, recv))
    return out


def build(h=0.006, seed=0):
    rng = np.random.default_rng(seed)
    nx, ny, nz = [int(2 * b / h) + 1 for b in BOX]
    gx = np.linspace(-BOX[0], BOX[0], nx); gy = np.linspace(-BOX[1], BOX[1], ny); gz = np.linspace(-BOX[2], BOX[2], nz)
    GX, GY, GZ = np.meshgrid(gx, gy, gz); P = np.c_[GX.ravel(), GY.ravel(), GZ.ravel()]
    interior = (np.abs(P[:, 0]) < BOX[0]-1e-9) & (np.abs(P[:, 1]) < BOX[1]-1e-9) & (np.abs(P[:, 2]) < BOX[2]-1e-9)
    P[interior] += h * 0.3 * rng.standard_normal(P[interior].shape)
    # densify along each canal
    amp = ampullae(); extra = []
    for _, pore, recv in amp:
        for s in np.linspace(0, 1, 12):
            extra.append(pore + s * (recv - pore) + rng.normal(0, R_CAN*0.5, 3))
    P = np.vstack([P, np.array(extra)])
    P = np.unique(np.round(P, 5), axis=0)
    m = TetMesh(P, Delaunay(P).simplices)
    with np.errstate(over="ignore", invalid="ignore"):
        cen = m.points[m.tets].mean(1)
    cen[~np.isfinite(cen)] = 0.0
    ell = (cen[:, 0]/A)**2 + (cen[:, 1]/B)**2 + (cen[:, 2]/C)**2
    sig = np.where(ell < 1.0, SIG_TIS, SIG_SW)                  # head tissue / seawater
    for _, pore, recv in amp:                                   # carve high-sigma canals
        d = recv - pore; L2 = d @ d
        t = np.clip(((cen - pore) @ d) / L2, 0, 1)
        proj = pore + t[:, None] * d
        sig = np.where(np.linalg.norm(cen - proj, axis=1) < R_CAN, SIG_GEL, sig)
    return m, sig, amp


def solve_angles(m, sig, amp, thetas, E0=1.0):
    K = assemble_primal_3d(m, sig)
    bnd = np.where((np.abs(np.abs(m.points[:, 0]) - BOX[0]) < 1e-6) |
                   (np.abs(np.abs(m.points[:, 1]) - BOX[1]) < 1e-6) |
                   (np.abs(np.abs(m.points[:, 2]) - BOX[2]) < 1e-6))[0]
    free = np.setdiff1d(np.arange(m.np), bnd)
    lu = spla.splu(K[free][:, free].tocsc())
    def nidx(pt): return int(np.argmin(np.sum((m.points - pt)**2, axis=1)))
    ridx = [nidx(r) for _, _, r in amp]; refi = nidx(np.array([0, 0, 0.0]))
    Vr = np.zeros((len(thetas), NAMP))
    for j, th in enumerate(thetas):
        phi = np.zeros(m.np)
        phi[bnd] = -E0 * (m.points[bnd, 0]*np.cos(th) + m.points[bnd, 1]*np.sin(th))
        rhs = -(K[free][:, bnd] @ phi[bnd])
        phi[free] = lu.solve(rhs)
        Vr[j] = phi[ridx] - phi[refi]
    return Vr, ridx


def figure():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    m, sig, amp = build()
    th = np.linspace(0, 2*np.pi, 49)
    Vr, ridx = solve_angles(m, sig, amp, th)
    cols = plt.cm.hsv(np.linspace(0, 1, NAMP, endpoint=False))
    fig = plt.figure(figsize=(13.5, 5.2))
    # (a) 3D head + canals
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    ax.plot_surface(A*np.cos(u)*np.sin(v)*1e3, B*np.sin(u)*np.sin(v)*1e3, C*np.cos(v)*1e3,
                    color="#cdbb9a", alpha=.25, linewidth=0)
    for k, (ph, pore, recv) in enumerate(amp):
        ax.plot(*[[pore[i]*1e3, recv[i]*1e3] for i in range(3)], color=cols[k], lw=3)
        ax.scatter(*[pore[i]*1e3 for i in range(3)], color=cols[k], s=25)
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.set_zlabel("z (mm)")
    ax.set_box_aspect((A, B, C*2.5)); ax.view_init(elev=22, azim=-60)
    ax.set_title("(a) Simplified ray head (flattened ellipsoid)\nwith an array of gel ampulla canals")
    # (b) directional rosette
    axp = fig.add_subplot(1, 2, 2, projection="polar")
    for k in range(NAMP):
        axp.plot(th, np.abs(Vr[:, k])*1e3, color=cols[k], lw=1.8, label=f"{np.degrees(amp[k][0]):.0f} deg")
    axp.set_title("(b) Directional tuning |Vr| (mV per V/m)\neach ampulla peaks along its canal", pad=18)
    axp.legend(fontsize=7, loc="upper right", bbox_to_anchor=(1.22, 1.1), title="canal azimuth")
    fig.suptitle("3D electroreception on a ray head: an array of certified Lorenzini ampullae (dual DGM conduction)", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "Fig_lorenzini_3d.png"), dpi=125)
    print("-> saved Fig_lorenzini_3d.png")
    for k in range(NAMP):
        best = np.degrees(th[np.argmax(np.abs(Vr[:, k]))])
        print(f"    ampulla {k}: canal={np.degrees(amp[k][0]):3.0f}  peak {best:3.0f}  |Vr|max={np.abs(Vr[:,k]).max()*1e3:.2f} mV/(V/m)")


if __name__ == "__main__":
    figure()
