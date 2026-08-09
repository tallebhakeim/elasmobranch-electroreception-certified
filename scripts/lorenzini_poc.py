"""Ampulla of Lorenzini as a certified electroquasistatic conduction problem (proof of concept).

The electroreceptor is a high-conductivity gel canal running from a surface pore into the body
to a sensory epithelium. In conductive surroundings (seawater) a weak external electric field
sets up a potential; the gel canal, almost as conductive as seawater, carries the pore
potential to the receptor, so the sensory membrane sees the external potential drop along the
canal while the resistive body stays nearly equipotential. This funnelling, and its strong
directionality, is what makes the organ exquisitely sensitive. We model it with the dual DGM
conduction solver: seawater, body tissue and gel canal as three conductivities, a uniform far
field as the stimulus, and the receptor voltage as the output.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.spatial import Delaunay
from dgm.mesh2d import TriMesh
from dgm import assemble_primal, solve_dirichlet

# Figure output directory. Defaults to the article folder, overridable so the repository
# runs anywhere (a hard-coded absolute path is not reproducible for a third party).
FIGDIR = os.environ.get("LORENZINI_FIGDIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures"))
os.makedirs(FIGDIR, exist_ok=True)


SIG_SW = 4.0          # seawater conductivity (S/m)
SIG_TIS = 0.3         # body tissue (S/m)
SIG_GEL = 4.0         # ampulla gel (highly conductive, ~ seawater)
LX, LY = 0.08, 0.06   # domain (m)
BODY_TOP = 0.040      # body (tissue) fills y < BODY_TOP; seawater above
PORE = (0.040, 0.040) # pore on the body surface
RECEPTOR = (0.040, 0.018)   # receptor depth (canal runs down from the pore)
CANAL_W = 0.0015      # canal half-width (m)


def build(h=0.0009, with_canal=True):
    rng = np.random.default_rng(0)
    gx = np.arange(0, LX + h, h); gy = np.arange(0, LY + h, h)
    GX, GY = np.meshgrid(gx, gy); P = np.c_[GX.ravel(), GY.ravel()]
    it = (P[:, 0] > 1e-9) & (P[:, 0] < LX-1e-9) & (P[:, 1] > 1e-9) & (P[:, 1] < LY-1e-9)
    P[it] += h*0.3*rng.standard_normal(P[it].shape)
    # densify along the canal line
    yc = np.linspace(RECEPTOR[1], PORE[1], 40)
    P = np.vstack([P, np.c_[np.full(40, PORE[0]), yc],
                   np.c_[np.full(40, PORE[0])-CANAL_W, yc], np.c_[np.full(40, PORE[0])+CANAL_W, yc]])
    P = np.unique(np.round(P, 6), axis=0)
    m = TriMesh(P, Delaunay(P).simplices); cen = m.points[m.tris].mean(1)
    sig = np.where(cen[:, 1] < BODY_TOP, SIG_TIS, SIG_SW)        # tissue body / seawater
    if with_canal:
        incanal = (np.abs(cen[:, 0]-PORE[0]) < CANAL_W) & (cen[:, 1] > RECEPTOR[1]) & (cen[:, 1] < PORE[1]+1e-9)
        sig = np.where(incanal, SIG_GEL, sig)
    return m, sig


def receptor_voltage(m, sig, theta, E0=1.0):
    """Uniform far field E0 at angle theta: phi = -E0 (x cos + y sin) on the boundary."""
    K = assemble_primal(m, sig)
    bnd = np.where((m.points[:, 0] < 1e-6) | (m.points[:, 0] > LX-1e-6) |
                   (m.points[:, 1] < 1e-6) | (m.points[:, 1] > LY-1e-6))[0]
    bc = {int(i): -E0*(m.points[i, 0]*np.cos(theta)+m.points[i, 1]*np.sin(theta)) for i in bnd}
    v = solve_dirichlet(K, bc)
    def nearest(pt): return int(np.argmin(np.hypot(m.points[:, 0]-pt[0], m.points[:, 1]-pt[1])))
    ref = nearest((0.020, 0.010))                  # deep body interior reference
    return v, v[nearest(RECEPTOR)] - v[ref], nearest(RECEPTOR), ref


def figure():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    m, sig = build(with_canal=True)
    v, Vr, ir, iref = receptor_voltage(m, sig, theta=np.pi/2, E0=1.0)
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    cf = ax[0].tricontourf(m.points[:, 0]*1e3, m.points[:, 1]*1e3, m.tris, v*1e3, levels=26, cmap="RdBu_r")
    ax[0].axhline(BODY_TOP*1e3, color="k", lw=1.2, ls="--")
    ax[0].text(LX*1e3*0.5, BODY_TOP*1e3+1, "seawater (sigma=4)", ha="center", fontsize=8)
    ax[0].text(LX*1e3*0.5, BODY_TOP*1e3-4, "body tissue (sigma=0.3)", ha="center", fontsize=8)
    ax[0].plot([PORE[0]*1e3, PORE[0]*1e3], [RECEPTOR[1]*1e3, PORE[1]*1e3], color="#3B6D11", lw=4, alpha=.7)
    ax[0].plot(PORE[0]*1e3, PORE[1]*1e3, "ko", ms=6); ax[0].text(PORE[0]*1e3+1, PORE[1]*1e3, "pore", fontsize=8)
    ax[0].plot(RECEPTOR[0]*1e3, RECEPTOR[1]*1e3, "r*", ms=12); ax[0].text(RECEPTOR[0]*1e3+1, RECEPTOR[1]*1e3, "receptor", fontsize=8, color="r")
    ax[0].set_aspect("equal"); ax[0].set_xlabel("x (mm)"); ax[0].set_ylabel("y (mm)")
    ax[0].set_title("(a) Potential under a vertical field: the gel canal\ncarries the surface potential to the receptor")
    fig.colorbar(cf, ax=ax[0], shrink=.8, label="phi (mV per V/m)")
    th = np.linspace(0, 2*np.pi, 49)
    Vt = np.array([receptor_voltage(m, sig, theta=t, E0=1.0)[1] for t in th])
    axp = fig.add_subplot(1, 2, 2, projection="polar")
    axp.plot(th, np.abs(Vt)*1e3, color="#A32D2D", lw=2)
    axp.set_title("(b) Directional tuning |Vr| (mV per V/m)\nvs external field angle", pad=18)
    fig.suptitle("Ampulla of Lorenzini as a certified conduction problem: funnelling and directionality (dual DGM)", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "Fig_lorenzini.png"), dpi=125)
    print(f"-> saved  funnelling Vr(along)={Vr*1e3:.1f} mV/(V/m)")


if __name__ == "__main__":
    figure()
    E0 = 1.0   # 1 V/m far field (then scale to the nV/cm regime)
    m_c, sig_c = build(with_canal=True)
    m_n, sig_n = build(with_canal=False)
    _, Vr_canal, _, _ = receptor_voltage(m_c, sig_c, theta=np.pi/2, E0=E0)   # field along the canal (vertical)
    _, Vr_none, _, _ = receptor_voltage(m_n, sig_n, theta=np.pi/2, E0=E0)
    print(f"Lorenzini ampulla PoC (sigma sw/tis/gel = {SIG_SW}/{SIG_TIS}/{SIG_GEL} S/m)")
    print(f"  receptor voltage, field along canal: with canal = {Vr_canal*1e3:.3f} mV/(V/m), "
          f"without = {Vr_none*1e3:.3f} mV/(V/m)  -> funnelling x{Vr_canal/Vr_none:.1f}")
    # directional tuning
    print("  directional tuning |Vr(theta)| (mV per V/m):")
    for th in np.linspace(0, np.pi, 7):
        _, Vr, _, _ = receptor_voltage(m_c, sig_c, theta=th, E0=E0)
        print(f"    theta={np.degrees(th):3.0f} deg: {abs(Vr)*1e3:.3f}")
