"""Real anatomical meshes (CC BY-NC, DigitalLife3D / Jer Bot) for the electroreception study:
a great hammerhead shark (Sphyrna mokarran) and a manta ray (Mobula birostris).

Both raw meshes share an orientation: long axis = z (anterior at +z), dorsal = +y, lateral = x.
We reorient to a body frame X=anterior(+ forward), Y=lateral, Z=dorsoventral(+ dorsal), scale to a
head-sized model, and place ampullae anatomically: across the ventral cephalofoil for the
hammerhead (the wide 'hammer' that spreads the electrosensory array), and ventral around the mouth
for the manta. The certified DGM conduction engine then gives the electric (prey dipole) and
magnetic (motional v x B) receptor levels, on the real geometry.
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import trimesh

# Surface models are looked up in the article folder first, so the results stay reproducible
# once the meshes are archived next to the manuscript, and only then in ~/Downloads. Set
# LORENZINI_MESHES to override. Source: DigitalLife3D (digitallife3d.org), CC BY-NC 4.0.
MESH_DIRS = [
    os.environ.get("LORENZINI_MESHES", ""),
    # meshes/ next to this script's parent: this is where they sit in the public repository,
    # so a reviewer who clones it finds them without setting anything.
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "meshes"),
    os.path.expanduser("~/Downloads"),
]
PATTERNS = {"hammerhead": ["*hammerhead*.obj", "*hammerhead*/*.gltf", "*hammerhead*/*.glb"],
            "manta": ["*manta*.obj", "*manta*/*.gltf", "*manta*/*.glb", "*manta*.glb"]}


def find_mesh(species):
    for d in MESH_DIRS:
        if not d or not os.path.isdir(d):
            continue
        for pat in PATTERNS[species]:
            hits = sorted(glob.glob(os.path.join(d, pat)))
            if hits:
                return hits[0]
    raise FileNotFoundError(
        f"no {species} surface model found. Looked in: "
        + ", ".join(d for d in MESH_DIRS if d)
        + ". Download it from DigitalLife3D (digitallife3d.org, CC BY-NC 4.0) and put it in "
          "the article folder, or point LORENZINI_MESHES at it.")


HAMMER = None           # resolved lazily by load(), so importing this file never fails
MANTA = None
L_TARGET = 1.0          # scale longest body extent to this (m); head features then ~ dm


def reorient(m):
    """Permute raw (x,y,z) -> body frame (X=anterior, Y=lateral, Z=dorsoventral)."""
    v = m.vertices.copy()
    V = np.c_[v[:, 2], v[:, 0], v[:, 1]]      # X=old z (anterior), Y=old x (lateral), Z=old y (dorsal)
    out = trimesh.Trimesh(vertices=V, faces=m.faces, process=False)
    s = L_TARGET / out.extents.max()
    out.apply_scale(s)
    out.apply_translation(-out.bounds.mean(0))   # centre
    return out


def load(species):
    f = find_mesh(species)
    print(f"  [{species}] {f}")
    return reorient(trimesh.load(f, force="mesh"))


CANAL_LEN = 0.06*L_TARGET


def ampullae(species, mesh):
    """Anatomical ampullae spread laterally across the front-ventral skin. For each target (X,Y) we
    read the real ventral surface height (the lowest mesh vertices in that column) as the pore, and
    run the canal straight up into the body (the receptor). Clean lateral spread, real ventral skin."""
    b = mesh.bounds; Xf = b[1, 0]; V = mesh.vertices; n = 9
    if species == "hammerhead":
        span = 0.42*(b[1, 1]-b[0, 1]); rad = 0.05*L_TARGET
        targ = [(Xf - 0.08*(b[1, 0]-b[0, 0]) - 0.12*abs(y), y) for y in np.linspace(-span, span, n)]
    else:
        span = 0.10*(b[1, 1]-b[0, 1]); rad = 0.04*L_TARGET
        targ = [(Xf - 0.05*(b[1, 0]-b[0, 0]), y) for y in np.linspace(-span, span, n)]
    out = []
    for (xx, yy) in targ:
        near = (np.abs(V[:, 0]-xx) < rad) & (np.abs(V[:, 1]-yy) < rad)
        zv = np.percentile(V[near, 2], 8) if near.sum() > 3 else b[0, 2]   # ventral skin height
        pore = np.array([xx, yy, zv])
        recv = pore + np.array([0, 0, CANAL_LEN])     # canal runs up into the body
        out.append((yy, pore, recv))
    return out


def head_box(species, mesh):
    """A sub-box around the head, where we build the fine conduction mesh (keeps np bounded)."""
    b = mesh.bounds
    amp = ampullae(species, mesh)
    pores = np.array([p for _, p, _ in amp] + [r for _, _, r in amp])
    c = pores.mean(0)
    half = np.array([0.22, 0.30, 0.18])*L_TARGET
    if species == "manta":
        half = np.array([0.18, 0.18, 0.14])*L_TARGET
    return c, half, amp


import lorenzini_shark as LS
from scipy.spatial import Delaunay
from dgm.mesh3d import TetMesh

# Figure output directory. Defaults to the article folder, overridable so the repository
# runs anywhere (a hard-coded absolute path is not reproducible for a third party).
FIGDIR = os.environ.get("LORENZINI_FIGDIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures"))
os.makedirs(FIGDIR, exist_ok=True)


SIG_SW, SIG_TIS, SIG_GEL = LS.SIG_SW, LS.SIG_TIS, LS.SIG_GEL
R_CAN = 0.020*L_TARGET
THRESH = LS.THRESH; B_EARTH = LS.B_EARTH; V_SWIM = LS.V_SWIM
PREY_DIPOLE = LS.PREY_DIPOLE


def build_real(species, mesh, amp, h=0.020, seed=0):
    rng = np.random.default_rng(seed)
    c, half, _ = head_box(species, mesh)
    lo, hi = c-half, c+half
    grid = [np.arange(lo[i], hi[i]+h, h) for i in range(3)]
    GX, GY, GZ = np.meshgrid(*grid); P = np.c_[GX.ravel(), GY.ravel(), GZ.ravel()]
    inbox = np.all((P > lo+1e-9) & (P < hi-1e-9), axis=1)    # keep face nodes exactly on the box faces
    P[inbox] += h*0.25*rng.standard_normal(P[inbox].shape)
    extra = []
    for _, pore, recv in amp:
        for s in np.linspace(0, 1, 12):
            extra.append(pore + s*(recv-pore) + rng.normal(0, R_CAN*0.4, 3))
    P = np.unique(np.round(np.vstack([P, np.array(extra)]), 5), axis=0)
    in_node = mesh.contains(P)
    m = TetMesh(P, Delaunay(P).simplices)
    inside = in_node[m.tets].mean(1) >= 0.5
    cen = m.points[m.tets].mean(1)
    sig = np.where(inside, SIG_TIS, SIG_SW)
    for _, pore, recv in amp:
        dd = recv-pore; L2 = float(dd@dd)
        t = np.clip(((cen-pore)@dd)/L2, 0, 1); proj = pore + t[:, None]*dd
        sig = np.where(np.linalg.norm(cen-proj, axis=1) < R_CAN, SIG_GEL, sig)
    return m, sig, (lo, hi), inside


def solve_species(species, h=0.020):
    mesh = load(species); amp = ampullae(species, mesh)
    m, sig, box, inside = build_real(species, mesh, amp, h=h)
    print(f"  {species}: np={m.np} nt={m.nt}", flush=True)
    K, bnd, free, lu = LS.factor(m, sig, box)
    # reference: a deep tissue node near the body interior
    incen = m.points[m.tets][inside].mean(1).mean(0)
    ref = int(np.argmin(np.sum((m.points-incen)**2, axis=1)))

    def vr(phi): return np.array([phi[LS.nidx(m, r)] - phi[ref] for _, _, r in amp])
    Emot = V_SWIM*B_EARTH
    VrB = vr(LS.solve_uniform(m, K, bnd, free, lu, np.array([0.0, Emot, 0.0])))
    # prey: a current dipole in the water just in front of (and below) the ampulla array
    ac = np.array([p for _, p, _ in amp]).mean(0)
    prey = ac + np.array([0.08*L_TARGET, 0, -0.04*L_TARGET])
    dd = 0.05*L_TARGET; I = PREY_DIPOLE/dd
    VrE = vr(LS.solve_dipole(m, K, bnd, free, lu, prey+[0, dd/2, 0], prey-[0, dd/2, 0], I))
    pores = np.array([p for _, p, _ in amp])
    return dict(species=species, ys=np.array([a[0] for a in amp]), VrB=VrB, VrE=VrE, pores=pores,
                surf=mesh, m=m, sig=sig, amp=amp, Emot=Emot)


ART = "/Users/talleb/Desktop/à faire/article/article 7"
COL = {"nav": "#185FA5", "prey": "#A32D2D"}


def fig_eb(d):  # Figure 5: E+B on the real hammerhead
    import matplotlib.pyplot as plt
    VrB, VrE, y = np.abs(d["VrB"])*1e9, np.abs(d["VrE"])*1e9, d["ys"]
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    x = np.arange(len(y)); w = 0.38
    ax[0].bar(x-w/2, VrB, w, color=COL["nav"], label="B: navigation (v x B, uniform)")
    ax[0].bar(x+w/2, VrE, w, color=COL["prey"], label="E: prey dipole (localised)")
    ax[0].axhline(THRESH*0.05*1e9, color="k", ls="--", lw=1, label="~receptor floor")
    ax[0].set_yscale("log"); ax[0].set_xlabel("ampulla (lateral position on the cephalofoil)")
    ax[0].set_ylabel("|receptor voltage| (nV)")
    ax[0].set_title("(a) Levels per ampulla on the hammerhead cephalofoil")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, axis="y")
    ax[1].plot(y, VrB, "o-", color=COL["nav"], label="navigation (uniform)")
    ax[1].plot(y, VrE, "s-", color=COL["prey"], label="prey (localised)")
    ax[1].set_xlabel("ampulla lateral position Y (m)"); ax[1].set_ylabel("|receptor voltage| (nV)")
    ax[1].set_title("(b) The wide cephalofoil spreads the array:\nthe uniform field gives a V-shaped baseline, the prey a local peak")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    fig.suptitle("Electric and magnetic electroreception on a real great hammerhead (certified DGM conduction)", fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "Fig_shark_EB.png"), dpi=125); plt.close(fig)
    print("-> Fig_shark_EB.png (hammerhead)")


def fig_compare(R):  # Figure 6: real hammerhead vs real manta
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(13.5, 9))
    for col, sp in enumerate(("hammerhead", "manta")):
        d = R[sp]; VrB, VrE = np.abs(d["VrB"])*1e9, np.abs(d["VrE"])*1e9; n = len(VrB)
        x = np.arange(n); w = 0.38
        a = ax[0, col]
        a.bar(x-w/2, VrE, w, color=COL["prey"], label="prey (dipole)")
        a.bar(x+w/2, VrB, w, color=COL["nav"], label="navigation (v x B)")
        a.set_yscale("log"); a.set_xlabel("ampulla"); a.set_ylabel("|receptor voltage| (nV)")
        a.set_title(f"({'a' if col==0 else 'b'}) {sp.upper()}: per-ampulla levels\n"
                    f"prey max {VrE.max():.0f} nV, navigation max {VrB.max():.0f} nV")
        a.legend(fontsize=8); a.grid(alpha=.3, axis="y")
        g = ax[1, col]; V = d["surf"].vertices
        g.scatter(V[::6, 0], V[::6, 1], s=0.4, c="#b9c4cc", alpha=.4)   # ventral (X-Y) outline
        vmax = max(VrE.max(), 1e-6)
        for (xx, yy), ve in zip(d["pores"][:, :2], VrE):
            g.scatter(xx, yy, s=30+260*ve/vmax, c=COL["prey"], zorder=5, edgecolor="k", linewidths=.3)
        g.set_aspect("equal"); g.set_xlabel("X anterior (m)"); g.set_ylabel("Y lateral (m)")
        g.set_title(f"({'c' if col==0 else 'd'}) {sp.upper()} ventral view, ampullae sized by prey response\n"
                    + ("wide cephalofoil spreads the array" if sp == "hammerhead" else "ventral cluster around the mouth"))
    fig.suptitle("Hammerhead vs manta on real anatomical meshes: the body plan sets the detection geometry (certified DGM)", fontsize=12.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "Fig_ray_vs_shark.png"), dpi=125); plt.close(fig)
    print("-> Fig_ray_vs_shark.png (hammerhead vs manta)")


def fig_mesh(d):  # Figure 2: the real tetrahedral mesh on the hammerhead
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    surf, m, sig, amp = d["surf"], d["m"], d["sig"], d["amp"]
    P, T = m.points, m.tets; cen = P[T].mean(1)
    fig = plt.figure(figsize=(14, 6.2))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.add_collection3d(Poly3DCollection(surf.vertices[surf.faces], alpha=0.16, facecolor="#9aa7b0", edgecolor="none"))
    for _, pore, recv in amp:
        ax.plot(*np.c_[pore, recv], color="#A32D2D", lw=2.2); ax.scatter(*pore, s=30, c="#A32D2D", depthshade=False)
    b = surf.bounds
    for i, s in enumerate("xyz"):
        getattr(ax, f"set_{s}lim")(b[0, i], b[1, i])
    try: ax.set_box_aspect(b[1]-b[0])
    except Exception: pass
    ax.view_init(elev=70, azim=-90); ax.set_axis_off()
    ax.set_title("(a) Great hammerhead immersed in seawater\nventral ampulla array on the cephalofoil", fontsize=10)
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    canal = np.zeros(m.nt, bool)
    for _, pore, recv in amp:
        dd = recv-pore; L2 = float(dd@dd); t = np.clip(((cen-pore)@dd)/L2, 0, 1)
        canal |= np.linalg.norm(cen-(pore+t[:, None]*dd), axis=1) < R_CAN
    mat = np.full(m.nt, 0); mat[np.isclose(sig, SIG_TIS)] = 1; mat[canal] = 2
    fc = {0: "#cfe0ea", 1: "#c2a878", 2: "#d6452f"}
    zc = np.median(cen[:, 2]); slab = np.abs(cen[:, 2]-zc) < 0.02*L_TARGET
    F = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]); polys, cols = [], []
    for ti in np.where(slab)[0]:
        for f in F: polys.append(P[T[ti]][f]); cols.append(fc[mat[ti]])
    ax2.add_collection3d(Poly3DCollection(polys, facecolor=cols, edgecolor="0.35", linewidths=0.2, alpha=0.92))
    for i in range(2):
        getattr(ax2, f"set_{'xyz'[i]}lim")(b[0, i], b[1, i])
    ax2.set_zlim(zc-0.1, zc+0.1); ax2.view_init(elev=88, azim=-90); ax2.set_axis_off()
    from matplotlib.patches import Patch
    ax2.legend(handles=[Patch(facecolor=fc[1], label="tissue (0.3 S/m)"), Patch(facecolor=fc[0], label="seawater (4 S/m)"),
                        Patch(facecolor=fc[2], label="gel canal (4 S/m)")], fontsize=8, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
    ax2.set_title("(b) The real tetrahedral DGM mesh (horizontal slab)\nresolving tissue, seawater and the gel canals", fontsize=10)
    fig.suptitle("Three-dimensional tetrahedral discretisation on a real great hammerhead head", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(os.path.join(FIGDIR, "Fig_mesh3d.png"), dpi=130); plt.close(fig)
    print("-> Fig_mesh3d.png (hammerhead)")


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg")
    R = {sp: solve_species(sp) for sp in ("hammerhead", "manta")}
    for sp, d in R.items():
        print(f"{sp}: nav {d['Emot']*1e6:.0f} uV/m -> |Vr| {np.round(np.abs(d['VrB'])*1e9,0)} nV ; "
              f"prey max {np.abs(d['VrE']).max()*1e9:.0f} nV")
    fig_mesh(R["hammerhead"]); fig_eb(R["hammerhead"]); fig_compare(R)
