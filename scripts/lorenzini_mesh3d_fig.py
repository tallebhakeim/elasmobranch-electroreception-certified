"""Figure: the actual 3D tetrahedral DGM mesh on the shark, with the ampulla canals resolved.

Panel (a): the shark surface immersed in seawater, with the snout ampulla array (pores and gel
canals). Panel (b): a mid-plane slab of the real tetrahedral mesh, coloured by material
(tissue / seawater / gel canal), showing the discretisation that carries the certified solve.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import lorenzini_shark as LS

SIG_TIS, SIG_SW, SIG_GEL = LS.SIG_TIS, LS.SIG_SW, LS.SIG_GEL
mesh = LS.load_shark()
print("building viz mesh ...", flush=True)
m, sig, amp, box = LS.build(mesh, h=0.020)
print(f"  np={m.np} nt={m.nt}", flush=True)
P = m.points; T = m.tets; cen = P[T].mean(1)

fig = plt.figure(figsize=(14, 6.2))

# (a) shark surface + ampulla array
ax = fig.add_subplot(1, 2, 1, projection="3d")
faces = mesh.vertices[mesh.faces]
pc = Poly3DCollection(faces, alpha=0.18, facecolor="#9aa7b0", edgecolor="none")
ax.add_collection3d(pc)
for ph, pore, recv in amp:
    ax.plot(*np.c_[pore, recv], color="#A32D2D", lw=2.5)
    ax.scatter(*pore, s=40, c="#A32D2D", depthshade=False)
b = mesh.bounds
ax.set_xlim(b[0, 0]-0.02, b[1, 0]+0.02); ax.set_ylim(b[0, 1]-0.02, b[1, 1]+0.02)
ax.set_zlim(b[0, 2]-0.02, b[1, 2]+0.02)
try: ax.set_box_aspect((b[1, 0]-b[0, 0], b[1, 1]-b[0, 1], b[1, 2]-b[0, 2]))
except Exception: pass
ax.view_init(elev=22, azim=-65); ax.set_axis_off()
ax.set_title("(a) Shark immersed in seawater\nwith the snout ampulla array (pores and gel canals)", fontsize=10)

# (b) real tetrahedral mesh, mid-plane slab, coloured by material
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
slab = (np.abs(cen[:, 2]) < 0.018)
# gel canal and seawater share sigma=4, so split the canal by geometry (distance to a canal axis)
canal = np.zeros(m.nt, bool)
for _, pore, recv in amp:
    dseg = recv-pore; L2 = float(dseg@dseg)
    t = np.clip(((cen-pore)@dseg)/L2, 0, 1); proj = pore + t[:, None]*dseg
    canal |= np.linalg.norm(cen-proj, axis=1) < LS.R_CAN
mat = np.full(m.nt, 0)                       # default seawater
mat[np.isclose(sig, SIG_TIS)] = 1            # tissue (uniquely sigma=0.3)
mat[canal] = 2                               # gel canal
fc = {0: "#cfe0ea", 1: "#c2a878", 2: "#d6452f"}
F = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
polys, cols = [], []
for ti in np.where(slab)[0]:
    verts = P[T[ti]]
    for f in F:
        polys.append(verts[f]); cols.append(fc[mat[ti]])
pc2 = Poly3DCollection(polys, facecolor=cols, edgecolor="0.35", linewidths=0.2, alpha=0.92)
ax2.add_collection3d(pc2)
ax2.set_xlim(b[0, 0]-0.05, b[1, 0]+0.05); ax2.set_ylim(b[0, 1]-0.05, b[1, 1]+0.05)
ax2.set_zlim(-0.1, 0.1)
ax2.view_init(elev=88, azim=-90); ax2.set_axis_off()
from matplotlib.patches import Patch

# Figure output directory. Defaults to the article folder, overridable so the repository
# runs anywhere (a hard-coded absolute path is not reproducible for a third party).
FIGDIR = os.environ.get("LORENZINI_FIGDIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures"))
os.makedirs(FIGDIR, exist_ok=True)

ax2.legend(handles=[Patch(facecolor=fc[1], label="tissue (0.3 S/m)"),
                    Patch(facecolor=fc[0], label="seawater (4 S/m)"),
                    Patch(facecolor=fc[2], label="gel canal (4 S/m)")],
           fontsize=8, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
ax2.set_title("(b) The real tetrahedral DGM mesh (mid-plane slab)\nresolving tissue, seawater and the gel canals", fontsize=10)

fig.suptitle("Three-dimensional tetrahedral discretisation behind the certified electroreception solve", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(FIGDIR, "Fig_mesh3d.png"), dpi=130)
print("-> saved Fig_mesh3d.png")
