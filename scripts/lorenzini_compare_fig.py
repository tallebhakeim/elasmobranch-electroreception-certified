"""Figure: ray vs shark electroreception comparison (reads /tmp/lorenzini_compare.npz)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.patches as mp
import lorenzini_compare as LC

# Figure output directory. Defaults to the article folder, overridable so the repository
# runs anywhere (a hard-coded absolute path is not reproducible for a third party).
FIGDIR = os.environ.get("LORENZINI_FIGDIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures"))
os.makedirs(FIGDIR, exist_ok=True)


d = np.load("/tmp/lorenzini_compare.npz", allow_pickle=True)
def arr(k): return np.array(d[k].tolist(), float)
data = {k: dict(VrB=arr(f"{k}_VrB"), VrE=arr(f"{k}_VrE")) for k in ("shark", "ray")}

fig, ax = plt.subplots(2, 2, figsize=(13.5, 9))
cols = {"prey": "#A32D2D", "nav": "#185FA5"}
for col, kind in enumerate(("shark", "ray")):
    VrE = np.abs(data[kind]["VrE"])*1e9; VrB = np.abs(data[kind]["VrB"])*1e9
    n = len(VrE); x = np.arange(n); w = 0.38
    a = ax[0, col]
    a.bar(x-w/2, VrE, w, color=cols["prey"], label="prey (dipole)")
    a.bar(x+w/2, VrB, w, color=cols["nav"], label="navigation (v x B)")
    a.set_yscale("log"); a.set_xlabel("ampulla"); a.set_ylabel("|receptor voltage| (nV)")
    a.set_title(f"({'a' if col==0 else 'b'}) {kind.upper()}: per-ampulla levels\n"
                f"prey max {VrE.max():.0f} nV (contrast {VrE.max()/max(VrE.min(),1e-9):.0f}x), "
                f"navig max {VrB.max():.0f} nV")
    a.legend(fontsize=8); a.grid(alpha=.3, axis="y")
    # geometry schematic (side view x-z)
    g = ax[1, col]; abc = LC.BODIES[kind]["abc"]; amp = LC.ampullae(kind, abc)
    th = np.linspace(0, 2*np.pi, 100)
    g.plot(abc[0]*np.cos(th)*100, abc[2]*np.sin(th)*100, "k", lw=1.5)   # body outline (x-z)
    g.fill(abc[0]*np.cos(th)*100, abc[2]*np.sin(th)*100, color="#cdbb9a", alpha=.4)
    vmax = max(VrE.max(), 1e-6)
    for (pore, recv), ve in zip(amp, VrE):
        g.plot([pore[0]*100, recv[0]*100], [pore[2]*100, recv[2]*100], color="0.4", lw=2)
        g.scatter(pore[0]*100, pore[2]*100, s=30+200*ve/vmax, c=[cols["prey"]], zorder=5)
    if kind == "shark":
        g.scatter((abc[0]+0.05)*100, 0, marker="*", s=160, c="#3B6D11"); g.text((abc[0]+0.05)*100, 1, "prey", color="#3B6D11", fontsize=8)
    else:
        g.scatter(0, (-abc[2]-0.05)*100, marker="*", s=160, c="#3B6D11"); g.text(1, (-abc[2]-0.05)*100, "prey", color="#3B6D11", fontsize=8)
    g.set_aspect("equal"); g.set_xlabel("x (cm)"); g.set_ylabel("z (cm)")
    g.set_title(f"({'c' if col==0 else 'd'}) {kind.upper()} body and ampullae\n"
                f"(pore size = prey response); {'fusiform, snout, detects ahead' if kind=='shark' else 'flat disc, ventral, detects below'}")
fig.suptitle("Ray vs shark electroreception: morphology sets the detection geometry (certified DGM conduction)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(os.path.join(FIGDIR, "Fig_ray_vs_shark.png"), dpi=125)
print("-> saved Fig_ray_vs_shark.png")
