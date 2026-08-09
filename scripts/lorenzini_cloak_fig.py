"""Figure: electroquasistatic cloak that neutralises the electric signature (reads npz)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import lorenzini_cloak as LC

# Figure output directory. Defaults to the article folder, overridable so the repository
# runs anywhere (a hard-coded absolute path is not reproducible for a third party).
FIGDIR = os.environ.get("LORENZINI_FIGDIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures"))
os.makedirs(FIGDIR, exist_ok=True)


d = np.load("/tmp/lorenzini_cloak.npz", allow_pickle=True)
pts, tris, region = d["pts"], d["tris"], d["region"]
names = list(d["names"]); pmags = d["pmags"]; p_bare = float(d["p_bare"])
sweep, pcurve, snull, s_neut = d["sweep"], d["pcurve"], float(d["snull"]), float(d["s_neutral"])
tri = mtri.Triangulation(pts[:, 0], pts[:, 1], tris)
A, B, T = LC.A_BODY, LC.B_BODY, LC.T_SUIT
th = np.linspace(0, 2*np.pi, 200)


def outline(ax):
    ax.plot(A*np.cos(th), B*np.sin(th), "k", lw=1.4)
    ax.plot((A+T)*np.cos(th), (B+T)*np.sin(th), color="0.5", lw=1, ls="--")


PMAX = None
def fieldmap(ax, phi, title):
    global PMAX
    pert = (phi - (-LC.E0*pts[:, 0]))*1e3      # perturbation potential (mV), removes the uniform field
    if PMAX is None: PMAX = np.abs(pert).max()
    lv = np.linspace(-PMAX, PMAX, 25)
    cf = ax.tricontourf(tri, pert, levels=lv, cmap="RdBu_r", extend="both")
    ax.tricontour(tri, pert, levels=lv, colors="k", linewidths=0.25, alpha=.5)
    outline(ax)
    ax.scatter(*LC.PROBE, marker="v", s=70, c="#1a1a1a", zorder=6)
    ax.text(LC.PROBE[0], LC.PROBE[1]+0.05, "shark\nampulla", color="#1a1a1a", fontsize=7, ha="center")
    ax.set_xlim(-0.5, 0.5); ax.set_ylim(-0.4, 0.4); ax.set_aspect("equal")
    ax.set_title(title, fontsize=10); ax.set_xticks([]); ax.set_yticks([])
    return cf


fig = plt.figure(figsize=(13.5, 9))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])

# (a) bare body: perturbation potential is a clear dipole the shark reads
ax_a = fig.add_subplot(gs[0, 0])
fieldmap(ax_a, d["phi_0"], "(a) Bare body: perturbation is a dipole\n(the signature the shark reads)")
# (b) cloaked: perturbation collapses -> invisible
ax_b = fig.add_subplot(gs[0, 1])
cf = fieldmap(ax_b, d["phi_3"], "(b) Neutral-inclusion suit: perturbation collapses\n(electrically invisible)")
cb = fig.colorbar(cf, ax=[ax_a, ax_b], shrink=0.8, pad=0.02)
cb.set_label("perturbation potential (mV)", fontsize=8)

# (c) sweep: dipole source vs coating conductivity, with the cloaking null
ax_c = fig.add_subplot(gs[1, 0])
ax_c.semilogy(sweep, pcurve/p_bare, color="#444", lw=2)
ax_c.axhline(1.0, color="0.6", lw=.8, ls=":")
ax_c.scatter([snull], [pcurve.min()/p_bare], s=80, c="#185FA5", zorder=6,
             label=f"cloaking null sigma~{snull:.0f} S/m ({20*np.log10(pcurve.min()/p_bare):.0f} dB)")
ax_c.axvline(LC.SIG_SW, color="#2C7BB6", lw=.8, ls="--")
ax_c.text(LC.SIG_SW, 1.6, " seawater", color="#2C7BB6", fontsize=8, rotation=90, va="bottom")
ax_c.axvline(LC.SIG_TIS, color="#A32D2D", lw=.8, ls="--")
ax_c.text(LC.SIG_TIS, 1.6, " tissue", color="#A32D2D", fontsize=8, rotation=90, va="bottom")
ax_c.set_xscale("log"); ax_c.set_xlabel("suit conductivity (S/m)")
ax_c.set_ylabel("dipole signature / bare body")
ax_c.set_title("(c) Tuning the suit conductivity: a deep cloaking null\n(an insulating wetsuit, far left, is worse than bare skin)")
ax_c.grid(alpha=.3, which="both"); ax_c.legend(fontsize=8, loc="upper right")

# (d) suppression bar for the four cases
ax_d = fig.add_subplot(gs[1, 1])
sup = [20*np.log10(max(p, 1e-30)/p_bare) for p in pmags]
short = ["bare\nbody", "neoprene\nwetsuit", "seawater-\nmatched", "neutral\ncloak"]
cols = ["#888780", "#A32D2D", "#E08214", "#185FA5"]
b = ax_d.bar(short, sup, color=cols)
ax_d.axhline(0, color="k", lw=.8)
for bi, s in zip(b, sup):
    ax_d.text(bi.get_x()+bi.get_width()/2, s + (0.4 if s >= 0 else -0.8),
              f"{s:+.1f} dB", ha="center", fontsize=9, fontweight="bold")
ax_d.set_ylabel("signature vs bare body (dB)")
ax_d.set_title("(d) The insulating wetsuit increases the signature;\nonly a conductive, tuned suit suppresses it")
ax_d.grid(alpha=.3, axis="y")

fig.suptitle("Neutralising the electric signature in seawater: an electroquasistatic cloak (certified DGM conduction)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(os.path.join(FIGDIR, "Fig_cloak.png"), dpi=125)
print("-> saved Fig_cloak.png")
