# -*- coding: utf-8 -*-
"""Regenerate every figure and every number of the article from scratch.

This is the entry point a reviewer should run. It executes each script in order, reports
what passed, what was skipped and why, and writes reproduce_all.log.

    python3 reproduce_all.py

Scripts that need the two surface models (DigitalLife3D, CC BY-NC 4.0) are skipped with an
explicit message if those files are absent, rather than failing obscurely. Put them in
meshes/ or point LORENZINI_MESHES at them.
"""
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
MEGA_EX = os.path.join(HERE, "scripts")
LOG = os.path.join(HERE, "reproduce_all.log")

# (script path, working directory, what it produces, needs the surface models?)
STEPS = [
    (os.path.join(MEGA_EX, "lorenzini_certified.py"), MEGA_EX,
     "guaranteed brackets on the access conductance", False),
    (os.path.join(MEGA_EX, "lorenzini_poc.py"), MEGA_EX,
     "canal funnelling, directional sweep, Fig_lorenzini.png", False),
    (os.path.join(MEGA_EX, "lorenzini_freq.py"), MEGA_EX,
     "band-pass corners, Fig_lorenzini_freq.png", False),
    (os.path.join(MEGA_EX, "lorenzini_3d.py"), MEGA_EX,
     "directional rosette, Fig_lorenzini_3d.png", False),
    (os.path.join(MEGA_EX, "lorenzini_cloak.py"), MEGA_EX,
     "electrical visibility sweep", False),
    (os.path.join(MEGA_EX, "lorenzini_cloak_fig.py"), MEGA_EX,
     "Fig_cloak.png", False),
    # Fig_mesh3d.png used to come from lorenzini_mesh3d_fig.py, which reads an old shark STL
    # that is not part of this article. lorenzini_real.py now draws it from the DigitalLife3D
    # surface model, so the figure and the receptor levels come from the same mesh.
    (os.path.join(MEGA_EX, "lorenzini_real.py"), MEGA_EX,
     "receptor levels on both body plans, Fig_mesh3d.png, Fig_shark_EB.png, Fig_ray_vs_shark.png",
     True),
    ("__sync_figures__", HERE,
     "copy regenerated PNGs from the engine directory into figures/", False),
    (os.path.join(HERE, "sensitivity_lorenzini.py"), HERE,
     "sensitivity campaign, sensitivity_results.json", False),
    (os.path.join(HERE, "make_sensitivity_figures.py"), HERE,
     "Fig_S1 to Fig_S4", False),
    (os.path.join(HERE, "build_si.py"), HERE,
     "Supplementary Information", False),
    (os.path.join(HERE, "build_jtb.py"), HERE,
     "manuscript, highlights, cover letter", "figures"),
]


def meshes_present():
    sys.path.insert(0, MEGA_EX)
    try:
        import lorenzini_real as LR
        for sp in ("hammerhead", "manta"):
            LR.find_mesh(sp)
        return True, ""
    except FileNotFoundError as e:
        return False, str(e).split(". Download")[0]
    except Exception as e:                       # import-time problem, report it as-is
        return False, f"{type(e).__name__}: {e}"


def sync_figures():
    """The engine scripts write their PNGs next to themselves; the manuscript reads figures/.
    Copy across so a full run really does rebuild the document from regenerated figures."""
    dst = os.path.join(HERE, "figures")
    os.makedirs(dst, exist_ok=True)
    n = 0
    for f in sorted(os.listdir(MEGA_EX)):
        if f.startswith("Fig_") and f.endswith(".png"):
            shutil.copy2(os.path.join(MEGA_EX, f), os.path.join(dst, f))
            n += 1
    return n


def main():
    have_meshes, why = meshes_present()
    lines, npass, nskip, nfail = [], 0, 0, 0
    t0 = time.time()
    print(f"reproducing from {HERE}\n")
    if not have_meshes:
        print(f"NOTE: surface models absent ({why}); the two body-plan steps will be skipped.\n")

    mesh_figs = ["Fig_mesh3d.png", "Fig_shark_EB.png", "Fig_ray_vs_shark.png"]

    def figs_present():
        d = os.path.join(HERE, "figures")
        return all(os.path.exists(os.path.join(d, f)) for f in mesh_figs)

    for path, cwd, what, needs_mesh in STEPS:
        name = os.path.basename(path)
        # "figures": needs the mesh-derived figures on disk, which a previous run may have
        # left there even when the surface models themselves are gone.
        if needs_mesh == "figures":
            if not figs_present():
                print(f"SKIP {name:32s} {what} (Figures 1, 5, 6 absent)")
                lines.append(f"SKIP {name}: mesh-derived figures absent")
                nskip += 1
                continue
        elif needs_mesh and not have_meshes:
            print(f"SKIP {name:32s} {what}")
            lines.append(f"SKIP {name}: surface models absent")
            nskip += 1
            continue
        if path == "__sync_figures__":
            n = sync_figures()
            print(f"ok   {'sync figures':32s} {n} PNG copied into figures/")
            lines.append(f"ok   sync figures: {n} copied")
            npass += 1
            continue
        t = time.time()
        env = dict(os.environ, LORENZINI_FIGDIR=os.path.join(HERE, "figures"))
        r = subprocess.run([sys.executable, path], cwd=cwd, env=env,
                           capture_output=True, text=True)
        dt = time.time() - t
        if r.returncode == 0:
            print(f"ok   {name:32s} {what}  ({dt:.1f}s)")
            lines.append(f"ok   {name} ({dt:.1f}s)\n{r.stdout}")
            npass += 1
        else:
            print(f"FAIL {name:32s} exit {r.returncode}  ({dt:.1f}s)")
            print((r.stderr or "").strip()[-700:])
            lines.append(f"FAIL {name} exit {r.returncode}\n{r.stdout}\n{r.stderr}")
            nfail += 1

    summary = (f"\n{npass} ok, {nskip} skipped, {nfail} failed "
               f"in {time.time()-t0:.0f} s")
    print(summary)
    if nskip:
        print("Skipped steps need the DigitalLife3D surface models (CC BY-NC 4.0). "
              "Put them in meshes/ or set LORENZINI_MESHES.")
    with open(LOG, "w") as f:
        f.write("\n".join(lines) + summary + "\n")
    print(f"log: {os.path.basename(LOG)}")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
