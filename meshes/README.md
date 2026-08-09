# Surface models

Two photogrammetric surface models are needed by `scripts/lorenzini_real.py` and
`scripts/lorenzini_mesh3d_fig.py`:

- a great hammerhead shark, *Sphyrna mokarran*
- a manta ray, *Mobula birostris*

Both are by DigitalLife3D (Jer Bot, digitallife3d.org) and distributed under the Creative
Commons Attribution-NonCommercial 4.0 licence (CC BY-NC 4.0). Download them and place them in
this directory.

Files are located by pattern, so exact filenames do not matter: anything matching
`*hammerhead*` or `*manta*` with extension `.obj`, `.gltf` or `.glb` is found. Set the
environment variable `LORENZINI_MESHES` to use a different directory.

The scripts reorient each model to a body frame (X anterior, Y lateral, Z dorsoventral),
scale the longest extent to one metre, and centre it; that transformation is in
`lorenzini_real.py` and is applied to the distributed models unchanged.
