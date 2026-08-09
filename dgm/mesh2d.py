"""2D primal/dual mesh container for the Discrete Geometric Method (DGM).

Implements the interlocked primal (Delaunay triangle) / dual (circumcenter,
i.e. Voronoi) meshes described in:

  Ren Dan, Xu Xiaoyu, Qu Hui, Ren Zhuoxiang,
  "Two-dimensional parasitic capacitance extraction for integrated circuit
   with dual discrete geometric methods", J. Semicond. 36(4), 045008 (2015).

Geometry conventions (paper Fig. 3-4):
  - For a triangle, the edge "opposite" a local vertex carries the primal
    capacitance of that triangle's contribution.
  - l_P : length of a primal edge.
  - l_D : signed distance from the triangle circumcenter to that edge,
          positive when the circumcenter lies on the same side as the
          opposite vertex (the "interior" side). For a well-centred
          (acute / Delaunay) mesh l_D >= 0.
  - The ratio l_D / l_P = 1/2 * cot(angle opposite the edge).

This module is deliberately self-contained and dimension-tagged ("2d") so a
sibling mesh3d.py (tetra primal / Voronoi dual) can mirror the same API later.
"""
from __future__ import annotations

import numpy as np


def _circumcenter(a, b, c):
    """Circumcenter of triangle (a, b, c). Each arg is an (...,2) array."""
    ax, ay = a[..., 0], a[..., 1]
    bx, by = b[..., 0], b[..., 1]
    cx, cy = c[..., 0], c[..., 1]
    d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    a2 = ax * ax + ay * ay
    b2 = bx * bx + by * by
    c2 = cx * cx + cy * cy
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    return np.stack([ux, uy], axis=-1)


class TriMesh:
    """Primal triangle mesh with derived dual (circumcenter) geometry.

    Parameters
    ----------
    points : (Np, 2) float array
        Node coordinates of the primal mesh.
    tris : (Nt, 3) int array
        Triangle connectivity (vertex indices).
    """

    def __init__(self, points, tris):
        self.points = np.asarray(points, dtype=float)
        self.tris = np.asarray(tris, dtype=int)
        self.np = len(self.points)
        self.nt = len(self.tris)
        self._build_geometry()
        self._build_edges()

    # ------------------------------------------------------------------ #
    # Geometry
    # ------------------------------------------------------------------ #
    def _build_geometry(self):
        p = self.points
        a = p[self.tris[:, 0]]
        b = p[self.tris[:, 1]]
        c = p[self.tris[:, 2]]
        # signed area (positive if CCW)
        self.area = 0.5 * ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
                           - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))
        self.circumcenter = _circumcenter(a, b, c)

    def _local_edge_geometry(self, t, i_local):
        """Return (l_P, l_D) for the edge of triangle t opposite local vertex
        i_local.  l_D is the SIGNED circumcenter-to-edge distance (interior +).
        """
        v = self.tris[t]
        opp = v[i_local]                       # opposite vertex
        e0 = v[(i_local + 1) % 3]
        e1 = v[(i_local + 2) % 3]
        p0, p1, pp = self.points[e0], self.points[e1], self.points[opp]
        edge = p1 - p0
        lp = np.hypot(*edge)
        # unit normal oriented toward the opposite vertex
        n = np.array([-edge[1], edge[0]])
        n /= np.hypot(*n)
        if np.dot(pp - p0, n) < 0:
            n = -n
        o = self.circumcenter[t]
        ld = np.dot(o - p0, n)                  # signed distance
        return lp, ld

    # ------------------------------------------------------------------ #
    # Edge topology
    # ------------------------------------------------------------------ #
    def _build_edges(self):
        """Build unique undirected edges and their incident triangles.

        edges          : (Ne, 2) int, sorted node pairs
        edge_tris      : list of [(tri, local_opposite_vertex), ...]
        boundary_edge  : (Ne,) bool, True if the edge has a single triangle
        """
        edict = {}
        edge_tris = []
        for t, v in enumerate(self.tris):
            for i_local in range(3):
                a = v[(i_local + 1) % 3]
                b = v[(i_local + 2) % 3]
                key = (a, b) if a < b else (b, a)
                if key not in edict:
                    edict[key] = len(edge_tris)
                    edge_tris.append([])
                edge_tris[edict[key]].append((t, i_local))
        self.edges = np.array(list(edict.keys()), dtype=int)
        self.ne = len(self.edges)
        self.edge_tris = edge_tris
        self.boundary_edge = np.array([len(et) == 1 for et in edge_tris])

    # ------------------------------------------------------------------ #
    # Boundary nodes
    # ------------------------------------------------------------------ #
    def boundary_nodes(self):
        nodes = set()
        for e, is_b in zip(self.edges, self.boundary_edge):
            if is_b:
                nodes.update(e.tolist())
        return np.array(sorted(nodes), dtype=int)


def annulus_mesh(r_in, r_out, n_in=40, n_out=120, n_radial=18, seed=0):
    """Unstructured Delaunay mesh of an annulus r_in < r < r_out.

    Boundary nodes lie exactly on the two circles (clean Dirichlet data).
    Triangles whose centroid falls in the hole are removed.

    Returns
    -------
    mesh : TriMesh
    inner_nodes, outer_nodes : index arrays of boundary nodes on each circle
    """
    rng = np.random.default_rng(seed)
    pts = []
    # inner circle
    th = np.linspace(0, 2 * np.pi, n_in, endpoint=False)
    pts.append(np.c_[r_in * np.cos(th), r_in * np.sin(th)])
    # outer circle
    th = np.linspace(0, 2 * np.pi, n_out, endpoint=False)
    pts.append(np.c_[r_out * np.cos(th), r_out * np.sin(th)])
    # interior jittered concentric rings
    for k in range(1, n_radial):
        r = r_in + (r_out - r_in) * k / n_radial
        m = int(n_in + (n_out - n_in) * k / n_radial)
        th = np.linspace(0, 2 * np.pi, m, endpoint=False) + rng.uniform(0, 2 * np.pi)
        jr = r + (r_out - r_in) / n_radial * 0.25 * rng.standard_normal(m)
        pts.append(np.c_[jr * np.cos(th), jr * np.sin(th)])
    points = np.vstack(pts)

    from scipy.spatial import Delaunay
    tri = Delaunay(points)
    # keep triangles whose centroid radius is inside the annulus
    cent = points[tri.simplices].mean(axis=1)
    rc = np.hypot(cent[:, 0], cent[:, 1])
    keep = (rc > r_in) & (rc < r_out)
    tris = tri.simplices[keep]

    # drop unused points and reindex
    used = np.unique(tris)
    remap = -np.ones(len(points), dtype=int)
    remap[used] = np.arange(len(used))
    points = points[used]
    tris = remap[tris]

    mesh = TriMesh(points, tris)
    rr = np.hypot(points[:, 0], points[:, 1])
    inner = np.where(np.abs(rr - r_in) < 1e-9)[0]
    outer = np.where(np.abs(rr - r_out) < 1e-9)[0]
    return mesh, inner, outer
