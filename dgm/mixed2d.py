"""2D complementary lower bound via RT0-P0 mixed FEM (robust on any triangulation).

The 2D counterpart of mixed3d: lowest-order Raviart-Thomas flux (one DoF per
edge) + piecewise-constant potential.  Unlike the circumcentric "diagonal-Hodge"
dual, this is a genuine H(div) discretization, so it stays valid under
aggressive (adaptive) refinement where slivers would break the circumcentric one.

RT0 basis (edge opposite vertex i):  phi_i(x) = (x - p_i)/(2A), unit flux through
edge i, div phi_i = 1/A.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _edges_and_signs(mesh):
    ge = np.zeros((mesh.nt, 3), int)
    sgn = np.zeros((mesh.nt, 3))
    eindex = {}
    for ei, (a, b) in enumerate(mesh.edges.tolist()):
        eindex[(a, b)] = ei
    for t, v in enumerate(mesh.tris):
        for il in range(3):
            a, b = int(v[(il + 1) % 3]), int(v[(il + 2) % 3])
            key = (a, b) if a < b else (b, a)
            ei = eindex[key]; ge[t, il] = ei
            sgn[t, il] = 1.0 if mesh.edge_tris[ei][0][0] == t else -1.0
    return ge, sgn


def assemble_rt0_2d(mesh, sigma):
    sigma = np.broadcast_to(np.asarray(sigma, float), (mesh.nt,))
    ge, sgn = _edges_and_signs(mesh)
    P = mesh.points[mesh.tris]                       # (Nt,3,2)
    area = np.abs(mesh.area)
    c = P.mean(axis=1)
    s2 = (((P - c[:, None, :]) ** 2).sum(-1)).sum(-1) / 12.0
    rows_a, cols_a, vals_a, rows_b, cols_b, vals_b = [], [], [], [], [], []
    for t in range(mesh.nt):
        d = c[t][None, :] - P[t]                     # (3,2) c - p_i
        dot = d @ d.T
        Ae = (s2[t] + dot) / (4.0 * sigma[t] * area[t])
        Ae = (sgn[t][:, None] * sgn[t][None, :]) * Ae
        for a in range(3):
            for b in range(3):
                rows_a.append(ge[t, a]); cols_a.append(ge[t, b]); vals_a.append(Ae[a, b])
            rows_b.append(t); cols_b.append(ge[t, a]); vals_b.append(sgn[t, a])
    A = sp.coo_matrix((vals_a, (rows_a, cols_a)), shape=(mesh.ne, mesh.ne)).tocsr()
    B = sp.coo_matrix((vals_b, (rows_b, cols_b)), shape=(mesh.nt, mesh.ne)).tocsr()
    return A, B


def energy_lower_2d(mesh, sigma, dirichlet):
    """Complementary lower bound of the field energy (W = 1/2 v^T K v)."""
    L = float(np.ptp(mesh.points, axis=0).max())
    work = mesh
    A, B = assemble_rt0_2d(work, sigma)
    ne, nt = work.ne, work.nt
    g = np.zeros(ne); pinned = np.zeros(ne, bool)
    for ei in range(ne):
        if not work.boundary_edge[ei]:
            continue
        a, b = work.edges[ei]
        if a in dirichlet and b in dirichlet:
            g[ei] = -0.5 * (dirichlet[int(a)] + dirichlet[int(b)])
        else:
            pinned[ei] = True                        # insulated edge: J.n = 0
    keep = ~pinned
    A = A[keep][:, keep]; B = B[:, keep]; g = g[keep]; ne = int(keep.sum())
    K = sp.bmat([[A, B.T], [B, None]], format="csc")
    sol = spla.spsolve(K, np.concatenate([g, np.zeros(nt)]))
    D = sol[:ne]
    return 0.5 * float(np.sum(D * (A @ D)))
