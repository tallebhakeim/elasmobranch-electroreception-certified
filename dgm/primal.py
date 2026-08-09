"""Primal DGM formulation (scalar potential on primal-mesh nodes).

This is the formulation  G^T M_eps^P G v = q  of Ren et al. (2015), Eq. (4).
The paper's Appendix A proves the elementary matrix is identical to the nodal
(Galerkin) FEM cotangent stiffness, so this single assembly serves both
"primal DGM" and "FEM primal".  It yields the UPPER bound of the field energy.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

EPS0 = 8.8541878128e-12  # vacuum permittivity [F/m]


def assemble_primal(mesh, eps_tri):
    """Assemble the primal stiffness K (Np x Np), sparse symmetric.

    Parameters
    ----------
    mesh : TriMesh
    eps_tri : (Nt,) array of absolute permittivities (= EPS0 * eps_r) per triangle.

    Element entry on the edge opposite local vertex i:  w_i = eps * l_D_i / l_P_i
    (= eps/2 * cot(angle_i)), matching Ren et al. Eq. (6),(10).
    """
    eps_tri = np.broadcast_to(np.asarray(eps_tri), (mesh.nt,))   # keep dtype (complex OK)
    dt = eps_tri.dtype
    rows, cols, vals = [], [], []
    for t in range(mesh.nt):
        v = mesh.tris[t]
        w = np.empty(3, dtype=dt)
        for i_local in range(3):
            lp, ld = mesh._local_edge_geometry(t, i_local)
            w[i_local] = eps_tri[t] * ld / lp
        # off-diagonal couplings: edge opposite vertex i couples the other two
        ke = np.zeros((3, 3), dtype=dt)
        for i_local in range(3):
            j, k = (i_local + 1) % 3, (i_local + 2) % 3
            ke[j, k] -= w[i_local]
            ke[k, j] -= w[i_local]
            ke[j, j] += w[i_local]
            ke[k, k] += w[i_local]
        for a in range(3):
            for b in range(3):
                rows.append(v[a]); cols.append(v[b]); vals.append(ke[a, b])
    K = sp.coo_matrix((vals, (rows, cols)), shape=(mesh.np, mesh.np)).tocsr()
    return K


def solve_dirichlet(K, dirichlet, load=None):
    """Solve K v = load with prescribed nodal values.

    Parameters
    ----------
    K : (N, N) sparse
    dirichlet : dict {node_index: value}
    load : (N,) array or None. Right-hand side (e.g. a heat source). Defaults to 0.

    Returns
    -------
    v : (N,) full solution vector.
    """
    n = K.shape[0]
    fixed = np.array(sorted(dirichlet), dtype=int)
    vals = np.array([dirichlet[i] for i in fixed], float)
    free = np.setdiff1d(np.arange(n), fixed)

    b = np.zeros(n) if load is None else np.asarray(load, float)
    v = np.zeros(n)
    v[fixed] = vals
    Kff = K[free][:, free]
    Kfx = K[free][:, fixed]
    rhs = b[free] - Kfx @ vals
    v[free] = spla.spsolve(Kff.tocsc(), rhs)
    return v


def energy(K, v):
    """Electrostatic field energy  W = 1/2 v^T K v  [J/m in 2D].

    Uses sum(v * (K@v)) rather than the `@` quadratic form: NumPy 2.x matmul
    spuriously raises floating-point flags from BLAS internals on this path,
    even though the result is exact. The elementwise form avoids that.
    """
    return 0.5 * float(np.sum(v * (K @ v)))


def capacitance_two_conductor(K, v):
    """For a two-electrode system solved at 1 V / 0 V, C = 2 W / V^2 = 2 W."""
    return 2.0 * energy(K, v)
