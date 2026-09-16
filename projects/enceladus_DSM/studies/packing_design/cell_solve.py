#!/usr/bin/env python3
"""Periodic homogenization cell problem on a raster, in numpy.

A small finite-volume twin of what src/keff_cell.c does with PetIGA. It exists
so the BIAS questions -- how much does a spurious solid bridge move k_eff, how
much does pore fragmentation cost vapour transport -- can be answered in
seconds on a laptop, at t = 0, without queueing a solver run. It is not a
replacement for the real solver and is not used by it.

For direction m it solves, on a periodic N x N grid of cell size h,

    div( K (grad t + e_m) ) = 0

by the standard finite-volume discretisation with HARMONIC face
conductivities -- the right average for a flux crossing two cells in series,
and the one that stays sane when the contrast is ~100 as it is for ice/air.
The resulting matrix is the weighted graph Laplacian, singular with the
constant null space, and the load is orthogonal to it (the face terms
telescope on a torus), so CG on the consistent singular system is well posed
and the arbitrary additive constant is simply never determined.

Then

    K_eff[m][n] = < K ( dt_m/dx_n + delta_mn ) >

averaged over faces. Used for two different coefficients:

    K(phi)          -> effective thermal conductivity
    phi_air         -> effective vapour diffusivity, i.e. how well the pore
                       network actually transports, which is the continuous
                       replacement for a percolation test whose answer is
                       always 'disconnected'.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _harmonic(a, b):
    """Harmonic mean, the series average, guarding a+b = 0."""
    s = a + b
    return np.where(s > 0, 2.0 * a * b / np.where(s > 0, s, 1.0), 0.0)


def _faces(K):
    """(Kx, Ky) harmonic face conductivities; Kx[i,j] joins (i,j)-(i,j+1)."""
    return (_harmonic(K, np.roll(K, -1, axis=1)),
            _harmonic(K, np.roll(K, -1, axis=0)))


def _operator(Kx, Ky):
    """Weighted periodic graph Laplacian as a sparse matrix."""
    ny, nx = Kx.shape
    n = nx * ny
    idx = np.arange(n).reshape(ny, nx)
    rows, cols, vals = [], [], []
    for K, shift, axis in ((Kx, -1, 1), (Ky, -1, 0)):
        nb = np.roll(idx, shift, axis=axis)
        w = K.ravel()
        c, d = idx.ravel(), nb.ravel()
        # off-diagonals both ways, and the matching diagonal debits
        rows += [c, d, c, d]
        cols += [d, c, c, d]
        vals += [w, w, -w, -w]
    A = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n, n)).tocsr()
    return A


def k_eff(K, h):
    """Effective conductivity tensor of the periodic coefficient field K.

    K is (ny, nx) and strictly positive. Returns a 2x2 array.

    DIRECT SOLVE, not CG. Jacobi-preconditioned CG handles the thermal problem
    (contrast ~114) but fails outright on the vapour one, where the coefficient
    is phi_air and so is ~0 over the two-thirds of the domain that is ice --
    it did not converge in 4000 iterations at any floor value. Factoring once
    and reusing for both directions costs ~2 s at N = 384 and is exact, which
    matters because the whole point here is to resolve small biases.

    The operator is singular with the constant null space, so one unknown is
    pinned. The load is orthogonal to the null space (the face terms telescope
    on a torus), so pinning selects a particular solution and the gradients --
    the only thing used below -- are unaffected.
    """
    ny, nx = K.shape
    Kx, Ky = _faces(K)
    A = _operator(Kx, Ky).tolil()
    A[0, :] = 0.0
    A[0, 0] = 1.0
    lu = spla.splu(A.tocsc())

    out = np.zeros((2, 2))
    for m in (0, 1):
        Kf = Kx if m == 0 else Ky
        ax = 1 if m == 0 else 0
        b = -h * (Kf - np.roll(Kf, 1, axis=ax)).ravel()
        b -= b.mean()
        b[0] = 0.0
        t = lu.solve(b).reshape(ny, nx)
        gx = (np.roll(t, -1, axis=1) - t) / h          # at x-faces
        gy = (np.roll(t, -1, axis=0) - t) / h          # at y-faces
        out[m, 0] = float(np.mean(Kx * (gx + (1.0 if m == 0 else 0.0))))
        out[m, 1] = float(np.mean(Ky * (gy + (1.0 if m == 1 else 0.0))))
    return out
