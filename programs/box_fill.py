"""Laplace fill with Dirichlet data on BOTH the measured set and the box edge.

    grad^2 u = 0   in the unknown region
    u = f          on the measured pixels
    u = mu         on the edge of the box

mu is the mean of the measured pixels. No free parameter: unlike the screened
form there is no relaxation length to choose, because the box edge anchors the
far field directly. Forcing every edge of the tile to the same constant also
makes the tile continuous when the FFT treats it as periodic.

Differs from periodic_harmonic_fill.py in one place: a neighbour that falls
outside the array is not wrapped, it is a known boundary value mu.
"""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import cg


def box_dirichlet_fill(image, measured, *, mu=None, use_amg=False,
                       rtol=1e-8, maxiter=None):
    image = np.asarray(image)
    measured = np.asarray(measured)
    if image.ndim != 2 or measured.shape != image.shape or measured.dtype != bool:
        raise ValueError("Expected a 2D image and a same-shape Boolean mask.")
    if not measured.any():
        raise ValueError("At least one measured pixel is required.")

    result = image.astype(np.float64, copy=True)
    missing = ~measured
    if not missing.any():
        return result
    if mu is None:
        mu = float(result[measured].mean())

    height, width = image.shape
    y, x = np.nonzero(missing)
    unknown = np.arange(y.size)
    index = np.full(image.shape, -1, dtype=np.intp)
    index[y, x] = unknown
    rhs = np.zeros(y.size)
    rows = [unknown]; cols = [unknown]; weights = [np.full(y.size, 4.0)]

    for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        ny, nx = y + dy, x + dx
        outside = (ny < 0) | (ny >= height) | (nx < 0) | (nx >= width)
        rhs[outside] += mu                      # the box edge is held at mu
        iny, inx = np.clip(ny, 0, height-1), np.clip(nx, 0, width-1)
        neighbor = np.where(outside, -1, index[iny, inx])
        known = (~outside) & (neighbor < 0)     # a measured neighbour
        rhs[known] += result[iny[known], inx[known]]
        unk = (~outside) & (neighbor >= 0)
        rows.append(unknown[unk]); cols.append(neighbor[unk])
        weights.append(np.full(int(unk.sum()), -1.0))

    matrix = coo_matrix(
        (np.concatenate(weights), (np.concatenate(rows), np.concatenate(cols))),
        shape=(y.size, y.size)).tocsr()

    M = None
    if use_amg:
        from pyamg import smoothed_aggregation_solver
        M = smoothed_aggregation_solver(matrix, symmetry="symmetric").aspreconditioner(cycle="V")
    filled, info = cg(matrix, rhs, M=M, rtol=rtol, atol=0.0, maxiter=maxiter)
    if info != 0:
        raise RuntimeError(f"Box fill failed to converge (CG info={info}).")
    result[missing] = filled
    return result
