"""Fill missing pixels of a periodic 2D FFT tile by harmonic inpainting."""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import cg


def periodic_harmonic_fill(image, measured, *, x0=None, use_amg=False,
                           rtol=1e-8, maxiter=None):
    """Return a float64 tile with measured pixels copied and gaps filled.

    measured is a Boolean mask, True wherever a pixel is valid, including
    valid zeros. Values outside the mask are ignored. An optional full-tile
    x0 supplies the starting guess in the gaps. Set use_amg=True to use
    PyAMG as a preconditioner for wide gaps. Requires SciPy >= 1.12.
    """
    image = np.asarray(image)
    measured = np.asarray(measured)
    if image.ndim != 2 or measured.shape != image.shape or measured.dtype != bool:
        raise ValueError("Expected a 2D image and a same-shape Boolean mask.")
    if image.dtype.kind not in "buif":
        raise ValueError("Expected a real numeric image.")
    if not measured.any():
        raise ValueError("At least one measured pixel is required.")
    if not np.isfinite(image[measured]).all():
        raise ValueError("Measured pixels must be finite.")

    result = image.astype(np.float64, copy=True)
    missing = ~measured
    if not missing.any():
        return result

    height, width = image.shape
    y, x = np.nonzero(missing)
    unknown = np.arange(y.size)
    index = np.full(image.shape, -1, dtype=np.intp)
    index[y, x] = unknown
    rhs = np.zeros(y.size)
    rows = [unknown]
    cols = [unknown]
    weights = [np.full(y.size, 4.0)]

    for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        ny, nx = (y + dy) % height, (x + dx) % width
        neighbor = index[ny, nx]
        fixed = neighbor < 0
        rhs[fixed] += result[ny[fixed], nx[fixed]]
        rows.append(unknown[~fixed])
        cols.append(neighbor[~fixed])
        weights.append(np.full(np.count_nonzero(~fixed), -1.0))

    # COO sums duplicate entries, including on one- or two-pixel axes.
    matrix = coo_matrix(
        (np.concatenate(weights), (np.concatenate(rows), np.concatenate(cols))),
        shape=(y.size, y.size),
    ).tocsr()

    guess = None
    if x0 is not None:
        x0 = np.asarray(x0, dtype=np.float64)
        if x0.shape != image.shape or not np.isfinite(x0[missing]).all():
            raise ValueError("x0 must match the tile and be finite in the gaps.")
        guess = x0[missing]

    preconditioner = None
    if use_amg:
        from pyamg import smoothed_aggregation_solver

        preconditioner = smoothed_aggregation_solver(
            matrix, symmetry="symmetric"
        ).aspreconditioner(cycle="V")

    filled, info = cg(matrix, rhs, x0=guess, M=preconditioner,
                      rtol=rtol, atol=0.0, maxiter=maxiter)
    if info != 0:
        raise RuntimeError(f"Harmonic fill failed to converge (CG info={info}).")

    result[missing] = filled
    return result
