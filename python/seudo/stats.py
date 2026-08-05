"""Small statistics helpers, ports of seudo.robustSTD (@seudo/seudo.m,
lines ~1351-1394) and extras/correlationVectorMatrix.m."""

import numpy as np


def robust_std(data):
    """Robust standard deviation: median absolute deviation, rescaled so it
    matches the ordinary std for normally-distributed data.

    data: 1D array (NaNs ignored).
    """
    data = np.asarray(data, dtype=float)
    if data.size == 0:
        return np.array([])
    centered = data - np.nanmedian(data)
    return np.nanmedian(np.abs(centered)) / 0.6741891400433162


def correlation_vector_matrix(v, m):
    """Pearson correlation of vector v (length T) against each column of
    matrix m (T x C), ignoring NaNs pairwise (per column).

    Returns a length-C array of correlations.
    """
    v = np.asarray(v, dtype=float).reshape(-1)
    m = np.asarray(m, dtype=float)
    assert v.shape[0] == m.shape[0]

    v = v - np.nanmean(v)
    m = m - np.nanmean(m, axis=0, keepdims=True)

    numers = v[:, np.newaxis] * m
    denom_v = np.nansum(v ** 2)
    denom_m = np.nansum(m ** 2, axis=0)

    with np.errstate(invalid='ignore', divide='ignore'):
        corrs = np.nansum(numers, axis=0) / np.sqrt(denom_v * denom_m)

    return corrs
