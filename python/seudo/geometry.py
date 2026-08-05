"""Port of the seudo.computeRoiCOMs static method (@seudo/seudo.m, lines
~1176-1217).

Coordinates here are 0-indexed (Pythonic), unlike the 1-indexed MATLAB
original -- the formulas are equivalent, just shifted.
"""

import numpy as np


def compute_roi_coms(rois):
    """Compute center-of-mass and bounding box for one or more ROIs.

    rois: (Y, X) or (Y, X, R) array

    Returns:
        coms: (R, 2) array, each row [comX, comY], 0-indexed
        outer_bounds: (R, 4) array, each row [y0, y1, x0, x1], 0-indexed
                       inclusive bounds of the nonzero region. If an ROI is
                       all zero, bounds default to [0, 0, 0, 0].
    """
    rois = np.asarray(rois)
    if rois.ndim == 2:
        rois = rois[:, :, np.newaxis]

    n_y, n_x, n_rois = rois.shape
    yy, xx = np.mgrid[0:n_y, 0:n_x]

    coms = np.full((n_rois, 2), np.nan)
    outer_bounds = np.zeros((n_rois, 4), dtype=int)

    for rr in range(n_rois):
        roi = np.abs(rois[:, :, rr])

        total = roi.sum()
        if total > 0:
            com_x = np.sum(xx * roi) / total
            com_y = np.sum(yy * roi) / total
        else:
            com_x = com_y = np.nan
        coms[rr] = [com_x, com_y]

        row_has = np.any(roi != 0, axis=1)
        col_has = np.any(roi != 0, axis=0)
        if row_has.any():
            y0 = int(np.argmax(row_has))
            y1 = n_y - 1 - int(np.argmax(row_has[::-1]))
            x0 = int(np.argmax(col_has))
            x1 = n_x - 1 - int(np.argmax(col_has[::-1]))
        else:
            y0 = y1 = x0 = x1 = 0
        outer_bounds[rr] = [y0, y1, x0, x1]

    return coms, outer_bounds
