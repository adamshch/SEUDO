"""Gaussian blob basis function, matching the blob built in
@seudo/estimateTimeCoursesWithSEUDO.m (lines ~262-270)."""

import numpy as np


def matlab_fspecial_gaussian(hsize, sigma):
    """Reimplementation of MATLAB's fspecial('gaussian', hsize, sigma)
    for a square, odd-sized kernel."""
    siz = (hsize - 1) / 2
    y, x = np.mgrid[-siz:siz + 1, -siz:siz + 1]
    h = np.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2))
    h[h < np.finfo(float).eps * h.max()] = 0.0
    total = h.sum()
    if total != 0:
        h = h / total
    return h


def make_seudo_blob(blob_radius, clip_height=0.01):
    """Build the single normalized Gaussian blob used as the basis function
    for unmodeled ("blob") activity in SEUDO."""
    crop_rad = int(np.ceil(blob_radius * 2.5 + np.finfo(float).eps))
    hsize = crop_rad * 2 + 1
    blob = matlab_fspecial_gaussian(hsize, blob_radius)
    blob = blob * (blob > clip_height * blob.max())
    blob = blob / np.sqrt(np.sum(blob ** 2))
    return blob
