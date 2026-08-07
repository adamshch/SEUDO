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


def make_smoothing_kernel(radius):
    """Sum-normalized (not L2-normalized) Gaussian smoothing kernel, for
    denoising a frame via convolution -- preserves the frame's overall
    pixel-intensity scale (mirrors matlab_fspecial_gaussian's own
    normalization), unlike make_seudo_blob's L2-normalized (sum of squares
    = 1) SEUDO basis-function convention, which would badly distort a
    frame's amplitude if used as a plain smoothing filter instead."""
    crop_rad = int(np.ceil(radius * 2.5 + np.finfo(float).eps))
    hsize = crop_rad * 2 + 1
    return matlab_fspecial_gaussian(hsize, radius)
