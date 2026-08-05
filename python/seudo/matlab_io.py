"""Loading helpers for MATLAB v7.3 (HDF5) .mat files, as saved by this
toolbox's demo data.

MATLAB stores arrays column-major; h5py reads them back with every axis
reversed relative to how MATLAB indexes them. So a MATLAB array of size
(movY, movX, movF) appears in h5py as shape (movF, movX, movY). This was
confirmed empirically (see python/scripts/run_on_demo_data.py) by comparing
a profile's bounding box, after transposing back, against the padding
radius baked into a precomputed transientInfo.window field.
"""

import numpy as np


class Hdf5Movie:
    """Lazy frame access into an h5py Dataset holding a MATLAB (movY, movX,
    movF) movie array (stored/read as (movF, movX, movY) -- see module
    docstring). Duck-type compatible with what Seudo expects: `.shape` as
    (movY, movX, movF) and `.get_frame(i)` / `.get_frames(idx)`.
    """

    def __init__(self, dataset):
        self.dataset = dataset
        f, x, y = dataset.shape
        self.shape = (y, x, f)

    def get_frame(self, frame_index):
        return self.dataset[frame_index, :, :].T.astype(float)

    def get_frames(self, frame_indices):
        frame_indices = list(frame_indices)
        raw = self.dataset[frame_indices, :, :]
        return np.transpose(raw, (2, 1, 0))


def load_profiles(h5file):
    """h5file: an open h5py.File. Returns profiles as (movY, movX, nCells)."""
    raw = h5file['P']  # (nCells, movX, movY)
    return np.transpose(raw[:], (2, 1, 0))


def load_time_courses(h5file, key='T'):
    """Returns (movF, nCells)."""
    raw = h5file[key]  # (nCells, movF)
    return np.transpose(raw[:], (1, 0))


def load_transient_frames(h5file, key='TF'):
    raw = h5file[key]  # (nCells, movF), uint8
    return np.transpose(raw[:], (1, 0)).astype(bool)
