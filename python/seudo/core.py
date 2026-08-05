"""Minimal Python port of @seudo/seudo.m's Seudo class.

Implements enough to run estimate_time_courses_with_seudo and
compute_transient_info. Other MATLAB seudo functionality (tif movie
sources, LSQ time course computation, param search, classification/review
GUIs, ...) is not yet ported.
"""

import numpy as np

from .auto_classify import auto_classify_transients as _auto_classify_transients
from .estimate import estimate_time_courses_with_seudo as _estimate_time_courses_with_seudo
from .run_seudo_on_transients import run_seudo_restricted_to_transients as _run_seudo_restricted_to_transients
from .transients import compute_transient_info as _compute_transient_info


class Seudo:
    def __init__(self, movie, profiles, name='untitled movie', zero_level=0.0, time_courses=None):
        """movie: (Y, X, F) array, or an object exposing `.shape == (Y, X, F)`
        and `.get_frame(i)` for lazy access (see matlab_io.Hdf5Movie).
        profiles: (Y, X, nCells) array.
        time_courses: optional (F, nCells) array of externally-provided
        time courses, stored as tc_default['tc'] -- mirrors the MATLAB
        constructor's 'timeCourses' parameter.
        """
        if hasattr(movie, 'get_frame') and hasattr(movie, 'shape'):
            self.mov_y, self.mov_x, self.mov_f = movie.shape
            self._lazy_movie = True
        else:
            movie = np.asarray(movie)
            if movie.ndim != 3:
                raise ValueError('movie must be a (Y, X, F) array')
            self.mov_y, self.mov_x, self.mov_f = movie.shape
            self._lazy_movie = False

        profiles = np.asarray(profiles)
        if profiles.shape[:2] != (self.mov_y, self.mov_x):
            raise ValueError(
                f'size of profiles {profiles.shape[:2]} does not match size of movie '
                f'({self.mov_y}, {self.mov_x})'
            )

        self.movie = movie
        self.profiles = profiles
        self.n_cells = profiles.shape[2]
        self.zero_level = zero_level
        self.name = name

        self.tc_default = None
        self.tc_contam = None
        self.tc_seudo = []

        if time_courses is not None:
            time_courses = np.asarray(time_courses)
            if time_courses.shape != (self.mov_f, self.n_cells):
                raise ValueError(
                    f'time_courses shape {time_courses.shape} does not match '
                    f'(mov_f={self.mov_f}, n_cells={self.n_cells})'
                )
            self.tc_default = {'tc': time_courses}

    def get_frame(self, frame_index):
        if self._lazy_movie:
            return self.movie.get_frame(frame_index) - self.zero_level
        return self.movie[:, :, frame_index].astype(float) - self.zero_level

    def get_frames(self, frame_indices):
        if self._lazy_movie:
            return self.movie.get_frames(frame_indices) - self.zero_level
        return self.movie[:, :, frame_indices].astype(float) - self.zero_level

    def _resolve_tc_struct(self, which_struct):
        """Port of @seudo/pickTcStruct.m. which_struct is one of:
        'default', 'contam', 'seudo', or ('seudo', index) with a (possibly
        negative) index into tc_seudo. Returns the actual dict (mutable, so
        callers can fill in fields like 'transient_info' in place).
        """
        if isinstance(which_struct, str):
            key = which_struct.lower()
        elif isinstance(which_struct, (tuple, list)) and len(which_struct) == 2:
            key = which_struct[0].lower()
        else:
            raise ValueError(f'time course struct specification not recognized: {which_struct!r}')

        if key in ('default', 'tcdefault'):
            if self.tc_default is None:
                raise ValueError('tc_default has not been set (pass time_courses to the constructor)')
            return self.tc_default
        if key in ('contam', 'tccontam'):
            if self.tc_contam is None:
                raise ValueError('tc_contam has not been set')
            return self.tc_contam
        if key in ('seudo', 'tcseudo'):
            if not self.tc_seudo:
                raise ValueError('tc_seudo is empty')
            index = which_struct[1] if isinstance(which_struct, (tuple, list)) else -1
            return self.tc_seudo[index]

        raise ValueError(f'time course struct specification not recognized: {which_struct!r}')

    def estimate_time_courses_with_seudo(self, **kwargs):
        result = _estimate_time_courses_with_seudo(
            self.movie, self.profiles, zero_level=self.zero_level, **kwargs
        )
        self.tc_seudo.append(result)
        return result

    def compute_transient_info(self, which_struct='default', **kwargs):
        return _compute_transient_info(self, which_struct, **kwargs)

    def auto_classify_transients(self, which_struct='default', **kwargs):
        return _auto_classify_transients(self, which_struct, **kwargs)

    def run_seudo_restricted_to_transients(self, which_struct='default', **kwargs):
        return _run_seudo_restricted_to_transients(self, which_struct, **kwargs)
