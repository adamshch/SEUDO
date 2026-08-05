"""Validate compute_transient_info against real MATLAB-computed reference
numbers from code/demoData1.mat.

demo.m calls:
    se = seudo(M,P,'timeCourses',T,'name','demo dataset')
    se.computeTransientInfo('default','transientFrames',TF)

i.e. transientFrames is supplied directly (bypassing identify_transients),
using default winRadius/tPre/tPost/useCOM. tDefaultDemo.transientInfo in
the file holds MATLAB's own computed times/window/corrWithProfile for the
*classified* run, which used the same TF and defaults -- so for cells where
classification didn't change the transient boundaries, this is a genuine
numeric-parity check, not just a behavioral one.

Restricted to a single cell (cell 0) and only the movie frames its own
transients touch, using the lazy Hdf5Movie loader -- avoids materializing
the full 8GB+ movie.
"""

import sys

import h5py
import numpy as np

sys.path.insert(0, '.')
from seudo.core import Seudo  # noqa: E402
from seudo.matlab_io import Hdf5Movie, load_profiles, load_time_courses, load_transient_frames  # noqa: E402

DEMO_PATH = '../code/demoData1.mat'
CELL_IDX = 0


def load_reference_transient_info(f, cell_idx):
    ti = f['tDefaultDemo']['transientInfo']
    times = np.array(f[ti['times'][0, cell_idx]][()]).T.astype(int) - 1  # -> 0-indexed [start, end]
    window = np.array(f[ti['window'][0, cell_idx]][()]).ravel().astype(int) - 1
    corr = np.array(f[ti['corrWithProfile'][0, cell_idx]][()]).ravel()
    return times, window, corr


def main():
    with h5py.File(DEMO_PATH, 'r') as f:
        profiles_full = load_profiles(f)
        T_full = load_time_courses(f, 'T')
        TF_full = load_transient_frames(f, 'TF')

        ref_times, ref_window, ref_corr = load_reference_transient_info(f, CELL_IDX)
        print(f'MATLAB reference: {ref_times.shape[0]} transients, window={ref_window}')

        # single-cell subset -- compute_transient_info treats cells independently,
        # so this is equivalent to running on the full 52-cell dataset for cell 0.
        profiles = profiles_full[:, :, [CELL_IDX]]
        T = T_full[:, [CELL_IDX]]
        TF = TF_full[:, [CELL_IDX]]

        movie = Hdf5Movie(f['M'])
        se = Seudo(movie, profiles, time_courses=T)

        info = se.compute_transient_info('default', transient_frames=TF)
        cell_info = info[0]

        print(f'Python result:    {cell_info["times"].shape[0]} transients, window={cell_info["window"]}')

        times_match = np.array_equal(cell_info['times'], ref_times)
        window_match = tuple(cell_info['window']) == tuple(ref_window)
        print(f'\ntimes match MATLAB exactly:  {times_match}')
        print(f'window matches MATLAB exactly: {window_match}')

        if not times_match:
            n = min(len(cell_info['times']), len(ref_times))
            diff = cell_info['times'][:n] - ref_times[:n]
            print('first mismatches (python - matlab):')
            print(diff[np.any(diff != 0, axis=1)][:10])

        corr = cell_info['corr_with_profile']
        n = min(len(corr), len(ref_corr))
        corr_diff = np.abs(corr[:n] - ref_corr[:n])
        print(f'\ncorrWithProfile: max abs diff over first {n} transients = {np.nanmax(corr_diff):.6f}')
        print(f'  python corr[:5] = {np.round(corr[:5], 4)}')
        print(f'  matlab corr[:5] = {np.round(ref_corr[:5], 4)}')


if __name__ == '__main__':
    main()
