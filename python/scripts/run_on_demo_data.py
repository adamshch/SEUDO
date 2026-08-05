"""Try the Python SEUDO port on real cells/transients from code/demoData1.mat.

demoData1.mat is a MATLAB v7.3 (HDF5) file. h5py reads MATLAB's column-major
arrays with axes fully reversed relative to MATLAB's own indexing, so an
array MATLAB sees as (movY, movX, nCells) shows up here with shape
(nCells, movX, movY). This was confirmed empirically: transposing a raw
profile slice and comparing its bounding box to the padding radius baked
into demoData1's precomputed transientInfo('default').window field lines up
(padSpace=10, matching computeTransientInfo's default) only when profiles
are read as P_h5[cc, :, :].T -> (movY, movX). See memory / conversation notes.
"""

import sys

import h5py
import numpy as np

sys.path.insert(0, '.')
from seudo.estimate import estimate_time_courses_with_seudo  # noqa: E402

DEMO_PATH = '../code/demoData1.mat'


def load_profiles(f):
    P_raw = f['P']  # (nCells, movX, movY)
    return np.transpose(P_raw[:], (2, 1, 0))  # -> (movY, movX, nCells)


def load_movie_crop(f, frame_lo, frame_hi):
    M_raw = f['M'][frame_lo:frame_hi, :, :]  # (F, movX, movY)
    return np.transpose(M_raw, (2, 1, 0))  # -> (movY, movX, F)


def load_transient(f, cell_idx, classification_value):
    ti = f['tDefaultDemo']['transientInfo']
    cls = np.array(f[ti['classification'][0, cell_idx]][()]).ravel()
    times = np.array(f[ti['times'][0, cell_idx]][()])  # (2, nTrans), MATLAB 1-indexed [start; end]
    idx = np.where(cls == classification_value)[0]
    return times[:, idx]


def run_case(f, profiles, cell_idx, matlab_start, matlab_end, label, baseline=30,
             seudo_params=None):
    seudo_params = seudo_params or {}
    start0 = int(matlab_start) - 1  # to 0-indexed
    end0 = int(matlab_end) - 1

    frame_lo = max(0, start0 - baseline)
    frame_hi = min(f['M'].shape[0], end0 + baseline + 1)

    print(f'\n=== {label}: cell {cell_idx}, MATLAB frames [{matlab_start:.0f},{matlab_end:.0f}], '
          f'loading crop frames [{frame_lo},{frame_hi}) ===')

    movie = load_movie_crop(f, frame_lo, frame_hi)
    T = np.array(f['T'][cell_idx, frame_lo:frame_hi])  # provided default time course

    result = estimate_time_courses_with_seudo(
        movie, profiles, which_cells=[cell_idx],
        verbose=False, solver_max_iter=300, **seudo_params,
    )
    tc_seudo = result['tc'][:, 0]
    tc_lsq = result['tc_lsq'][:, 0]

    t_start_local = start0 - frame_lo
    t_end_local = end0 - frame_lo
    transient_slice = slice(t_start_local, t_end_local + 1)
    baseline_slice = slice(0, max(1, t_start_local - 5))

    print(f'{"frame":>6} {"T(provided)":>12} {"our_lsq":>10} {"our_seudo":>10}')
    for i in range(0, len(tc_seudo), max(1, len(tc_seudo) // 25)):
        marker = '  <-- transient' if t_start_local <= i <= t_end_local else ''
        print(f'{frame_lo + i:6d} {T[i]:12.3f} {tc_lsq[i]:10.3f} {tc_seudo[i]:10.3f}{marker}')

    print(f'-- during transient [{t_start_local},{t_end_local}] --')
    print(f'  mean T           = {np.mean(T[transient_slice]):.3f}')
    print(f'  mean our_lsq     = {np.mean(tc_lsq[transient_slice]):.3f}')
    print(f'  mean our_seudo   = {np.mean(tc_seudo[transient_slice]):.3f}')
    print(f'-- baseline (pre-transient) --')
    print(f'  mean T           = {np.mean(T[baseline_slice]):.3f}')
    print(f'  mean our_lsq     = {np.mean(tc_lsq[baseline_slice]):.3f}')
    print(f'  mean our_seudo   = {np.mean(tc_seudo[baseline_slice]):.3f}')
    corr = np.corrcoef(T, tc_lsq)[0, 1]
    print(f'  corr(T, our_lsq) over whole crop = {corr:.4f}')


def main():
    with h5py.File(DEMO_PATH, 'r') as f:
        profiles = load_profiles(f)
        cell_idx = 0

        seudo_params = dict(ds_time=3, sigma2=0.0020, lambda_blob=10, blob_radius=3, pad_space=5)

        false_times = load_transient(f, cell_idx, -1.0)
        true_times = load_transient(f, cell_idx, 1.0)

        # a false (contaminant) transient -- SEUDO should suppress it
        ms, me = false_times[:, 0]
        run_case(f, profiles, cell_idx, ms, me, label='FALSE transient (should be suppressed)',
                  seudo_params=seudo_params)

        # a true transient -- SEUDO should preserve it
        ms, me = true_times[:, 0]
        run_case(f, profiles, cell_idx, ms, me, label='TRUE transient (should be preserved)',
                  seudo_params=seudo_params)


if __name__ == '__main__':
    main()
