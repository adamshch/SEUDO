"""Benchmark the native_l_mode parameter (see estimate_time_courses_with_seudo)
on real, human-labeled data from code/demoData1.mat.

Cross-checking the realSEUDO paper (Dmitrieva, Babkin & Charles, NeurIPS
2024, Sec 3.1) against the actual vendored C++ source found two FISTA
optimizations it describes -- per-dimension ("multi-L") Lipschitz constant
estimation ("~30% fewer steps") and momentum "fast brake" (dimensional
momentum reset on gradient sign change + the eta=1 simplification) -- both
already implemented (fista_v1_minimizer.hpp's MultiGradientMinimizer;
fista.cpp's Run::oneStep()) and already exposed via native_l_mode:
    0: dynamic L                   1: static multi-L
    2: dynamic L + fast brake      3: static multi-L + fast brake
Our own default everywhere is native_l_mode=2 -- fast brake, but *not*
multi-L -- so mode 3 may be leaving that ~30% on the table. This script
checks whether it actually helps on real data, the same way
benchmark_blob_spacing.py checked blob_spacing.

For several cells' real labeled TRUE and FALSE transients, sweeps
native_l_mode and measures, at each value:
  - wall-clock time per solve
  - FALSE-transient suppression: mean |activity| during a labeled
    contamination event (want this low)
  - TRUE-transient preservation: correlation of the SEUDO trace against the
    provided time course T over the whole crop (want this high)
"""

import sys
import time

import h5py
import numpy as np

sys.path.insert(0, '.')
from seudo.estimate import estimate_time_courses_with_seudo  # noqa: E402
from seudo.matlab_io import Hdf5Movie, load_profiles, load_time_courses  # noqa: E402

DEMO_PATH = '../code/demoData1.mat'
L_MODES = [0, 1, 2, 3]
L_MODE_LABELS = {0: 'dynamic L', 1: 'static multi-L', 2: 'dynamic L + fast brake', 3: 'multi-L + fast brake'}
N_CELLS = 8      # how many (false+true)-labeled cells to sample
BASELINE = 30    # frames of context loaded before/after each transient

SEUDO_PARAMS = dict(sigma2=0.0020, lambda_blob=10.0, blob_radius=3.0, pad_space=5,
                     solver_max_iter=300, use_native=True, verbose=False)


def load_transient_times(f, cell_idx, classification_value):
    ti = f['tDefaultDemo']['transientInfo']
    cls = np.array(f[ti['classification'][0, cell_idx]][()]).ravel()
    times = np.array(f[ti['times'][0, cell_idx]][()])  # (2, nTrans), MATLAB 1-indexed [start; end]
    idx = np.where(cls == classification_value)[0]
    return times[:, idx]


def run_one(movie, profiles, cell_idx, matlab_start, matlab_end, native_l_mode):
    start0, end0 = int(matlab_start) - 1, int(matlab_end) - 1
    frame_lo = max(0, start0 - BASELINE)
    frame_hi = min(movie.shape[2], end0 + BASELINE + 1)

    movie_crop = movie.get_frames(range(frame_lo, frame_hi))

    t0 = time.time()
    result = estimate_time_courses_with_seudo(
        movie_crop, profiles, which_cells=[cell_idx],
        native_l_mode=native_l_mode, **SEUDO_PARAMS,
    )
    elapsed = time.time() - t0

    tc_seudo = result['tc'][:, 0]
    transient_slice = slice(start0 - frame_lo, end0 - frame_lo + 1)
    return tc_seudo, frame_lo, frame_hi, transient_slice, elapsed


def main():
    with h5py.File(DEMO_PATH, 'r') as f:
        profiles = load_profiles(f)
        T = load_time_courses(f, 'T')
        movie = Hdf5Movie(f['M'])
        n_cells = profiles.shape[2]

        per_mode = {lm: dict(elapsed=[], false_suppression=[], true_correlation=[]) for lm in L_MODES}

        cells_tested = 0
        for cell_idx in range(n_cells):
            if cells_tested >= N_CELLS:
                break
            false_times = load_transient_times(f, cell_idx, -1.0)
            true_times = load_transient_times(f, cell_idx, 1.0)
            if false_times.shape[1] == 0 or true_times.shape[1] == 0:
                continue

            cells_tested += 1
            print(f'cell {cell_idx}: {false_times.shape[1]} false, {true_times.shape[1]} true transient(s)')

            ms, me = false_times[:, 0]
            for lm in L_MODES:
                tc_seudo, _lo, _hi, sl, elapsed = run_one(movie, profiles, cell_idx, ms, me, lm)
                per_mode[lm]['elapsed'].append(elapsed)
                per_mode[lm]['false_suppression'].append(float(np.mean(np.abs(tc_seudo[sl]))))

            ms, me = true_times[:, 0]
            for lm in L_MODES:
                tc_seudo, lo, hi, _sl, elapsed = run_one(movie, profiles, cell_idx, ms, me, lm)
                per_mode[lm]['elapsed'].append(elapsed)
                t_crop = T[lo:hi, cell_idx]
                corr = np.corrcoef(t_crop, tc_seudo)[0, 1]
                per_mode[lm]['true_correlation'].append(float(corr))

        print(f'\n{cells_tested} cells benchmarked, {len(L_MODES)} native_l_mode values each\n')
        print(f'{"l_mode":>7} {"label":>23} {"ms/solve":>10} {"speedup":>9} '
              f'{"false |activity|":>17} {"true corr(T)":>13}')
        base_time = per_mode[2]['elapsed']  # our current default, as the comparison baseline
        base_ms = 1000 * np.mean(base_time)
        for lm in L_MODES:
            d = per_mode[lm]
            mean_ms = 1000 * np.mean(d['elapsed'])
            speedup = base_ms / mean_ms
            mean_false = np.mean(d['false_suppression'])
            mean_true_corr = np.mean(d['true_correlation'])
            print(f'{lm:>7} {L_MODE_LABELS[lm]:>23} {mean_ms:>10.1f} {speedup:>8.2f}x '
                  f'{mean_false:>17.4f} {mean_true_corr:>13.4f}')
        print('\n(speedup relative to l_mode=2, our current default)')


if __name__ == '__main__':
    main()
