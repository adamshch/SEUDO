"""Benchmark the native_nthreads parameter (see estimate_time_courses_with_seudo)
on real, human-labeled data from code/demoData1.mat.

native_nthreads controls parallelism *within* a single SEUDO solve -- the
realSEUDO paper's two-pass gradient computation is explicitly designed to
be parallelizable across pixels/dimensions via the vendored POSIX-thread
("TPOPP") machinery. This is a different axis from n_jobs (which
parallelizes *across* independent (frame, cell) solves) -- both already
exist as parameters, but native_nthreads has never been benchmarked, and
its interaction with n_jobs (oversubscription: n_jobs threads each trying
to spin up native_nthreads more threads) is untested.

Sweeps native_nthreads at a couple of n_jobs settings and measures, at each
combination:
  - wall-clock time per crop (a whole transient-window run, matching
    benchmark_blob_spacing.py/benchmark_native_l_mode.py's methodology)
  - FALSE-transient suppression and TRUE-transient correlation, included as
    a correctness sanity check -- native_nthreads should be a pure speed
    lever with numerically identical output; any real difference here would
    indicate a race condition, not an expected tradeoff.
"""

import sys
import time

import h5py
import numpy as np

sys.path.insert(0, '.')
from seudo.estimate import estimate_time_courses_with_seudo  # noqa: E402
from seudo.matlab_io import Hdf5Movie, load_profiles, load_time_courses  # noqa: E402

DEMO_PATH = '../code/demoData1.mat'
N_JOBS_VALUES = [1, 8]
NTHREADS_VALUES = [1, 2, 4, 8]
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


def run_one(movie, profiles, cell_idx, matlab_start, matlab_end, n_jobs, native_nthreads):
    start0, end0 = int(matlab_start) - 1, int(matlab_end) - 1
    frame_lo = max(0, start0 - BASELINE)
    frame_hi = min(movie.shape[2], end0 + BASELINE + 1)

    movie_crop = movie.get_frames(range(frame_lo, frame_hi))

    t0 = time.time()
    result = estimate_time_courses_with_seudo(
        movie_crop, profiles, which_cells=[cell_idx],
        n_jobs=n_jobs, native_nthreads=native_nthreads, **SEUDO_PARAMS,
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

        combos = [(nj, nt) for nj in N_JOBS_VALUES for nt in NTHREADS_VALUES]
        per_combo = {c: dict(elapsed=[], false_suppression=[], true_correlation=[]) for c in combos}

        cells_tested = 0
        for cell_idx in range(n_cells):
            if cells_tested >= N_CELLS:
                break
            false_times = load_transient_times(f, cell_idx, -1.0)
            true_times = load_transient_times(f, cell_idx, 1.0)
            if false_times.shape[1] == 0 or true_times.shape[1] == 0:
                continue

            cells_tested += 1
            print(f'cell {cell_idx}: {false_times.shape[1]} false, {true_times.shape[1]} true transient(s)',
                  flush=True)

            ms, me = false_times[:, 0]
            for nj, nt in combos:
                tc_seudo, _lo, _hi, sl, elapsed = run_one(movie, profiles, cell_idx, ms, me, nj, nt)
                per_combo[(nj, nt)]['elapsed'].append(elapsed)
                per_combo[(nj, nt)]['false_suppression'].append(float(np.mean(np.abs(tc_seudo[sl]))))

            ms, me = true_times[:, 0]
            for nj, nt in combos:
                tc_seudo, lo, hi, _sl, elapsed = run_one(movie, profiles, cell_idx, ms, me, nj, nt)
                per_combo[(nj, nt)]['elapsed'].append(elapsed)
                t_crop = T[lo:hi, cell_idx]
                corr = np.corrcoef(t_crop, tc_seudo)[0, 1]
                per_combo[(nj, nt)]['true_correlation'].append(float(corr))

        print(f'\n{cells_tested} cells benchmarked, {len(combos)} (n_jobs, native_nthreads) combos each\n')
        print(f'{"n_jobs":>7} {"nthreads":>9} {"ms/crop":>10} {"speedup":>9} '
              f'{"false |activity|":>17} {"true corr(T)":>13}')
        base_ms = 1000 * np.mean(per_combo[(1, 1)]['elapsed'])
        for nj, nt in combos:
            d = per_combo[(nj, nt)]
            mean_ms = 1000 * np.mean(d['elapsed'])
            speedup = base_ms / mean_ms
            mean_false = np.mean(d['false_suppression'])
            mean_true_corr = np.mean(d['true_correlation'])
            print(f'{nj:>7} {nt:>9} {mean_ms:>10.1f} {speedup:>8.2f}x '
                  f'{mean_false:>17.4f} {mean_true_corr:>13.4f}')
        print('\n(speedup relative to n_jobs=1, native_nthreads=1)')


if __name__ == '__main__':
    main()
