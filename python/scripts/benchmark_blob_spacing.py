"""Benchmark the blob_spacing parameter (see estimate_time_courses_with_seudo)
on real, human-labeled data from code/demoData1.mat: does a coarser blob
grid actually give the speed/quality tradeoff the realSEUDO paper
(Dmitrieva, Babkin & Charles, NeurIPS 2024, Sec 3.1) reports -- "more than
quadratic" FISTA speedup, over 100x for a 30px-diameter kernel, "without
substantial degradation" of false-transient rejection or true-signal
recognition?

For several cells' real labeled TRUE and FALSE transients (see
run_on_demo_data.py for the same crop-around-a-transient pattern this
reuses), sweeps blob_spacing and measures, at each value:
  - wall-clock time per solve
  - FALSE-transient suppression: mean |activity| during a labeled
    contamination event (want this low -- SEUDO should still reject it)
  - TRUE-transient preservation: correlation of the SEUDO trace against the
    provided time course T over the whole crop (want this high -- SEUDO
    should still track genuine signal shape)
"""

import sys
import time

import h5py
import numpy as np

sys.path.insert(0, '.')
from seudo.estimate import estimate_time_courses_with_seudo  # noqa: E402
from seudo.matlab_io import Hdf5Movie, load_profiles, load_time_courses  # noqa: E402

DEMO_PATH = '../code/demoData1.mat'
BLOB_SPACINGS = [1, 2, 3, 4, 6, 8]
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


def run_one(movie, profiles, cell_idx, matlab_start, matlab_end, blob_spacing):
    start0, end0 = int(matlab_start) - 1, int(matlab_end) - 1
    frame_lo = max(0, start0 - BASELINE)
    frame_hi = min(movie.shape[2], end0 + BASELINE + 1)

    movie_crop = movie.get_frames(range(frame_lo, frame_hi))

    t0 = time.time()
    result = estimate_time_courses_with_seudo(
        movie_crop, profiles, which_cells=[cell_idx],
        blob_spacing=blob_spacing, **SEUDO_PARAMS,
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

        per_spacing = {bs: dict(elapsed=[], false_suppression=[], true_correlation=[]) for bs in BLOB_SPACINGS}

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
            for bs in BLOB_SPACINGS:
                tc_seudo, _lo, _hi, sl, elapsed = run_one(movie, profiles, cell_idx, ms, me, bs)
                per_spacing[bs]['elapsed'].append(elapsed)
                per_spacing[bs]['false_suppression'].append(float(np.mean(np.abs(tc_seudo[sl]))))

            ms, me = true_times[:, 0]
            for bs in BLOB_SPACINGS:
                tc_seudo, lo, hi, _sl, elapsed = run_one(movie, profiles, cell_idx, ms, me, bs)
                per_spacing[bs]['elapsed'].append(elapsed)
                t_crop = T[lo:hi, cell_idx]
                corr = np.corrcoef(t_crop, tc_seudo)[0, 1]
                per_spacing[bs]['true_correlation'].append(float(corr))

        print(f'\n{cells_tested} cells benchmarked, {len(BLOB_SPACINGS)} blob_spacing values each\n')
        print(f'{"blob_spacing":>12} {"ms/solve":>10} {"speedup":>9} '
              f'{"false |activity|":>17} {"true corr(T)":>13}')
        base_time = None
        for bs in BLOB_SPACINGS:
            d = per_spacing[bs]
            mean_ms = 1000 * np.mean(d['elapsed'])
            if base_time is None:
                base_time = mean_ms
            speedup = base_time / mean_ms
            mean_false = np.mean(d['false_suppression'])
            mean_true_corr = np.mean(d['true_correlation'])
            print(f'{bs:>12} {mean_ms:>10.1f} {speedup:>8.2f}x '
                  f'{mean_false:>17.4f} {mean_true_corr:>13.4f}')


if __name__ == '__main__':
    main()
