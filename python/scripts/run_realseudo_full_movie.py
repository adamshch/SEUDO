"""Run the streaming ROI-discovery pipeline (realSEUDOfit) on the FULL real
demo movie (code/demoData1.mat, all 41756 frames) starting from no prior
knowledge of cell locations, and save the discovered profiles/traces into
code/demoDataRealSEUDO.mat -- a fresh HDF5 file mirroring demoData1.mat's own
P/T dataset layout (same raw storage shapes/axis order, readable by the same
seudo.matlab_io.load_profiles/load_time_courses this codebase already uses),
so downstream tooling can treat our discovered result the same way it treats
the original CNMF ground truth.

Unlike compare_realseudo_to_cnmf.py (which processes a partial prefix of the
movie, N_FRAMES=8000, and only saves to a gitignored .npz for a later GUI
step), this script is meant to be the actual full-movie run -- now tractable
at a reasonable wall-clock time after this session's blob_spacing=3 default
and two-stage-detection performance work.
"""

import sys
import time

import h5py
import numpy as np

sys.path.insert(0, '.')
from seudo.matlab_io import Hdf5Movie, load_profiles, load_time_courses  # noqa: E402
from seudo.streaming import FitParams, StreamingState, realSEUDOfit  # noqa: E402

DEMO_PATH = '../code/demoData1.mat'
OUT_PATH = '../code/demoDataRealSEUDO.mat'

N_FRAMES = None    # None = the whole movie
PROGRESS_EVERY = 500

# same demo.m-derived defaults used elsewhere in this project for this
# dataset, plus the validated n_jobs/native_nthreads combo, the blob_spacing=3
# default, and ds_time=3 (a causal trailing 3-frame average, added to reduce
# per-pixel noise -- see FitParams.ds_time's docstring). blob_spacing and
# ds_time no longer need to be passed explicitly (both are FitParams
# defaults now), but spelled out here for a self-documenting record of
# what this run used.
SEUDO_PARAMS = dict(sigma2=0.0020, lambda_blob=10.0, blob_radius=3.0, pad_space=5,
                     n_jobs=8, native_nthreads=4, blob_spacing=3.0, ds_time=3)


def run_discovery(movie, n_frames):
    state = StreamingState(movie.shape[:2], fit=FitParams(**SEUDO_PARAMS))
    activity_by_frame = []

    t0 = time.time()
    for ff in range(n_frames):
        frame = movie.get_frame(ff)
        result = realSEUDOfit(frame, state)
        activity_by_frame.append(result.activity)

        if result.new_cells:
            print(f'  frame {ff}: discovered cell(s) {result.new_cells} '
                  f'(total so far: {state.profiles.shape[2]})', flush=True)

        if (ff + 1) % PROGRESS_EVERY == 0:
            elapsed = time.time() - t0
            rate = elapsed / (ff + 1)
            remaining = rate * (n_frames - (ff + 1))
            print(f'  ...{ff + 1}/{n_frames} frames, {rate * 1000:.1f} ms/frame, '
                  f'{state.profiles.shape[2]} cells so far, '
                  f'{remaining:.0f}s remaining', flush=True)

    state.close()

    n_cells = state.profiles.shape[2]
    discovered_tc = np.full((n_frames, n_cells), np.nan)
    for ff, activity in enumerate(activity_by_frame):
        for cell_id, value in activity.items():
            discovered_tc[ff, cell_id] = value

    total_elapsed = time.time() - t0
    return state, discovered_tc, total_elapsed


def save_mirroring_mat(path, profiles, discovered_tc, first_detected_frame):
    """Writes P/T datasets in the same raw HDF5 storage layout demoData1.mat
    uses (see matlab_io.py's module docstring: MATLAB stores column-major,
    so h5py sees every axis reversed) -- (2,1,0) and (1,0) are each their own
    inverse, so applying the same transpose the loaders use to go raw->Python
    also goes Python->raw. Readable back via the existing
    load_profiles(f)/load_time_courses(f, 'T') helpers, unmodified."""
    raw_p = np.transpose(profiles, (2, 1, 0))            # (nCells, movX, movY)
    raw_t = np.transpose(discovered_tc, (1, 0))           # (nCells, movF)

    with h5py.File(path, 'w') as f:
        f.create_dataset('P', data=raw_p)
        f.create_dataset('T', data=raw_t)
        f.create_dataset('firstDetectedFrame', data=np.array(
            [first_detected_frame.get(i, -1) if first_detected_frame.get(i) is not None else -1
             for i in range(profiles.shape[2])], dtype=np.float64))


def main():
    with h5py.File(DEMO_PATH, 'r') as f:
        cnmf_profiles = load_profiles(f)
        movie = Hdf5Movie(f['M'])

        n_frames = N_FRAMES if N_FRAMES is not None else movie.shape[2]
        print(f'movie: {movie.shape}, CNMF cells: {cnmf_profiles.shape[2]}, '
              f'processing {n_frames} frames (full movie)\n', flush=True)

        state, discovered_tc, total_elapsed = run_discovery(movie, n_frames)
        n_disc = state.profiles.shape[2]

        print(f'\ndiscovered {n_disc} cell(s) in {n_frames} frames (CNMF has {cnmf_profiles.shape[2]})')
        print(f'total time: {total_elapsed:.1f}s, {total_elapsed / n_frames * 1000:.2f} ms/frame average')

        save_mirroring_mat(OUT_PATH, state.profiles, discovered_tc, state.first_detected_frame)
        print(f'\nsaved discovered profiles/traces to {OUT_PATH}')


if __name__ == '__main__':
    main()
