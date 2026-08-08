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
from seudo.streaming import DetectionParams, FitParams, PromotionParams, StreamingState, realSEUDOfit  # noqa: E402

DEMO_PATH = '../code/demoData1.mat'
OUT_PATH = '../code/demoDataRealSEUDO.mat'

N_FRAMES = None    # None = the whole movie
PROGRESS_EVERY = 500

# same demo.m-derived defaults used elsewhere in this project for this
# dataset, plus the validated n_jobs/native_nthreads combo and the
# blob_spacing=3 default. Temporal smoothing now comes from
# lookahead_frames=3 (a FitParams default -- see its docstring), which
# superseded the earlier ds_time=3 (causal trailing average) after a
# direct real-data A/B test found the forward-looking version measurably
# better (34/52 vs. 32/52 matched CNMF cells over the same 3000-frame
# window) -- ds_time is no longer set here since lookahead_frames takes
# over the averaging entirely whenever it's > 1. This also means
# realSEUDOfit can return None early on and reports on a lagged
# frame_index -- see run_discovery's handling of that below.
#
# spatial_denoise_radius=2.0: the paper-vs-code audit's item 3 -- the paper
# preprocesses every incoming frame with a spatial Gaussian filter AND a
# temporal running average before any fitting, we only had the temporal
# half until now. Real-data sweep on a 3000-frame window (None/1.0/1.5/2.0/
# 3.0) found radius=2.0 the clear standout: matched CNMF cells jumped
# 20->32 (+60%) for only 1->3 new false positives -- by far the single
# biggest recall improvement found in this whole audit. radius=3.0 is
# dramatically slower (165ms/frame vs 35ms) for worse recall, so not
# just "bigger is better" -- this value was chosen from real measurement,
# not guessed.
SEUDO_PARAMS = dict(sigma2=0.0020, lambda_blob=10.0, blob_radius=3.0, pad_space=5,
                     n_jobs=8, native_nthreads=4, blob_spacing=3.0,
                     spatial_denoise_radius=2.0)

# real-data one-at-a-time sweep (scripts/benchmark_detection_params.py) found
# these three -- none with a reference value from the paper or rois_params.m,
# unlike min_avg_px/blobify_radius -- measurably improve recall on a
# 3000-frame window with no/minimal false-positive cost: exclude_radius=2
# (was 5) finds cells packed closer to already-known ones; noise_grid_shape
# =(2,2) (was (1,1)) adapts to spatial noise variation the same way already
# validated on synthetic data; consecutive_frames_required=3 (was 5) is both
# faster to promote AND more forgiving of a real cell that doesn't sustain
# 5 straight causal hits. Combined: 20/21 matched vs. baseline's 16/16 in the
# same window (one new false positive) -- still net better than any pairwise
# subset of the three (see the benchmark script's combined-sweep follow-up).
DETECTION_PARAMS = dict(exclude_radius_known_cells=2, noise_grid_shape=(2, 2))
PROMOTION_PARAMS = dict(consecutive_frames_required=3)


def run_discovery(movie, n_frames):
    state = StreamingState(
        movie.shape[:2], fit=FitParams(**SEUDO_PARAMS),
        detection=DetectionParams(**DETECTION_PARAMS), promotion=PromotionParams(**PROMOTION_PARAMS),
    )
    # keyed by result.frame_index, not loop position -- with
    # FitParams.lookahead_frames > 1 (the default), realSEUDOfit returns
    # None for the first lookahead_frames-1 calls (not enough future
    # context yet) and every reported frame_index lags the raw input frame
    # ff by lookahead_frames-1, so the two are no longer interchangeable.
    activity_by_frame = {}

    t0 = time.time()
    for ff in range(n_frames):
        frame = movie.get_frame(ff)
        result = realSEUDOfit(frame, state)
        if result is None:
            continue
        activity_by_frame[result.frame_index] = result.activity

        if result.new_cells:
            print(f'  frame {result.frame_index}: discovered cell(s) {result.new_cells} '
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
    for frame_idx, activity in activity_by_frame.items():
        for cell_id, value in activity.items():
            discovered_tc[frame_idx, cell_id] = value

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
