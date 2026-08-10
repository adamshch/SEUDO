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
# blob_spacing=3 default. Temporal smoothing comes from lookahead_frames
# (a FitParams field -- see its docstring), which superseded the earlier
# ds_time=3 (causal trailing average) after a direct real-data A/B test
# found the forward-looking version measurably better (34/52 vs. 32/52
# matched CNMF cells over the same 3000-frame window at lookahead_frames=3)
# -- ds_time is no longer set here since lookahead_frames takes over the
# averaging entirely whenever it's > 1. This also means realSEUDOfit can
# return None early on and reports on a lagged frame_index -- see
# run_discovery's handling of that below.
#
# lookahead_frames=5 (was 3), user-directed, kept despite a full-movie A/B
# showing a real precision cost at these merge-threshold settings: matched
# stayed at exactly 51 (identical to lookahead=3's full-movie result) while
# unmatched climbed 27->35 -- zero extra real cells found, 8 more false
# positives, over the whole 41756-frame movie (a 3000-frame sweep alone had
# looked like a clean win, 29/29/0 vs. lookahead=3's 24/24/0 -- see
# streaming.py's project memory for the full writeup of why the short
# window didn't generalize). User elected to keep 5 and tune
# eq8/eq9_merge_threshold further to recover precision instead of
# reverting -- see PROMOTION_PARAMS below for the current sweep result.
#
# spatial_denoise_radius: REMOVED (was 2.0), user-directed. The paper-vs-
# code audit's item 3 originally motivated this (paper preprocesses every
# incoming frame with a spatial Gaussian filter AND a temporal running
# average) and an isolated 3000-frame sweep found radius=2.0 the clear
# standout for matched-cell recall. But building the streaming debug GUI
# (seudo/gui/debug_streaming.py) made the actual failure mode directly
# visible: a 2D Gaussian blur turns raw pixel noise into smooth, spatially-
# correlated local maxima -- exactly the shape the threshold+connected-
# components detector is looking for -- manufacturing ROI-shaped Xtemp
# candidates out of pure noise, not just smoothing real signal. Only the
# existing 3-frame temporal averaging (FitParams.lookahead_frames=3,
# unchanged) remains; no spatial filtering at all now.
SEUDO_PARAMS = dict(sigma2=0.0020, lambda_blob=10.0, blob_radius=3.0, pad_space=5,
                     n_jobs=8, native_nthreads=4, blob_spacing=3.0,
                     spatial_denoise_radius=None, lookahead_frames=5)

# real-data one-at-a-time sweep (scripts/benchmark_detection_params.py) found
# noise_grid_shape=(2,2) (was (1,1)) and consecutive_frames_required=3 (was
# 5, see PROMOTION_PARAMS below) measurably improve recall on a 3000-frame
# window with no/minimal false-positive cost: noise_grid_shape adapts to
# spatial noise variation the same way already validated on synthetic data;
# consecutive_frames_required=3 is both faster to promote AND more forgiving
# of a real cell that doesn't sustain 5 straight causal hits.
#
# exclude_radius_known_cells: fully disabled (-1), user-directed. This
# dataset has many genuinely overlapping cells, and exclude_radius_known_cells
# is not paper/rois_params.m-derived (confirmed by direct grep -- invented
# for this port). Even r=0 (tried first) still unconditionally excludes a
# known cell's own footprint pixels from ever being used by a new candidate
# -- not actually "off" -- so a second, real cell sharing pixels with an
# already-known one could never be detected. A full-movie r=0 benchmark
# still showed heavy fragmentation (66 discovered/43 matched/23 unmatched)
# and a fully-disabled sweep recovered far more real recall (up to 50/52
# matched) at a real, substantial precision cost (up to 138 unmatched) --
# not a clean win, but the user's explicit direction is to prioritize not
# structurally blocking overlapping-cell detection over precision here.
# See DetectionParams.exclude_radius_known_cells's own docstring for the
# r<0 "real off" semantics this required adding to streaming.py.
#
# A first full-movie attempt at exclude_radius_known_cells=-1 with
# eq8/eq9_merge_threshold left at their 0.75 default was aborted early:
# with no spatial exclusion protecting a promoted cell's own territory,
# nothing stopped the SAME region from re-spawning a fresh Xtemp candidate
# on essentially every subsequent frame, and 0.75 (containment-ratio-based,
# only triggers on real pixel OVERLAP) turned out not to catch most of
# these -- 175 "cells" discovered by 32% through the movie, still climbing,
# with per-frame cost growing right along with cell count (each fit gets
# more expensive as state.profiles grows), compounding into a runaway: ETA
# was INCREASING call over call, not shrinking. eq8_merge_threshold/
# eq9_merge_threshold lowered to 0.4 (from 0.75) as the intended
# replacement defense per PromotionParams' own docstring ("these merge
# checks become the primary remaining defense instead of a secondary
# safety net" when exclusion is reduced/disabled) -- NOTE an earlier
# threshold sweep in this session found lowering these thresholds a null
# result specifically at exclude_radius=0, because THAT config's duplicates
# were spatially DISJOINT (no overlap for Eq. 8/9 to ever catch) -- but
# exclude_radius=-1 (this config) is a different regime: a cell's own
# footprint is no longer excluded at all, so a re-spawned duplicate now
# lands squarely ON TOP of the original (maximal overlap), which IS what
# Eq. 8/9's containment-ratio check is designed to catch. Verify on a
# 3000-frame run before trusting this at full-movie scale again.
#
# candidate_profile_threshold: zeroes out any pixel below this fraction of
# a promoted candidate's own peak (see DetectionParams' docstring) -- cleans
# up the jagged, low-amplitude mask-union fringe a wobbling raw detection
# leaves behind (confirmed directly: removed a visible staircase artifact
# on a known illumination-artifact "cell"). Originally tuned to 0.3 via a
# 3000-frame sweep while spatial_denoise_radius=2.0 was still active (33/52
# matched, 4 unmatched vs. baseline's 32/52 matched, 5 unmatched) -- but
# after removing spatial denoising, 0.3 was found to chop out a lot of
# genuinely-real interior/edge pixels (denoising used to smooth those above
# the threshold; without it, raw per-frame amplitude is noisier and more of
# a real cell's true footprint legitimately dips below any fixed fraction
# of peak on a given confirmation frame), producing visibly holey/sparse
# ROI shapes. Lowered to 0.05, which visibly restored solid, contiguous
# shapes in the debug GUI without changing cell/track counts (this
# parameter only reshapes what gets stored, it doesn't affect detection or
# promotion decisions). Requires the exclude-mask fix in _promote_candidate
# (footprint_source) to avoid regressing precision -- thresholding the
# profile used for R2/exclude-mask purposes too (tried first) measurably
# worsened results.
#
# xtemp_smooth_sigma: purely a shape-quality request (make discovered ROI
# shapes visually smoother/rounder, not a detection-accuracy change) -- see
# DetectionParams.xtemp_smooth_sigma's docstring for the mechanism and why
# it's safe (can't manufacture new candidates, unlike the removed
# spatial_denoise_radius). First tried at 1.0: a 12000-frame visual
# comparison (scripts/output/realseudo_rois_grid_{no_smooth,smooth1}_
# 12000frames.png) confirmed a clear qualitative win, shapes went from
# jagged/speckled to smooth and round, confirmed again at full-movie scale
# (scripts/output/*_xtempsmooth1.png) -- but cost more recall than
# expected at full scale (46/52 matched -> 41/52 matched vs. no smoothing).
# User asked for something more subtle: lowered to 0.5. Visual check
# (scripts/output/realseudo_rois_grid_sigma05_merge03_12000frames.png)
# confirms it's genuinely more subtle -- still visibly cleaner than no
# smoothing, but retains more texture/graininess than 1.0's fully-rounded
# look, as intended.
DETECTION_PARAMS = dict(exclude_radius_known_cells=-1, noise_grid_shape=(2, 2),
                         candidate_profile_threshold=0.05, xtemp_smooth_sigma=0.5)
# eq8/eq9_merge_threshold: tried 0.3 (from 0.2) after lowering
# xtemp_smooth_sigma to 0.5 -- smoothing changes the profiles Eq.8/9's
# containment-ratio checks operate on, so the optimal threshold isn't
# guaranteed to carry over from one smoothing level to another. A 12000-
# frame sweep (0.15, 0.2, 0.25, 0.3, 0.4) at sigma=0.5 suggested the
# optimum DID shift: 0.3 strictly dominated the old 0.2 (37/52 matched,
# 3 unmatched vs. 34/52 matched, 4 unmatched in that window) and dominated
# 0.25 too (same precision, more recall); 0.4 finds more cells but costs
# more false positives. DID NOT hold at full-movie scale, though: 0.3 gave
# 56 discovered/46 matched/10 unmatched -- visibly more fragmented than
# either sigma=1.0/merge=0.2 (46/41/5) or no-smoothing/merge=0.2 (52/46/6)
# in the contour overlay (scripts/output/realseudo_vs_cnmf_contours_
# sigma05_merge03.png) -- same "short window doesn't reliably predict
# full-movie precision" lesson as lookahead_frames=5 earlier this session.
# Reverted to 0.2 pending a full-movie confirmation of sigma=0.5/merge=0.2
# specifically (never directly tested -- 0.2 was only validated at
# sigma=1.0, 0.3 was only validated in the 12000-frame window), to isolate
# whether 0.5 itself or the 0.3 retune is the actual regression.
PROMOTION_PARAMS = dict(consecutive_frames_required=3, eq8_merge_threshold=0.2, eq9_merge_threshold=0.2)


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
