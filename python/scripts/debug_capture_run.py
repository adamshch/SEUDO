"""Runs realSEUDOfit over a prefix of the real demo movie, capturing
everything scripts/open_debug_streaming_gui.py's three-panel GUI needs to
inspect what's actually happening frame by frame:

1. Every denoised input frame (avg_frame -- post lookahead-average + spatial
   denoise, exactly what detection/fitting actually see) -- reconstructed
   independently from the raw movie using the same recipe realSEUDOfit
   itself uses internally (verified against the real internal computation
   earlier in this project's development), since avg_frame is a local
   variable inside realSEUDOfit and isn't otherwise exposed.
2. Every Xtemp candidate track EVER created -- its LATEST profile snapshot,
   creation frame, and the last frame it was actually MATCHED (not merely
   tolerated as a gap -- both timestamps come from track.history, the same
   internal, lookahead-lagged clock _update_candidate_tracks itself stamps
   entries with, not the raw movie-frame loop counter, which runs
   ~lookahead_frames-1 frames ahead of it), PLUS its actual fate (promoted /
   merged into an existing cell / dropped via gap), captured by hooking
   _promote_candidate directly rather than inferred after the fact -- a
   track that reaches promotion is popped from state.candidate_tracks and
   promoted within the SAME realSEUDOfit call, so its snapshot here is
   always its state one frame before that, never its final one.
3. Every promoted Xstab cell's final profile and promotion frame (already
   tracked by state.profiles / state.first_detected_frame -- no extra
   instrumentation needed).

Saves everything to a single HDF5 file for the GUI to load and browse
interactively, so the (potentially slow) discovery run only has to happen
once per debugging session.
"""

import sys

import h5py
import numpy as np
from scipy.signal import convolve2d

sys.path.insert(0, '.')
from seudo.blob import make_smoothing_kernel  # noqa: E402
from seudo.matlab_io import Hdf5Movie  # noqa: E402
from seudo import streaming as streaming_module  # noqa: E402
from seudo.streaming import (DetectionParams, FitParams, PromotionParams, StreamingState,  # noqa: E402
                              _build_promoted_profile, realSEUDOfit)

# outcome codes saved as xtemp_outcome_codes -- see seudo/gui/debug_streaming.py
OUTCOME_DROPPED = 0
OUTCOME_PROMOTED = 1
OUTCOME_MERGED = 2

DEMO_PATH = '../code/demoData1.mat'
OUT_PATH = 'scripts/output/debug_capture.h5'

N_FRAMES = 3000  # how many raw input frames to feed in

# same production defaults as scripts/run_realseudo_full_movie.py, EXCEPT
# exclude_radius_known_cells=0 (production default: 2) -- deliberately
# using the more permissive setting here, since that's specifically where
# this session's investigation found heavy Xtemp fragmentation (66
# discovered/43 matched/23 unmatched vs. the radius=2 baseline's
# 37/33/4) -- the whole point of this capture is to see that fragmentation
# directly, not to reproduce the cleaner production run.
#
# spatial_denoise_radius removed entirely (was 2.0) -- user-directed, after
# the GUI made it visible that a 2D Gaussian blur turns raw pixel noise
# into smooth, spatially-correlated local maxima that the threshold+
# connected-components detector reads as genuine small ROIs. Only the
# existing 3-frame temporal averaging (FitParams.lookahead_frames=3,
# unchanged) remains.
SEUDO_PARAMS = dict(sigma2=0.0020, lambda_blob=10.0, blob_radius=3.0, pad_space=5,
                     n_jobs=8, native_nthreads=4, blob_spacing=3.0,
                     spatial_denoise_radius=None)
DETECTION_PARAMS = dict(exclude_radius_known_cells=0, noise_grid_shape=(2, 2),
                         candidate_profile_threshold=0.05)
PROMOTION_PARAMS = dict(consecutive_frames_required=3)


def main():
    with h5py.File(DEMO_PATH, 'r') as f:
        movie = Hdf5Movie(f['M'])
        y, x = movie.shape[:2]
        n_raw = min(N_FRAMES, movie.shape[2])

        print(f'loading {n_raw} raw frames...', flush=True)
        raw_frames = [movie.get_frame(ff) for ff in range(n_raw)]

    lookahead = SEUDO_PARAMS.get('lookahead_frames', 3)
    denoise_radius = SEUDO_PARAMS.get('spatial_denoise_radius')
    kernel = make_smoothing_kernel(denoise_radius) if denoise_radius else None

    def reconstruct_avg_frame(fidx):
        window = raw_frames[fidx:fidx + lookahead]
        avg = np.mean(window, axis=0)
        return convolve2d(avg, kernel, mode='same') if kernel is not None else avg

    state = StreamingState((y, x), fit=FitParams(**SEUDO_PARAMS),
                            detection=DetectionParams(**DETECTION_PARAMS),
                            promotion=PromotionParams(**PROMOTION_PARAMS))
    threshold = state.detection.candidate_profile_threshold

    avg_frame_list = []
    avg_frame_indices = []
    xtemp_snapshots = {}  # track_id -> dict(creation_frame, last_seen_frame, profile)
    track_fates = {}  # track_id -> dict(outcome_code, cell_id, frame)
    object_id_to_track_id = {}  # id(track) -> track_id, resolved for _promote_candidate's spy below

    _orig_promote_candidate = streaming_module._promote_candidate

    def _spy_promote_candidate(state, track, frame_index):
        tid = object_id_to_track_id.get(id(track))
        cell_id, is_new = _orig_promote_candidate(state, track, frame_index)
        if tid is not None:
            track_fates[tid] = dict(outcome_code=OUTCOME_PROMOTED if is_new else OUTCOME_MERGED,
                                     cell_id=cell_id, frame=frame_index)
        return cell_id, is_new

    streaming_module._promote_candidate = _spy_promote_candidate

    n_processable = n_raw - lookahead + 1
    print(f'running realSEUDOfit over {n_processable} processable frames...', flush=True)
    for ff in range(n_raw):
        # snapshot object identity -> track_id BEFORE this frame's
        # processing, so the _promote_candidate spy above can resolve
        # which track_id is being promoted/merged this frame -- promotion
        # pops and processes the bare track object within this SAME
        # realSEUDOfit call, before we'd otherwise get a chance to see it
        for tid, track in state.candidate_tracks.items():
            object_id_to_track_id[id(track)] = tid

        result = realSEUDOfit(raw_frames[ff], state)

        for tid, track in state.candidate_tracks.items():
            # both timestamps read from track.history (the internal,
            # lookahead-lagged frame clock _update_candidate_tracks itself
            # stamps entries with) -- NOT the raw loop counter `ff`, which
            # runs ~lookahead_frames-1 frames ahead of it and would silently
            # mix two different clocks together
            creation_frame = track.history[0][0] if track.history else ff
            last_seen_frame = track.history[-1][0] if track.history else ff
            profile, _bbox = _build_promoted_profile(track, (y, x), threshold)
            # match_count: the actual number of confirmed detections
            # (track.consecutive_frames), stored explicitly rather than
            # inferred from last_seen_frame - creation_frame -- a tolerated
            # gap (PromotionParams.max_track_gap) between two matches
            # inflates that frame-index difference without adding a match,
            # so the difference alone can't distinguish "3 consecutive
            # matches, one frame short of promotion" from "2 matches with
            # a gap between them."
            xtemp_snapshots[tid] = dict(creation_frame=creation_frame, last_seen_frame=last_seen_frame,
                                         match_count=track.consecutive_frames, profile=profile)

        if result is not None:
            fidx = result.frame_index
            avg_frame_list.append(reconstruct_avg_frame(fidx))
            avg_frame_indices.append(fidx)

        if (ff + 1) % 500 == 0:
            print(f'  ...{ff + 1}/{n_raw} frames, {len(state.candidate_tracks)} active Xtemp tracks, '
                  f'{state.profiles.shape[2]} Xstab cells so far', flush=True)

    state.close()
    streaming_module._promote_candidate = _orig_promote_candidate

    n_xstab = state.profiles.shape[2]
    n_promoted = sum(1 for fate in track_fates.values() if fate['outcome_code'] == OUTCOME_PROMOTED)
    n_merged = sum(1 for fate in track_fates.values() if fate['outcome_code'] == OUTCOME_MERGED)
    n_dropped = len(xtemp_snapshots) - len(track_fates)
    print(f'\ndone: {n_xstab} Xstab cells, {len(xtemp_snapshots)} Xtemp tracks ever created '
          f'({n_promoted} promoted, {n_merged} merged, {n_dropped} dropped), '
          f'{len(avg_frame_list)} denoised frames captured')

    xtemp_ids = sorted(xtemp_snapshots.keys())
    xtemp_profiles = np.stack([xtemp_snapshots[tid]['profile'] for tid in xtemp_ids], axis=0) if xtemp_ids \
        else np.zeros((0, y, x))
    xtemp_creation = np.array([xtemp_snapshots[tid]['creation_frame'] for tid in xtemp_ids], dtype=np.int64)
    xtemp_last_seen = np.array([xtemp_snapshots[tid]['last_seen_frame'] for tid in xtemp_ids], dtype=np.int64)
    xtemp_match_counts = np.array([xtemp_snapshots[tid]['match_count'] for tid in xtemp_ids], dtype=np.int64)
    xtemp_outcome_codes = np.array(
        [track_fates[tid]['outcome_code'] if tid in track_fates else OUTCOME_DROPPED for tid in xtemp_ids],
        dtype=np.int64)
    xtemp_outcome_cell_ids = np.array(
        [track_fates[tid]['cell_id'] if tid in track_fates else -1 for tid in xtemp_ids], dtype=np.int64)
    xtemp_outcome_frames = np.array(
        [track_fates[tid]['frame'] if tid in track_fates else -1 for tid in xtemp_ids], dtype=np.int64)

    xstab_promotion = np.array(
        [state.first_detected_frame.get(i, -1) if state.first_detected_frame.get(i) is not None else -1
         for i in range(n_xstab)], dtype=np.int64)

    with h5py.File(OUT_PATH, 'w') as f:
        f.create_dataset('avg_frames', data=np.stack(avg_frame_list, axis=0) if avg_frame_list
                          else np.zeros((0, y, x)))
        f.create_dataset('avg_frame_indices', data=np.array(avg_frame_indices, dtype=np.int64))
        f.create_dataset('xtemp_profiles', data=xtemp_profiles)
        f.create_dataset('xtemp_track_ids', data=np.array(xtemp_ids, dtype=np.int64))
        f.create_dataset('xtemp_creation_frames', data=xtemp_creation)
        f.create_dataset('xtemp_last_seen_frames', data=xtemp_last_seen)
        f.create_dataset('xtemp_match_counts', data=xtemp_match_counts)
        f.create_dataset('xtemp_outcome_codes', data=xtemp_outcome_codes)
        f.create_dataset('xtemp_outcome_cell_ids', data=xtemp_outcome_cell_ids)
        f.create_dataset('xtemp_outcome_frames', data=xtemp_outcome_frames)
        f.create_dataset('xstab_profiles', data=state.profiles)
        f.create_dataset('xstab_promotion_frames', data=xstab_promotion)

    print(f'saved to {OUT_PATH}')


if __name__ == '__main__':
    main()
