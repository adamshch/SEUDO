"""Run the streaming ROI-discovery pipeline (realSEUDOfit) on the real demo
movie (code/demoData1.mat) starting from NO prior knowledge of cell
locations, and compare what it discovers against the CNMF-derived ground
truth already bundled in that file (profiles P, time courses T).

This is a correctness/quality check for the streaming pipeline, not a claim
that either algorithm is "better" -- CNMF jointly optimizes spatial and
temporal factors offline over the whole movie, while realSEUDOfit discovers
ROIs online, causally, one frame at a time, with no lookahead. Reports, for
each CNMF cell: whether realSEUDO found a matching ROI (by centroid
proximity), how spatially similar the two profiles are (correlation over a
shared local window, not raw pixel values -- avoids the movie's known raw
pixel-calibration mismatch against the demo's T/P scale), and how well
their time courses correlate over the frames actually processed.

Deliberately does NOT open the classification GUI -- per direction, that
stays a separate, later step. Instead this script saves the discovered
profiles + traces to an .npz file so that later step can load them into a
Seudo object without re-running discovery (which can take a while -- see
the per-frame timing this script prints).
"""

import os
import sys
import time

import h5py
import numpy as np

sys.path.insert(0, '.')
from seudo.geometry import compute_roi_coms  # noqa: E402
from seudo.matlab_io import Hdf5Movie, load_profiles, load_time_courses  # noqa: E402
from seudo.stats import correlation_vector_matrix  # noqa: E402
from seudo.streaming import FitParams, StreamingState, realSEUDOfit  # noqa: E402

DEMO_PATH = '../code/demoData1.mat'
OUT_PATH = 'scripts/output/realseudo_discovery.npz'

N_FRAMES = 8000    # how many frames of the movie to process; None = the whole movie
PROGRESS_EVERY = 200
MATCH_MAX_DIST = 10.0  # px, centroid-distance cutoff for calling two ROIs "the same cell"

# same demo.m-derived defaults used elsewhere in this project for this dataset,
# plus the best (n_jobs, native_nthreads) combo found by benchmark_native_nthreads.py
SEUDO_PARAMS = dict(sigma2=0.0020, lambda_blob=10.0, blob_radius=3.0, pad_space=5,
                     n_jobs=8, native_nthreads=4)


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
                  f'(total so far: {state.profiles.shape[2]})')

        if (ff + 1) % PROGRESS_EVERY == 0:
            elapsed = time.time() - t0
            rate = elapsed / (ff + 1)
            print(f'  ...{ff + 1}/{n_frames} frames, {rate * 1000:.1f} ms/frame, '
                  f'{state.profiles.shape[2]} cells so far, '
                  f'{rate * n_frames:.0f}s estimated for this range')

    state.close()  # release the per-cell fitting thread pool (n_jobs > 1)

    n_cells = state.profiles.shape[2]
    discovered_tc = np.full((n_frames, n_cells), np.nan)
    for ff, activity in enumerate(activity_by_frame):
        for cell_id, value in activity.items():
            discovered_tc[ff, cell_id] = value

    return state, discovered_tc


def match_rois(cnmf_profiles, discovered_profiles, max_dist):
    """Greedy nearest-centroid matching, closest pairs first."""
    cnmf_coms, _ = compute_roi_coms(cnmf_profiles)
    disc_coms, _ = compute_roi_coms(discovered_profiles)

    pairs = []
    for i in range(cnmf_profiles.shape[2]):
        if np.any(np.isnan(cnmf_coms[i])):
            continue
        for j in range(discovered_profiles.shape[2]):
            if np.any(np.isnan(disc_coms[j])):
                continue
            d = float(np.hypot(*(cnmf_coms[i] - disc_coms[j])))
            if d <= max_dist:
                pairs.append((d, i, j))
    pairs.sort(key=lambda p: p[0])

    matched_cnmf, matched_disc = {}, set()
    for d, i, j in pairs:
        if i in matched_cnmf or j in matched_disc:
            continue
        matched_cnmf[i] = (j, d)
        matched_disc.add(j)

    return matched_cnmf, matched_disc


def profile_correlation(cnmf_profiles, discovered_profiles, cnmf_idx, disc_idx, pad=10):
    """Pearson correlation of the two profiles, restricted to a padded
    window around their combined footprint -- correlating over the whole
    frame would be dominated by the huge shared all-zero background."""
    _, cnmf_bounds = compute_roi_coms(cnmf_profiles[:, :, cnmf_idx])
    _, disc_bounds = compute_roi_coms(discovered_profiles[:, :, disc_idx])
    mov_y, mov_x = cnmf_profiles.shape[:2]

    y0 = max(0, int(min(cnmf_bounds[0, 0], disc_bounds[0, 0])) - pad)
    y1 = min(mov_y - 1, int(max(cnmf_bounds[0, 1], disc_bounds[0, 1])) + pad)
    x0 = max(0, int(min(cnmf_bounds[0, 2], disc_bounds[0, 2])) - pad)
    x1 = min(mov_x - 1, int(max(cnmf_bounds[0, 3], disc_bounds[0, 3])) + pad)

    a = cnmf_profiles[y0:y1 + 1, x0:x1 + 1, cnmf_idx].ravel()
    b = discovered_profiles[y0:y1 + 1, x0:x1 + 1, disc_idx].ravel()
    return float(correlation_vector_matrix(a, b[:, np.newaxis])[0])


def trace_correlation(T, discovered_tc, cnmf_idx, disc_idx, first_detected_frame):
    """Correlate over the frames the discovered cell actually has a value
    for (from its promotion frame onward) -- comparing against the NaN
    pre-promotion gap would be meaningless."""
    start = first_detected_frame.get(disc_idx) or 0
    n = discovered_tc.shape[0]
    t_true = T[start:n, cnmf_idx]
    t_disc = discovered_tc[start:, disc_idx]
    valid = ~np.isnan(t_disc)
    if valid.sum() < 5:
        return np.nan
    return float(np.corrcoef(t_true[valid], t_disc[valid])[0, 1])


def main():
    with h5py.File(DEMO_PATH, 'r') as f:
        profiles = load_profiles(f)
        T = load_time_courses(f, 'T')
        movie = Hdf5Movie(f['M'])

        n_frames = N_FRAMES if N_FRAMES is not None else movie.shape[2]
        print(f'movie: {movie.shape}, CNMF cells: {profiles.shape[2]}, processing {n_frames} frames\n')

        state, discovered_tc = run_discovery(movie, n_frames)
        n_disc = state.profiles.shape[2]
        print(f'\ndiscovered {n_disc} cell(s) in {n_frames} frames (CNMF has {profiles.shape[2]})')

        matched_cnmf, matched_disc = match_rois(profiles, state.profiles, MATCH_MAX_DIST)

        print(f'\n{"cnmf_id":>7} {"disc_id":>7} {"dist(px)":>9} {"prof_corr":>9} {"trace_corr":>10}')
        for cc in range(profiles.shape[2]):
            if cc in matched_cnmf:
                dd, dist = matched_cnmf[cc]
                pc = profile_correlation(profiles, state.profiles, cc, dd)
                tcorr = trace_correlation(T, discovered_tc, cc, dd, state.first_detected_frame)
                print(f'{cc:7d} {dd:7d} {dist:9.2f} {pc:9.3f} {tcorr:10.3f}')
            else:
                print(f'{cc:7d} {"--":>7} {"--":>9} {"--":>9} {"--":>10}   (not found by realSEUDO)')

        unmatched_disc = [dd for dd in range(n_disc) if dd not in matched_disc]
        if unmatched_disc:
            print(f'\n{len(unmatched_disc)} discovered cell(s) with no CNMF match '
                  f'(possible false positives, or real cells CNMF missed): {unmatched_disc}')

        print(f'\nsummary: {len(matched_cnmf)}/{profiles.shape[2]} CNMF cells matched, '
              f'{len(unmatched_disc)}/{n_disc} discovered cells unmatched')

        os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
        np.savez(
            OUT_PATH,
            profiles=state.profiles,
            activity=discovered_tc,
            first_detected_frame=np.array(
                [state.first_detected_frame.get(i) for i in range(n_disc)], dtype=float),
            n_frames_processed=n_frames,
        )
        print(f'\nsaved discovered profiles/traces to {OUT_PATH}')


if __name__ == '__main__':
    main()
