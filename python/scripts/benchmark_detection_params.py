"""One-at-a-time real-data sensitivity sweep over streaming's remaining
UNVALIDATED detection/promotion knobs -- cutoff_multiplier,
exclude_radius_known_cells, max_track_gap, consecutive_frames_required,
noise_grid_shape, match_max_centroid_dist -- none of which have a reference
value from either the realSEUDO paper or rois_params.m (unlike min_avg_px/
blobify_radius, already checked and fixed against rois_params.m's own
defaults). Each value is tried in isolation (all other params held at
current FitParams/DetectionParams/PromotionParams defaults) over the same
5000-frame prefix of the real demo movie, reporting discovered-cell count
and match quality against the CNMF ground truth already bundled there.

3000 frames (matching this session's other quick real-data checks) to get
several discovery events in the window before this dataset's usual early
plateau, without paying full-movie cost for every sweep point.
"""

import sys

import h5py
import numpy as np

sys.path.insert(0, '.')
sys.path.insert(0, 'scripts')
from seudo.matlab_io import Hdf5Movie, load_profiles  # noqa: E402
from seudo.streaming import DetectionParams, FitParams, PromotionParams, StreamingState, realSEUDOfit  # noqa: E402
from compare_realseudo_to_cnmf import match_rois  # noqa: E402

DEMO_PATH = '../code/demoData1.mat'
N_FRAMES = 3000
MATCH_MAX_DIST = 10.0

# ds_time no longer set here -- FitParams.lookahead_frames defaults to 3 and
# takes over temporal smoothing entirely whenever it's > 1, making an
# explicit ds_time=3 alongside it a no-op (see lookahead_frames' docstring).
BASE_FIT = dict(sigma2=0.0020, lambda_blob=10.0, blob_radius=3.0, pad_space=5,
                 n_jobs=8, native_nthreads=4, blob_spacing=3.0)

SWEEPS = [
    ('cutoff_multiplier', 'detection', [2.5, 3.0, 4.0, 5.0]),
    ('exclude_radius_known_cells', 'detection', [2, 5, 8]),
    ('noise_grid_shape', 'detection', [(1, 1), (2, 2), (3, 3)]),
    ('max_track_gap', 'promotion', [1, 2, 4]),
    ('consecutive_frames_required', 'promotion', [3, 5, 8]),
    ('match_max_centroid_dist', 'promotion', [3.0, 5.0, 8.0]),
]


def run_once(movie, cnmf_profiles, detection, promotion):
    state = StreamingState(movie.shape[:2], fit=FitParams(**BASE_FIT), detection=detection, promotion=promotion)
    for ff in range(N_FRAMES):
        realSEUDOfit(movie.get_frame(ff), state)
    state.close()
    n_disc = state.profiles.shape[2]
    matched_cnmf, matched_disc = match_rois(cnmf_profiles, state.profiles, MATCH_MAX_DIST)
    n_unmatched = n_disc - len(matched_disc)
    return n_disc, len(matched_cnmf), n_unmatched


def main():
    with h5py.File(DEMO_PATH, 'r') as f:
        cnmf_profiles = load_profiles(f)
        movie = Hdf5Movie(f['M'])

        print(f'baseline (all defaults), {N_FRAMES} frames:')
        n_disc, n_matched, n_unmatched = run_once(movie, cnmf_profiles, DetectionParams(), PromotionParams())
        print(f'  discovered={n_disc} matched={n_matched} unmatched={n_unmatched}\n')

        for field, group, values in SWEEPS:
            print(f'--- sweeping {group}.{field} ---')
            for v in values:
                kwargs = {field: v}
                detection = DetectionParams(**kwargs) if group == 'detection' else DetectionParams()
                promotion = PromotionParams(**kwargs) if group == 'promotion' else PromotionParams()
                n_disc, n_matched, n_unmatched = run_once(movie, cnmf_profiles, detection, promotion)
                print(f'  {field}={v}: discovered={n_disc} matched={n_matched} unmatched={n_unmatched}')
            print()


if __name__ == '__main__':
    main()
