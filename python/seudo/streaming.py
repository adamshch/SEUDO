"""Streaming/online SEUDO: `realSEUDOfit(frame, state)` fits one frame at a
time against a *growing* set of cell profiles, detecting and promoting new
cells as they appear, instead of requiring the whole movie and a fixed,
pre-known cell set upfront (see estimate.py).

This is a fresh design informed by (but not a port of) the sibling realSEUDO
C++/MATLAB repo's real-time ideas -- its actual per-frame orchestration loop
turned out to be dead/prototype code with a load-bearing function
(find_still_rois_in, the candidate-blob detector) never defined anywhere in
that repo. Per explicit user direction this implementation is: strictly
causal (no forward-looking lookahead buffer, unlike realSEUDO's avg_frames),
low-latency (a short consecutive-frames promotion window, not realSEUDO's
up-to-10-frame one), a simple per-frame activity-value return contract (no
event-journal/replay subsystem), and spatially tiled for large fields of
view (a fresh per-frame design -- realSEUDO's patch-splitting was offline,
whole-movie-in-memory only).

Reuses the offline solver's actual fitting internals directly rather than
re-deriving the math: _setup_cell_window/_solve_one_frame_cell (estimate.py),
make_seudo_blob (blob.py), compute_roi_coms (geometry.py).

Candidate detection is two-stage, per the realSEUDO paper's design: raw
blobs are first detected on R1 (the residual after known-cell fitting) and
used to advance already-tracked candidates by centroid matching
(_update_candidate_tracks); each active track's own evolving profile is then
fit against R1 with a real second SEUDO solve and subtracted out
(_subtract_tracked_contributions) to build R2, and only R2 -- genuinely
unexplained by both known cells and tracked candidates -- gets scanned for
brand-new ones (_start_candidate_tracks).
"""

import math
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator
from scipy.signal import convolve2d

from .blob import make_seudo_blob
from .estimate import _setup_cell_window, _solve_one_frame_cell
from .geometry import compute_roi_coms
from ._native import NATIVE_AVAILABLE as _NATIVE_AVAILABLE


@dataclass
class TilingConfig:
    """tile_shape=None (default): one tile spanning the whole frame -- no
    behavior change from an untiled design. Tiling only ever scopes the
    *candidate-detection* pass (see _tile_threshold_mask/_detect_in_tile);
    known-cell fitting always uses each cell's own window regardless of
    tiling, since _setup_cell_window already builds that window independent
    of any grid."""
    tile_shape: tuple = None  # (tile_ht, tile_wd) or None
    overlap: int = 15         # halo pixels borrowed from neighboring tiles on each side

    def build_tiles(self, mov_y, mov_x):
        if self.tile_shape is None:
            core = (0, mov_y - 1, 0, mov_x - 1)
            return [Tile(core=core, halo=core)]

        tile_ht, tile_wd = self.tile_shape
        tiles = []
        for y0 in range(0, mov_y, tile_ht):
            y1 = min(mov_y, y0 + tile_ht) - 1
            for x0 in range(0, mov_x, tile_wd):
                x1 = min(mov_x, x0 + tile_wd) - 1
                core = (y0, y1, x0, x1)
                halo = (
                    max(0, y0 - self.overlap), min(mov_y - 1, y1 + self.overlap),
                    max(0, x0 - self.overlap), min(mov_x - 1, x1 + self.overlap),
                )
                tiles.append(Tile(core=core, halo=halo))
        return tiles


@dataclass
class Tile:
    core: tuple  # (y0, y1, x0, x1), inclusive, global frame coords
    halo: tuple  # (y0, y1, x0, x1), inclusive, core expanded by overlap and clamped


@dataclass
class DetectionParams:
    """Field names/semantics mirror realSEUDO's rois_params.m knobs where
    sensible (min_roi_size, min_avg_px, mask_blur_rad) -- the closest
    available spec for its candidate detector, which was never implemented
    in that repo. This is fresh code, not a port.

    noise_grid_shape: (n_sections_y, n_sections_x) for spatially-varying
    noise estimation -- splits the frame into this many sections, estimates
    the noise floor locally in each, then smoothly interpolates between
    section centers to build a full-resolution noise map, instead of one
    global scalar. Guards against a single region with unusually high or
    low background (uneven illumination, neuropil) skewing detection
    everywhere else in the frame -- per the realSEUDO paper's own
    noise-estimation approach ("kriging procedure or a cheap approximation
    to a local median evaluated independently for each pixel"). Default
    (1, 1) is exactly the old single-scalar behavior (bit-identical)."""
    cutoff_multiplier: float = 4.0
    mask_blur_rad: int = 1
    min_roi_size: int = 50
    min_avg_px: float = -2.0  # <0: |min_avg_px| * local noise level; >=0: absolute threshold
    exclude_radius_known_cells: int = 5
    noise_grid_shape: tuple = (1, 1)


@dataclass
class PromotionParams:
    """consecutive_frames_required=5 is deliberately tight relative to
    realSEUDO's design (an avg_frames=5 forward-looking buffer *before* up
    to recent_frames=10 more frames of promotion latency) -- per user
    direction for causal, low-latency detection. max_track_gap=1 tolerates
    one missed frame without resetting a track's progress."""
    consecutive_frames_required: int = 5
    max_track_gap: int = 1
    match_max_centroid_dist: float = 5.0


@dataclass
class FitParams:
    """Same knobs as estimate_time_courses_with_seudo -- no ds_time here,
    since streaming is strictly causal with no frame-averaging buffer
    (equivalent to ds_time=1).

    n_jobs: number of known cells to fit concurrently within a single frame,
    in a ThreadPoolExecutor (same lever, same rationale, as
    estimate_time_courses_with_seudo's own n_jobs -- each cell's fit is
    already fully independent, and the native solver releases the GIL for
    its whole call, so real multi-core speedup happens even under Python
    threads). n_jobs=1 (default) is sequential, bit-identical to n_jobs>1 --
    only wall-clock time changes. Kept separate from native_nthreads (which
    parallelizes *within* one cell's own solve) to avoid oversubscription.

    blob_spacing: only affects the native solver (the pure-Python fallback
    always places one blob per pixel) -- see estimate_time_courses_with_seudo's
    docstring for the speed/quality tradeoff (realSEUDO paper finding: a
    coarser blob grid gives more-than-quadratic FISTA speedup with little
    quality loss). Benchmark on real data before relying on it; see
    scripts/benchmark_blob_spacing.py."""
    p: float = 1e-5
    sigma2: float = 0.01
    lambda_blob: float = 20.0
    blob_radius: float = 1.2
    lambda_prof: float = 0.0
    min_pix_for_inclusion: int = 1
    pad_space: int = 10
    use_com: bool = False
    solver_tol: float = 0.01
    solver_max_iter: int = 1000
    use_native: object = 'auto'
    native_l_mode: int = 2
    native_nthreads: int = 1
    blob_spacing: float = 1.0
    n_jobs: int = 1


@dataclass
class FrameResult:
    frame_index: int
    activity: dict       # {cell_id: float}, every currently-known cell, including any promoted this frame
    new_cells: list       # cell_ids promoted this frame (subset of activity.keys())


@dataclass
class RawCandidate:
    centroid: tuple    # (y, x), float, global frame coords
    bbox: tuple         # (y0, y1, x0, x1), inclusive, global frame coords
    mask: np.ndarray     # local boolean mask, shape matches bbox


@dataclass
class CandidateTrack:
    centroid: tuple
    consecutive_frames: int
    gap: int
    history: list = field(default_factory=list)  # [(frame_index, bbox, mask, residual_crop), ...]


def _cell_window_bounds(prof, pad_space, mov_y, mov_x, use_com):
    """Same window-bounds logic as estimate_time_courses_with_seudo's
    per-cell setup (its own precompute loop, not currently factored into a
    reusable function there)."""
    coms, outer_bounds = compute_roi_coms(prof)
    if use_com:
        cy, cx = int(round(coms[0, 1])), int(round(coms[0, 0]))
        y0 = y1 = cy
        x0 = x1 = cx
    else:
        y0, y1, x0, x1 = outer_bounds[0]
    y0 = max(0, int(y0) - pad_space)
    y1 = min(mov_y - 1, int(y1) + pad_space)
    x0 = max(0, int(x0) - pad_space)
    x1 = min(mov_x - 1, int(x1) + pad_space)
    return y0, y1, x0, x1


def _bbox_overlap(a, b):
    ay0, ay1, ax0, ax1 = a
    by0, by1, bx0, bx1 = b
    return ay0 <= by1 and by0 <= ay1 and ax0 <= bx1 and bx0 <= ax1


class StreamingState:
    """Mutable state carried across realSEUDOfit() calls: the (growing) cell
    profile set, a cached per-cell fit setup for each known cell, the
    known-cell exclusion mask for candidate detection, and in-progress
    (not-yet-promoted) candidate tracks. Mutates in place across calls --
    avoids copying a growing profile array every frame."""

    def __init__(
        self, mov_shape, initial_profiles=None, zero_level=0.0,
        tiling=None, detection=None, promotion=None, fit=None,
    ):
        self.mov_y, self.mov_x = int(mov_shape[0]), int(mov_shape[1])
        self.zero_level = zero_level
        self.tiling = tiling or TilingConfig()
        self.detection = detection or DetectionParams()
        self.promotion = promotion or PromotionParams()
        self.fit = fit or FitParams()

        if initial_profiles is None:
            self.profiles = np.zeros((self.mov_y, self.mov_x, 0), dtype=float)
        else:
            self.profiles = np.array(initial_profiles, dtype=float)
            if self.profiles.shape[:2] != (self.mov_y, self.mov_x):
                raise ValueError(
                    f'initial_profiles shape {self.profiles.shape[:2]} does not match '
                    f'mov_shape {(self.mov_y, self.mov_x)}'
                )

        self.use_native = self.fit.use_native
        if self.use_native == 'auto':
            self.use_native = _NATIVE_AVAILABLE
        elif self.use_native and not _NATIVE_AVAILABLE:
            raise RuntimeError(
                'use_native=True but the compiled native SEUDO extension is not available; '
                'build it with seudo/_native/build_native.sh, or pass use_native=False'
            )

        self.one_blob = make_seudo_blob(self.fit.blob_radius)
        self.frame_index = 0
        self.first_detected_frame = {}
        self._cell_setups = {}
        self.known_cell_exclude_mask = np.zeros((self.mov_y, self.mov_x), dtype=bool)
        self.candidate_tracks = {}
        self._next_track_id = 0

        self.tiles = self.tiling.build_tiles(self.mov_y, self.mov_x)

        for cell_id in range(self.profiles.shape[2]):
            self.first_detected_frame[cell_id] = None
            y0, y1, x0, x1 = _cell_window_bounds(
                self.profiles[:, :, cell_id], self.fit.pad_space, self.mov_y, self.mov_x, self.fit.use_com)
            _add_cell_setup(self, cell_id, y0, y1, x0, x1)
            _update_known_cell_exclude_mask(self, cell_id, y0, y1, x0, x1)

        # created lazily, not per-frame -- pool spinup/teardown cost would
        # otherwise repeat on every realSEUDOfit() call
        self._executor = (
            ThreadPoolExecutor(max_workers=self.fit.n_jobs) if self.fit.n_jobs and self.fit.n_jobs > 1 else None
        )

    def close(self):
        """Shut down the per-cell fitting thread pool, if one was created
        (fit.n_jobs > 1). Not required for correctness -- Python's own
        atexit handling joins any live ThreadPoolExecutor at interpreter
        exit -- but call it (or use StreamingState as a context manager)
        to release the worker threads promptly in a long-lived process."""
        if self._executor is not None:
            self._executor.shutdown()
            self._executor = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def _add_cell_setup(state, cell_id, y0, y1, x0, x1):
    setup = _setup_cell_window(
        state.profiles, cell_id, y0, y1, x0, x1, state.fit.min_pix_for_inclusion,
        state.fit.lambda_prof, state.fit.lambda_blob, state.fit.sigma2, state.fit.p, state.one_blob,
    )
    setup['y0'], setup['y1'], setup['x0'], setup['x1'] = y0, y1, x0, x1
    state._cell_setups[cell_id] = setup


def _invalidate_overlapping_setups(state, new_id, y0, y1, x0, x1):
    new_bbox = (y0, y1, x0, x1)
    for cell_id, setup in list(state._cell_setups.items()):
        if cell_id == new_id:
            continue
        cell_bbox = (setup['y0'], setup['y1'], setup['x0'], setup['x1'])
        if _bbox_overlap(new_bbox, cell_bbox):
            _add_cell_setup(state, cell_id, *cell_bbox)


def _update_known_cell_exclude_mask(state, cell_id, y0, y1, x0, x1):
    footprint = state.profiles[y0:y1 + 1, x0:x1 + 1, cell_id] > 0
    r = state.detection.exclude_radius_known_cells
    if r > 0:
        structure = np.ones((2 * r + 1, 2 * r + 1), dtype=bool)
        footprint = ndimage.binary_dilation(footprint, structure=structure)
    state.known_cell_exclude_mask[y0:y1 + 1, x0:x1 + 1] |= footprint


def _estimate_noise_level(residual):
    return float(np.median(residual) - np.min(residual))


def _estimate_noise_map(residual, grid_shape):
    """Spatially-varying noise floor: split the frame into a coarse
    grid_shape=(n_sections_y, n_sections_x) grid, estimate the noise level
    locally in each section (same median-minus-min heuristic as
    _estimate_noise_level, just applied to a crop instead of the whole
    frame), then linearly interpolate between section centers to build a
    full-resolution map -- smooth, not blocky, and avoids one region's
    unusual background skewing detection everywhere else. grid_shape=(1,1)
    (the default) collapses to the old single-scalar behavior exactly,
    broadcast to the frame's shape."""
    mov_y, mov_x = residual.shape
    n_y, n_x = grid_shape

    if n_y == 1 and n_x == 1:
        return np.full((mov_y, mov_x), _estimate_noise_level(residual))

    y_edges = np.linspace(0, mov_y, n_y + 1).astype(int)
    x_edges = np.linspace(0, mov_x, n_x + 1).astype(int)
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2.0
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2.0

    coarse = np.empty((n_y, n_x))
    for i in range(n_y):
        for j in range(n_x):
            section = residual[y_edges[i]:y_edges[i + 1], x_edges[j]:x_edges[j + 1]]
            coarse[i, j] = _estimate_noise_level(section)

    interp = RegularGridInterpolator((y_centers, x_centers), coarse, method='linear',
                                      bounds_error=False, fill_value=None)
    yy, xx = np.mgrid[0:mov_y, 0:mov_x]
    points = np.stack([yy.ravel(), xx.ravel()], axis=-1)
    return interp(points).reshape(mov_y, mov_x)


def _tile_threshold_mask(smoothed, noise_map, tile, exclude_mask, detection):
    """The cheap part of tile detection: crop + threshold + exclude, nothing
    else. Meant to run sequentially over every tile before any thread-pool
    submission -- a plain boolean comparison is fast enough that doing it
    for every tile up front, in this thread, costs less than the overhead
    of farming it out, and it's what lets "blanking" (below) skip a tile
    without ever touching the executor at all.

    Returns (crop, noise_crop, mask) if mask has any True pixel, else None
    ("blanking": nothing here could possibly pass threshold, so the caller
    can skip this tile completely -- correctness-neutral, since the
    expensive dilate+label+filter passes below would find zero candidates
    here too)."""
    hy0, hy1, hx0, hx1 = tile.halo
    crop = smoothed[hy0:hy1 + 1, hx0:hx1 + 1]
    noise_crop = noise_map[hy0:hy1 + 1, hx0:hx1 + 1]
    excl = exclude_mask[hy0:hy1 + 1, hx0:hx1 + 1]

    cutoff = detection.cutoff_multiplier * noise_crop
    mask = (crop > cutoff) & ~excl

    return (crop, noise_crop, mask) if mask.any() else None


def _detect_in_tile(crop, noise_crop, mask, tile, detection):
    """The expensive part of tile detection: dilate + connected-component
    label + per-component size/brightness filtering. Only ever called for a
    tile _tile_threshold_mask already found non-blank -- this is the part
    worth parallelizing across tiles (see realSEUDOfit)."""
    hy0, _hy1, hx0, _hx1 = tile.halo
    cy0, cy1, cx0, cx1 = tile.core

    if detection.mask_blur_rad > 0:
        r = detection.mask_blur_rad
        structure = np.ones((2 * r + 1, 2 * r + 1), dtype=bool)
        mask = ndimage.binary_dilation(mask, structure=structure)

    labeled, n = ndimage.label(mask)

    candidates = []
    for label_id in range(1, n + 1):
        comp_mask = labeled == label_id
        size = int(comp_mask.sum())
        if size < detection.min_roi_size:
            continue
        # local noise level for THIS component, not a single frame-wide
        # scalar -- matters once noise_grid_shape makes noise_crop vary
        local_noise = float(noise_crop[comp_mask].mean())
        avg_thresh = (abs(detection.min_avg_px) * local_noise if detection.min_avg_px < 0
                      else detection.min_avg_px)
        if float(crop[comp_mask].mean()) < avg_thresh:
            continue

        ys, xs = np.nonzero(comp_mask)
        centroid_global = (float(ys.mean()) + hy0, float(xs.mean()) + hx0)

        # keep only if the centroid's pixel falls in THIS tile's core --
        # guarantees exactly-once attribution regardless of grid alignment.
        # Round to a pixel index before comparing against the (integer)
        # core bounds -- comparing the raw float directly would leave a gap
        # in continuous coordinate space between adjacent tiles' bounds
        # (e.g. a centroid of 19.96 satisfies neither "<=19" nor ">=20" and
        # would be dropped by every tile).
        centroid_px = (round(centroid_global[0]), round(centroid_global[1]))
        if not (cy0 <= centroid_px[0] <= cy1 and cx0 <= centroid_px[1] <= cx1):
            continue

        ly0, ly1, lx0, lx1 = int(ys.min()), int(ys.max()), int(xs.min()), int(xs.max())
        bbox_global = (ly0 + hy0, ly1 + hy0, lx0 + hx0, lx1 + hx0)
        local_mask = comp_mask[ly0:ly1 + 1, lx0:lx1 + 1]

        candidates.append(RawCandidate(centroid=centroid_global, bbox=bbox_global, mask=local_mask))

    return candidates


def _run_tile_detection(state, smoothed, noise_map, exclude_mask):
    """Shared phase-1/phase-2 tile scan (see _tile_threshold_mask /
    _detect_in_tile): cheap sequential blanking check over every tile, then
    dilate+label+filter -- in the thread pool when it's worth it -- only for
    tiles found non-blank. Used twice per frame: once against R1 (to find
    candidates for matching against existing tracks) and once against R2
    (to find genuinely new ones) -- see realSEUDOfit."""
    non_blank_tiles = []
    for tile in state.tiles:
        hit = _tile_threshold_mask(smoothed, noise_map, tile, exclude_mask, state.detection)
        if hit is not None:
            non_blank_tiles.append((tile, hit))

    raw_candidates = []
    if state._executor is not None and len(non_blank_tiles) > 1:
        futures = [state._executor.submit(_detect_in_tile, crop, noise_crop, mask, tile, state.detection)
                   for tile, (crop, noise_crop, mask) in non_blank_tiles]
        for fut in futures:
            raw_candidates.extend(fut.result())
    else:
        for tile, (crop, noise_crop, mask) in non_blank_tiles:
            raw_candidates.extend(_detect_in_tile(crop, noise_crop, mask, tile, state.detection))

    return raw_candidates


def _update_candidate_tracks(state, raw_candidates, residual, frame_index):
    """Match this frame's raw candidates (detected on R1, the residual after
    known-cell fitting) against existing tracks by nearest centroid, and
    advance (or drop) each track accordingly -- the same matching logic this
    module has always used. Deliberately does NOT create new tracks for
    unmatched raw candidates: see _start_candidate_tracks, which only ever
    sees raw candidates detected on R2 (R1 with every active track's own
    fitted contribution subtracted out -- see _subtract_tracked_contributions),
    so a track already being followed here can never spawn a duplicate."""
    pairs = []
    for track_id, track in state.candidate_tracks.items():
        for raw_idx, cand in enumerate(raw_candidates):
            d = math.hypot(track.centroid[0] - cand.centroid[0], track.centroid[1] - cand.centroid[1])
            if d <= state.promotion.match_max_centroid_dist:
                pairs.append((d, track_id, raw_idx))
    pairs.sort(key=lambda p: p[0])

    assigned_tracks, assigned_raw = set(), set()
    cap = state.promotion.consecutive_frames_required

    for _d, track_id, raw_idx in pairs:
        if track_id in assigned_tracks or raw_idx in assigned_raw:
            continue
        assigned_tracks.add(track_id)
        assigned_raw.add(raw_idx)

        track = state.candidate_tracks[track_id]
        cand = raw_candidates[raw_idx]
        track.centroid = cand.centroid
        track.consecutive_frames += 1
        track.gap = 0
        y0, y1, x0, x1 = cand.bbox
        track.history.append((frame_index, cand.bbox, cand.mask, residual[y0:y1 + 1, x0:x1 + 1]))
        if len(track.history) > cap:
            track.history.pop(0)

    for track_id in list(state.candidate_tracks.keys()):
        if track_id not in assigned_tracks:
            track = state.candidate_tracks[track_id]
            track.gap += 1
            if track.gap > state.promotion.max_track_gap:
                del state.candidate_tracks[track_id]
            # else: gap tolerated -- consecutive_frames left untouched, not reset


def _subtract_tracked_contributions(state, residual):
    """Two-stage candidate detection, stage 1: fit every currently-active
    track's own evolving profile (built from its accumulated history via
    _build_promoted_profile) against `residual` with a real second SEUDO
    solve -- the same fitting machinery used for known cells, just scoped to
    the track's own small window -- and subtract out whatever it explains.
    This is the realSEUDO paper's two-stage design (fit known cells -> R1 ->
    fit each tracked candidate's own profile against R1 -> R2 -> only R2
    gets scanned for genuinely new blobs): it stops an already-tracked
    candidate's own residual from being mistaken for a second, nearby,
    genuinely-new one in the brand-new-candidate scan below
    (_start_candidate_tracks), on top of (not instead of) the ordinary
    centroid matching that keeps advancing each track's own progress (see
    _update_candidate_tracks).

    Returns (residual2, exclude_mask): exclude_mask is the union of every
    active track's (dilated) footprint regardless of this frame's fitted
    weight -- a track's own territory should never spawn a duplicate new
    track even on a frame its regression fit happens to be weak."""
    residual2 = residual.copy()
    exclude_mask = np.zeros((state.mov_y, state.mov_x), dtype=bool)
    r = state.detection.exclude_radius_known_cells
    structure = np.ones((2 * r + 1, 2 * r + 1), dtype=bool) if r > 0 else None

    for track in state.candidate_tracks.values():
        profile, _history_bbox = _build_promoted_profile(track, (state.mov_y, state.mov_x))
        y0, y1, x0, x1 = _cell_window_bounds(profile, state.fit.pad_space, state.mov_y, state.mov_x, state.fit.use_com)
        footprint = profile[y0:y1 + 1, x0:x1 + 1] > 0

        setup = _setup_cell_window(
            profile[:, :, np.newaxis], 0, y0, y1, x0, x1, state.fit.min_pix_for_inclusion,
            state.fit.lambda_prof, state.fit.lambda_blob, state.fit.sigma2, state.fit.p, state.one_blob,
        )
        this_residual = residual[y0:y1 + 1, x0:x1 + 1].ravel()
        _tc_lsq, fit_fancy, _fit_x, _lsq_cost, _bob_cost = _solve_one_frame_cell(
            this_residual, setup['rois'], setup['rois_scaled'], setup['lambdas'], setup['norm_factors'],
            setup['k1'], setup['k2'], setup['n_y'], setup['n_x'], state.one_blob, setup['operators'],
            state.use_native, state.fit.native_l_mode, state.fit.native_nthreads,
            state.fit.solver_tol, state.fit.solver_max_iter, state.fit.blob_spacing,
        )
        weight = float(fit_fancy[setup['cell_index_within']])
        if weight > 0:
            residual2[y0:y1 + 1, x0:x1 + 1] -= profile[y0:y1 + 1, x0:x1 + 1] * weight

        dilated = ndimage.binary_dilation(footprint, structure=structure) if structure is not None else footprint
        exclude_mask[y0:y1 + 1, x0:x1 + 1] |= dilated

    return residual2, exclude_mask


def _start_candidate_tracks(state, raw_candidates, residual, frame_index):
    """Stage 2: every raw candidate detected on residual2 (see
    _subtract_tracked_contributions) already excludes every currently-active
    track's footprint, so each one is by construction a genuinely new
    candidate -- no matching against existing tracks needed, just start
    tracking it."""
    for cand in raw_candidates:
        track_id = state._next_track_id
        state._next_track_id += 1
        y0, y1, x0, x1 = cand.bbox
        state.candidate_tracks[track_id] = CandidateTrack(
            centroid=cand.centroid, consecutive_frames=1, gap=0,
            history=[(frame_index, cand.bbox, cand.mask, residual[y0:y1 + 1, x0:x1 + 1])],
        )


def _build_promoted_profile(track, mov_shape):
    """Mean of the (nonneg-clipped) residual crops over the track's
    confirmation window, masked to the union of its detected footprints,
    then peak-normalized. Uses the union bbox across history (rather than
    assuming a fixed bbox) so a slowly growing/shifting candidate is handled
    correctly.

    _setup_cell_window L2-normalizes every profile internally (and undoes it
    after solving), so the SOLVE doesn't care about profile scale -- but
    without peak-normalizing here, the built profile's magnitude would bake
    in whatever raw pixel intensity happened to be present during promotion
    (roughly true_amplitude * shape), making the fitted activity coefficient
    read back out as ~1.0 for a steady signal instead of tracking the true
    amplitude. Peak-normalizing to 1.0 matches the "profile is shape only,
    amplitude lives in the fitted coefficient" convention every externally
    supplied profile in this codebase already follows (e.g. test fixtures'
    gaussian_patch), so a promoted cell's activity is directly comparable to
    a pre-known cell's."""
    y0 = min(b[0] for _, b, _, _ in track.history)
    y1 = max(b[1] for _, b, _, _ in track.history)
    x0 = min(b[2] for _, b, _, _ in track.history)
    x1 = max(b[3] for _, b, _, _ in track.history)

    acc = np.zeros((y1 - y0 + 1, x1 - x0 + 1))
    count = np.zeros_like(acc)
    mask_union = np.zeros(acc.shape, dtype=bool)

    for _frame_idx, bbox, local_mask, local_crop in track.history:
        by0, _by1, bx0, _bx1 = bbox
        oy0, ox0 = by0 - y0, bx0 - x0
        h, w = local_crop.shape
        acc[oy0:oy0 + h, ox0:ox0 + w] += np.clip(local_crop, 0.0, None)
        count[oy0:oy0 + h, ox0:ox0 + w] += 1
        mask_union[oy0:oy0 + h, ox0:ox0 + w] |= local_mask

    avg = np.divide(acc, count, out=np.zeros_like(acc), where=count > 0)
    avg[~mask_union] = 0.0
    peak = avg.max()
    if peak > 0:
        avg = avg / peak

    profile = np.zeros(mov_shape[:2])
    profile[y0:y1 + 1, x0:x1 + 1] = avg
    return profile, (y0, y1, x0, x1)


def _promote_candidate(state, track, frame_index):
    profile, _bbox = _build_promoted_profile(track, (state.mov_y, state.mov_x))
    new_id = state.profiles.shape[2]
    state.profiles = np.concatenate([state.profiles, profile[:, :, np.newaxis]], axis=2)
    state.first_detected_frame[new_id] = frame_index

    y0, y1, x0, x1 = _cell_window_bounds(profile, state.fit.pad_space, state.mov_y, state.mov_x, state.fit.use_com)
    _add_cell_setup(state, new_id, y0, y1, x0, x1)
    _invalidate_overlapping_setups(state, new_id, y0, y1, x0, x1)
    _update_known_cell_exclude_mask(state, new_id, y0, y1, x0, x1)
    return new_id


def _fit_cell(state, frame, cell_id):
    setup = state._cell_setups[cell_id]
    y0, y1, x0, x1 = setup['y0'], setup['y1'], setup['x0'], setup['x1']
    this_frame = frame[y0:y1 + 1, x0:x1 + 1].ravel()
    _tc_lsq, fit_fancy, _fit_x, _lsq_cost, _bob_cost = _solve_one_frame_cell(
        this_frame, setup['rois'], setup['rois_scaled'], setup['lambdas'], setup['norm_factors'],
        setup['k1'], setup['k2'], setup['n_y'], setup['n_x'], state.one_blob, setup['operators'],
        state.use_native, state.fit.native_l_mode, state.fit.native_nthreads,
        state.fit.solver_tol, state.fit.solver_max_iter, state.fit.blob_spacing,
    )
    return float(fit_fancy[setup['cell_index_within']]), (y0, y1, x0, x1)


def realSEUDOfit(frame, state, frame_index=None, zero_level=None):
    """Fit one frame against state's (growing) known-cell set, detect and
    promote any new cells, and return every known cell's activity for this
    frame. Strictly causal: only ever looks at `frame` and state accumulated
    from prior calls -- no lookahead, no whole-movie access.

    frame: (mov_y, mov_x) array -- one movie frame.
    state: a StreamingState, mutated in place.
    frame_index/zero_level: default to state.frame_index/state.zero_level.

    Returns a FrameResult(frame_index, activity, new_cells).
    """
    if frame_index is None:
        frame_index = state.frame_index
    if zero_level is None:
        zero_level = state.zero_level

    frame = np.asarray(frame, dtype=float) - zero_level
    if frame.shape != (state.mov_y, state.mov_x):
        raise ValueError(f'frame shape {frame.shape} does not match state shape {(state.mov_y, state.mov_x)}')

    activity = {}
    reconstruction = np.zeros((state.mov_y, state.mov_x), dtype=float)

    cell_ids = list(state._cell_setups.keys())
    if state._executor is not None and len(cell_ids) > 1:
        # each cell's fit is fully independent (own read-only setup, no
        # shared mutable state -- see _fit_cell/_solve_one_frame_cell), so
        # safe to run concurrently; the native solver releases the GIL for
        # its whole call, so this gets real multi-core speedup even though
        # these are Python threads, not processes
        futures = {state._executor.submit(_fit_cell, state, frame, cell_id): cell_id for cell_id in cell_ids}
        fit_results = {futures[fut]: fut.result() for fut in futures}
    else:
        fit_results = {cell_id: _fit_cell(state, frame, cell_id) for cell_id in cell_ids}

    # accumulate into the shared reconstruction sequentially, in this
    # thread, after every fit has returned -- cells' windows can overlap in
    # pixel space, so concurrent += here would be a real data race
    for cell_id, (value, (y0, y1, x0, x1)) in fit_results.items():
        activity[cell_id] = value
        reconstruction[y0:y1 + 1, x0:x1 + 1] += state.profiles[y0:y1 + 1, x0:x1 + 1, cell_id] * value

    residual = frame - reconstruction
    noise_map = _estimate_noise_map(residual, state.detection.noise_grid_shape)

    # detect on R1 (raw residual after known-cell fitting) and advance every
    # existing track by centroid-matching against it -- unchanged from
    # before two-stage detection; see _update_candidate_tracks.
    smoothed = convolve2d(residual, state.one_blob, mode='same')
    raw_candidates = _run_tile_detection(state, smoothed, noise_map, state.known_cell_exclude_mask)
    _update_candidate_tracks(state, raw_candidates, residual, frame_index)

    # two-stage detection: fit every currently-active track's own evolving
    # profile against R1 and subtract out whatever it explains -- see
    # _subtract_tracked_contributions. Only what's left (R2) gets scanned
    # for brand-new candidates, so an already-tracked candidate's own
    # residual (already accounted for above) can never spawn a duplicate.
    residual2, candidate_exclude_mask = _subtract_tracked_contributions(state, residual)
    combined_exclude_mask = state.known_cell_exclude_mask | candidate_exclude_mask
    smoothed2 = convolve2d(residual2, state.one_blob, mode='same')
    raw_candidates_new = _run_tile_detection(state, smoothed2, noise_map, combined_exclude_mask)
    _start_candidate_tracks(state, raw_candidates_new, residual, frame_index)

    new_cells = []
    ready = [tid for tid, t in state.candidate_tracks.items()
             if t.consecutive_frames >= state.promotion.consecutive_frames_required]
    for track_id in ready:
        track = state.candidate_tracks.pop(track_id)
        new_id = _promote_candidate(state, track, frame_index)
        value, _bbox = _fit_cell(state, frame, new_id)
        activity[new_id] = value
        new_cells.append(new_id)

    state.frame_index = frame_index + 1
    return FrameResult(frame_index=frame_index, activity=activity, new_cells=new_cells)
