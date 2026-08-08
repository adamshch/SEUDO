"""Tests for the streaming/online SEUDO pipeline (seudo/streaming.py):
realSEUDOfit() fits one frame at a time against a growing cell-profile set,
detecting and promoting new cells as they appear.
"""

import numpy as np
import pytest

from seudo._native import NATIVE_AVAILABLE
from seudo.estimate import estimate_time_courses_with_seudo
from seudo.geometry import compute_roi_coms
from seudo.streaming import DetectionParams, FitParams, StreamingState, TilingConfig, realSEUDOfit

requires_native = pytest.mark.skipif(not NATIVE_AVAILABLE, reason='native extension not built')

COMMON_FIT = dict(p=1e-4, sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, pad_space=6,
                   solver_tol=0.01, solver_max_iter=500, use_native=False)


def gaussian_patch(shape, center, sigma, amplitude=1.0):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return amplitude * np.exp(-(((yy - center[0]) ** 2 + (xx - center[1]) ** 2) / (2 * sigma ** 2)))


def test_parity_with_offline_solver_frame_by_frame():
    # one known cell, one spatially-overlapping unmodeled contaminant --
    # same scenario shape as test_estimate_synthetic.py's, exercising the
    # same LSQ-vs-sparse decision branching, not just a trivial single-profile fit.
    rng = np.random.default_rng(0)
    y, x, f = 30, 30, 40
    cell_profile = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    cell_profile[cell_profile < 0.05 * cell_profile.max()] = 0.0
    contam_profile = gaussian_patch((y, x), center=(19, 19), sigma=2.0)
    contam_profile[contam_profile < 0.05 * contam_profile.max()] = 0.0

    cell_activity = np.zeros(f)
    contam_activity = np.zeros(f)
    cell_activity[5:10] = 4.0
    contam_activity[15:20] = 12.0
    cell_activity[25:30] = 4.0
    contam_activity[25:30] = 12.0

    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = (cell_profile * cell_activity[ff] + contam_profile * contam_activity[ff]
                            + 0.02 * rng.normal(size=(y, x)))

    profiles = cell_profile[:, :, np.newaxis]

    offline = estimate_time_courses_with_seudo(movie, profiles, ds_time=1, **COMMON_FIT)
    offline_tc = offline['tc'][:, 0]

    state = StreamingState((y, x), initial_profiles=profiles, fit=FitParams(**COMMON_FIT))
    streaming_tc = np.array([realSEUDOfit(movie[:, :, ff], state).activity[0] for ff in range(f)])

    np.testing.assert_allclose(streaming_tc, offline_tc, atol=1e-8)


def make_two_cell_movie(seed=1, onset=50, y=30, x=60, f=100):
    yy, xx = np.mgrid[0:y, 0:x]
    prof0 = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof0[prof0 < 0.05 * prof0.max()] = 0.0
    prof1 = gaussian_patch((y, x), center=(15, 45), sigma=2.0)
    prof1[prof1 < 0.05 * prof1.max()] = 0.0

    rng = np.random.default_rng(seed)
    act0 = np.full(f, 2.0)
    act1 = np.zeros(f)
    act1[onset:] = 3.0

    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof0 * act0[ff] + prof1 * act1[ff] + 0.01 * rng.normal(size=(y, x))

    return movie, prof0, prof1, act0, act1


def default_streaming_fit():
    return FitParams(sigma2=0.002, lambda_blob=10.0, blob_radius=3.0, pad_space=5, use_native=False)


def test_bounded_latency_new_cell_detection():
    onset = 50
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset)
    y, x, f = movie.shape

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=default_streaming_fit())

    promoted_at = None
    for ff in range(f):
        result = realSEUDOfit(movie[:, :, ff], state)
        if result.new_cells:
            assert promoted_at is None, 'should only promote once'
            promoted_at = ff
            assert result.new_cells == [1]
            assert 1 in result.activity

    assert promoted_at is not None, 'the new cell was never detected'
    assert onset <= promoted_at <= onset + 10, (
        f'promotion latency too high: onset={onset}, promoted at {promoted_at}')

    # activity tracks the true (steady) amplitude after promotion
    late_activity = [realSEUDOfit(movie[:, :, ff], state).activity[1] for ff in range(80, 90)]
    assert np.allclose(late_activity, 3.0, atol=0.5)


def test_contamination_never_promoted():
    y, x, f = 30, 60, 100
    prof0 = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof0[prof0 < 0.05 * prof0.max()] = 0.0
    prof1 = gaussian_patch((y, x), center=(15, 45), sigma=2.0)
    prof1[prof1 < 0.05 * prof1.max()] = 0.0

    rng = np.random.default_rng(2)
    act0 = np.full(f, 2.0)
    act1 = np.zeros(f)
    # flickers on for 2 frames, off for 3, repeating -- never sustains
    # PromotionParams.consecutive_frames_required (default 5) consecutive matches
    for i in range(0, f, 5):
        act1[i:i + 2] = 3.0

    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof0 * act0[ff] + prof1 * act1[ff] + 0.01 * rng.normal(size=(y, x))

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=default_streaming_fit())
    for ff in range(f):
        result = realSEUDOfit(movie[:, :, ff], state)
        assert result.new_cells == []

    assert state.profiles.shape[2] == 1


def test_promoted_profile_recovers_true_shape():
    onset = 50
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset)
    y, x, f = movie.shape

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=default_streaming_fit())
    for ff in range(f):
        realSEUDOfit(movie[:, :, ff], state)

    assert state.profiles.shape[2] == 2
    true_coms, _ = compute_roi_coms(prof1)
    built_coms, _ = compute_roi_coms(state.profiles[:, :, 1])
    np.testing.assert_allclose(built_coms[0], true_coms[0], atol=1.0)


def test_tiling_boundary_no_double_detect_or_miss():
    y, x, f = 40, 40, 60
    prof = gaussian_patch((y, x), center=(20, 20), sigma=2.0)  # exactly on a 20x20 tile-grid corner
    prof[prof < 0.05 * prof.max()] = 0.0

    rng = np.random.default_rng(3)
    act = np.zeros(f)
    act[10:] = 3.0
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof * act[ff] + 0.01 * rng.normal(size=(y, x))

    tiling = TilingConfig(tile_shape=(20, 20), overlap=15)
    state = StreamingState((y, x), tiling=tiling, fit=default_streaming_fit())
    for ff in range(f):
        realSEUDOfit(movie[:, :, ff], state)

    assert state.profiles.shape[2] == 1, 'should promote exactly once, not zero or two'


def test_tiling_vs_no_tiling_agreement():
    onset = 50
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset)
    y, x, f = movie.shape

    fit = default_streaming_fit()

    state_untiled = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=fit)
    untiled_activity0 = []
    promoted_untiled = None
    for ff in range(f):
        r = realSEUDOfit(movie[:, :, ff], state_untiled)
        untiled_activity0.append(r.activity[0])
        if r.new_cells:
            promoted_untiled = ff

    state_tiled = StreamingState(
        (y, x), initial_profiles=prof0[:, :, np.newaxis], fit=fit,
        tiling=TilingConfig(tile_shape=(15, 20), overlap=15))
    tiled_activity0 = []
    promoted_tiled = None
    for ff in range(f):
        r = realSEUDOfit(movie[:, :, ff], state_tiled)
        tiled_activity0.append(r.activity[0])
        if r.new_cells:
            promoted_tiled = ff

    assert promoted_untiled == promoted_tiled
    assert state_untiled.profiles.shape[2] == state_tiled.profiles.shape[2] == 2
    np.testing.assert_allclose(untiled_activity0, tiled_activity0, atol=1e-8)


def test_parallel_tile_detection_matches_sequential():
    # a cell centered exactly on a tile-grid corner exercises multiple
    # tiles' detection genuinely (see test_tiling_boundary_no_double_detect_or_miss)
    y, x, f = 40, 40, 60
    prof = gaussian_patch((y, x), center=(20, 20), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    rng = np.random.default_rng(3)
    act = np.zeros(f)
    act[10:] = 3.0
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof * act[ff] + 0.01 * rng.normal(size=(y, x))

    tiling = TilingConfig(tile_shape=(20, 20), overlap=15)

    def run(n_jobs):
        fit = FitParams(**{**vars(default_streaming_fit()), 'n_jobs': n_jobs})
        state = StreamingState((y, x), tiling=tiling, fit=fit)
        promoted_at = None
        for ff in range(f):
            r = realSEUDOfit(movie[:, :, ff], state)
            if r.new_cells:
                promoted_at = ff
        state.close()
        return promoted_at, state.profiles.shape[2]

    seq = run(1)
    par = run(4)
    assert seq == par


def test_blanking_skips_expensive_detection_for_empty_tiles(monkeypatch):
    # a tile with no signal at all should never reach ndimage.label --
    # verifies the "blanking" pre-filter is actually exercised, not just
    # present in the source with no effect
    import seudo.streaming as streaming_mod

    y, x, f = 40, 40, 15
    movie = np.zeros((y, x, f))  # pure background, no cell anywhere

    calls = []
    orig_label = streaming_mod.ndimage.label

    def spy_label(*args, **kwargs):
        calls.append(1)
        return orig_label(*args, **kwargs)

    monkeypatch.setattr(streaming_mod.ndimage, 'label', spy_label)

    state = StreamingState((y, x), fit=default_streaming_fit())
    for ff in range(f):
        realSEUDOfit(movie[:, :, ff], state)

    assert calls == [], 'ndimage.label should never be reached for an all-background frame'


def test_single_active_tile_never_reaches_the_executor(monkeypatch):
    # the sparse-activity regression this split was built to fix: with only
    # one tile ever non-blank, detection should never touch the thread pool
    # at all -- not "submit it and let it be fast," genuinely never call
    # submit() for detection, since submission/sync overhead on an
    # already-cheap tile is a net loss (see commit history/memory notes)
    y, x, f = 120, 120, 10
    # cell centered well inside one tile's core (>overlap px from every edge
    # of its own tile), so its signal can never leak into a neighboring
    # tile's halo -- see build_tiles: tile_shape=60/overlap=15 gives tile 0
    # core (0,59,0,59), and (30,30) is 30px from the nearest edge
    prof = gaussian_patch((y, x), center=(30, 30), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    rng = np.random.default_rng(6)
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof * 3.0 + 0.01 * rng.normal(size=(y, x))

    fit = FitParams(**{**vars(default_streaming_fit()), 'n_jobs': 4})
    tiling = TilingConfig(tile_shape=(60, 60), overlap=15)  # 2x2 = 4 tiles, only tile 0 ever near the cell
    state = StreamingState((y, x), tiling=tiling, fit=fit)  # no known cells -- isolates detection submissions

    submit_calls = []
    orig_submit = state._executor.submit

    def spy_submit(*args, **kwargs):
        submit_calls.append(args[0])  # the submitted function
        return orig_submit(*args, **kwargs)

    monkeypatch.setattr(state._executor, 'submit', spy_submit)

    for ff in range(f):
        realSEUDOfit(movie[:, :, ff], state)
    state.close()

    detect_submissions = [c for c in submit_calls if c.__name__ == '_detect_in_tile']
    assert detect_submissions == [], 'a single non-blank tile should never be submitted to the executor'


def make_multi_cell_movie(seed=4, n_cells=6, y=40, x=80, f=30):
    rng = np.random.default_rng(seed)
    xs = np.linspace(10, x - 10, n_cells)
    profiles = np.zeros((y, x, n_cells))
    for i, cx in enumerate(xs):
        p = gaussian_patch((y, x), center=(y // 2, cx), sigma=2.0)
        p[p < 0.05 * p.max()] = 0.0
        profiles[:, :, i] = p

    activity = rng.uniform(0.5, 4.0, size=(f, n_cells))
    movie = np.zeros((y, x, f))
    for ff in range(f):
        frame = np.zeros((y, x))
        for i in range(n_cells):
            frame += profiles[:, :, i] * activity[ff, i]
        movie[:, :, ff] = frame + 0.01 * rng.normal(size=(y, x))

    return movie, profiles


def test_n_jobs_gives_bit_identical_results_to_sequential():
    movie, profiles = make_multi_cell_movie()
    y, x, f = movie.shape

    fit_seq = default_streaming_fit()
    fit_par = FitParams(**{**vars(fit_seq), 'n_jobs': 4})

    state_seq = StreamingState((y, x), initial_profiles=profiles, fit=fit_seq)
    state_par = StreamingState((y, x), initial_profiles=profiles, fit=fit_par)
    try:
        for ff in range(f):
            r_seq = realSEUDOfit(movie[:, :, ff], state_seq)
            r_par = realSEUDOfit(movie[:, :, ff], state_par)
            assert r_seq.activity.keys() == r_par.activity.keys()
            for cell_id in r_seq.activity:
                assert r_seq.activity[cell_id] == r_par.activity[cell_id]
    finally:
        state_par.close()


def test_executor_lifecycle():
    y, x = 20, 20
    state_seq = StreamingState((y, x), fit=default_streaming_fit())
    assert state_seq._executor is None  # n_jobs=1 -- no pool created

    fit_par = FitParams(**{**vars(default_streaming_fit()), 'n_jobs': 4})
    with StreamingState((y, x), fit=fit_par) as state_par:
        assert state_par._executor is not None
    assert state_par._executor is None  # closed on context-manager exit


def test_state_contract_sanity():
    onset = 20
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset, f=40)
    y, x, f = movie.shape

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=default_streaming_fit())
    state_id = id(state)

    n_before = state.profiles.shape[2]
    for ff in range(f):
        result = realSEUDOfit(movie[:, :, ff], state)
        assert id(state) == state_id, 'state should mutate in place, not be replaced'

        if result.new_cells:
            assert state.profiles.shape[2] == n_before + len(result.new_cells)
            for cell_id in result.new_cells:
                assert cell_id in result.activity
            n_before = state.profiles.shape[2]
        else:
            assert state.profiles.shape[2] == n_before


def make_contaminated_frames(seed=5, y=30, x=30, n_frames=20):
    """A single known cell plus a withheld, overlapping 'contaminant' --
    unlike make_multi_cell_movie (clean, no interference, so LSQ always
    wins the cost comparison and the blob dictionary is never engaged
    regardless of spacing), this scenario actually exercises the blob
    branch, which is where blob_spacing has any effect at all."""
    yy, xx = np.mgrid[0:y, 0:x]
    prof0 = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof0[prof0 < 0.05 * prof0.max()] = 0.0
    contaminant = gaussian_patch((y, x), center=(19, 19), sigma=2.0)
    contaminant[contaminant < 0.05 * contaminant.max()] = 0.0

    rng = np.random.default_rng(seed)
    frames = [prof0 * 4.0 + contaminant * 12.0 + 0.02 * rng.normal(size=(y, x)) for _ in range(n_frames)]
    return frames, prof0[:, :, np.newaxis], (y, x)


@requires_native
def test_blob_spacing_changes_native_result():
    frames, profiles, (y, x) = make_contaminated_frames()

    def run(blob_spacing):
        fit = FitParams(sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, pad_space=6,
                         use_native=True, blob_spacing=blob_spacing, solver_max_iter=300)
        state = StreamingState((y, x), initial_profiles=profiles, fit=fit)
        return [realSEUDOfit(frame, state).activity[0] for frame in frames]

    fine = run(1)
    coarse = run(4)
    assert not np.allclose(fine, coarse)


def test_blob_spacing_is_a_noop_for_python_fallback():
    frames, profiles, (y, x) = make_contaminated_frames()

    def run(blob_spacing):
        fit = FitParams(sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, pad_space=6,
                         use_native=False, blob_spacing=blob_spacing, solver_max_iter=300)
        state = StreamingState((y, x), initial_profiles=profiles, fit=fit)
        return [realSEUDOfit(frame, state).activity[0] for frame in frames]

    fine = run(1)
    coarse = run(4)
    np.testing.assert_array_equal(fine, coarse)


def test_noise_grid_shape_default_matches_old_single_scalar_behavior():
    # noise_grid_shape=(1,1) (the default) must reproduce the pre-existing
    # single-global-scalar behavior exactly -- this is the backward
    # compatibility guarantee for every earlier test in this file
    onset = 20
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset, f=40)
    y, x, f = movie.shape

    fit = default_streaming_fit()
    default_detection = DetectionParams()
    assert default_detection.noise_grid_shape == (1, 1)

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=fit, detection=default_detection)
    activity = [realSEUDOfit(movie[:, :, ff], state).activity[0] for ff in range(f)]

    state_ref = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=fit)
    activity_ref = [realSEUDOfit(movie[:, :, ff], state_ref).activity[0] for ff in range(f)]

    np.testing.assert_array_equal(activity, activity_ref)


def test_spatially_varying_noise_detects_modest_cell_faster_in_quiet_region():
    # a modest-amplitude cell in a quiet region, next to an unrelated noisy
    # region -- a single global noise scalar gets dragged up by the noisy
    # region, delaying detection of the modest cell; per-section estimation
    # should detect it much faster since the quiet region's own noise floor
    # is correctly low from the start
    y, x, f = 60, 60, 30
    yy, xx = np.mgrid[0:y, 0:x]
    prof = gaussian_patch((y, x), center=(30, 15), sigma=2.0)  # in the quiet left half
    prof[prof < 0.05 * prof.max()] = 0.0

    rng = np.random.default_rng(2)
    act = np.zeros(f)
    act[10:] = 1.0  # modest amplitude

    frames = []
    for ff in range(f):
        frame = prof * act[ff]
        noise = np.zeros((y, x))
        noise[:, :30] = 0.01 * rng.normal(size=(y, 30))       # quiet left half
        noise[:, 30:] = 0.15 * rng.normal(size=(y, x - 30))   # noisy right half
        frames.append(frame + noise)

    def run(noise_grid_shape):
        fit = default_streaming_fit()
        detection = DetectionParams(noise_grid_shape=noise_grid_shape)
        state = StreamingState((y, x), fit=fit, detection=detection)
        promoted_at = None
        for ff in range(f):
            r = realSEUDOfit(frames[ff], state)
            if r.new_cells:
                promoted_at = ff
        return promoted_at, state.profiles.shape[2]

    global_promoted, global_n = run((1, 1))
    split_promoted, split_n = run((1, 2))  # left/right halves get separate noise estimates

    assert global_n == split_n == 1, 'both should eventually find the cell'
    assert split_promoted < global_promoted, (
        'spatially-varying noise estimation should detect the modest cell faster '
        'than a global scalar dragged up by the unrelated noisy region')


def test_subtract_tracked_contributions_removes_tracked_signal():
    # a single active (not-yet-promoted) candidate track whose accumulated
    # history already captures a clean cell shape -- stage 1 of two-stage
    # detection should recognize and remove most of its contribution from a
    # fresh residual, not leave it untouched (a no-op would defeat the whole
    # point: telling an already-tracked candidate's own signal apart from
    # something genuinely new -- see _subtract_tracked_contributions).
    from seudo.streaming import CandidateTrack, _subtract_tracked_contributions

    y, x = 30, 30
    prof = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    state = StreamingState((y, x), fit=default_streaming_fit())
    bbox = (10, 20, 10, 20)
    crop = prof[bbox[0]:bbox[1] + 1, bbox[2]:bbox[3] + 1]
    mask = crop > 0
    history_crop = crop * 3.0
    track = CandidateTrack(
        centroid=(15.0, 15.0), consecutive_frames=3, gap=0,
        history=[(0, bbox, mask, history_crop)] * 3,
    )
    state.candidate_tracks[0] = track

    rng = np.random.default_rng(0)
    residual = prof * 3.0 + 0.01 * rng.normal(size=(y, x))
    residual2, exclude_mask, track_footprints = _subtract_tracked_contributions(state, residual)

    footprint = prof > 0
    assert np.abs(residual2[footprint]).sum() < 0.5 * np.abs(residual[footprint]).sum(), (
        'stage 1 should explain away most of an already-tracked candidate\'s own signal')
    assert exclude_mask[15, 15], 'the tracked footprint should be excluded from new-candidate scanning'


def test_subtract_tracked_contributions_fits_overlapping_tracks_jointly():
    # two active tracks whose profiles genuinely overlap (spatially close
    # Gaussians) -- fitting each in total isolation against a residual that
    # actually contains BOTH signals would misattribute the overlap region
    # (each track's lone-column LSQ has no way to know the other profile's
    # contribution exists), leaving real residual behind. A joint fit
    # (both tracks sharing one window, like known cells already do) should
    # correctly decompose the two and drive the residual down everywhere.
    from seudo.streaming import CandidateTrack, _subtract_tracked_contributions

    y, x = 30, 30
    prof_a = gaussian_patch((y, x), center=(15, 13), sigma=2.0)
    prof_a[prof_a < 0.05 * prof_a.max()] = 0.0
    prof_b = gaussian_patch((y, x), center=(15, 17), sigma=2.0)  # 4px away -- real overlap
    prof_b[prof_b < 0.05 * prof_b.max()] = 0.0

    state = StreamingState((y, x), fit=default_streaming_fit())

    def make_track(prof, bbox):
        crop = prof[bbox[0]:bbox[1] + 1, bbox[2]:bbox[3] + 1]
        mask = crop > 0
        return CandidateTrack(centroid=(15.0, 15.0), consecutive_frames=3, gap=0,
                               history=[(0, bbox, mask, crop * 3.0)] * 3)

    bbox_a = (10, 20, 8, 18)
    bbox_b = (10, 20, 12, 22)
    state.candidate_tracks[0] = make_track(prof_a, bbox_a)
    state.candidate_tracks[1] = make_track(prof_b, bbox_b)

    true_weight_a, true_weight_b = 3.0, 5.0
    residual = prof_a * true_weight_a + prof_b * true_weight_b  # clean, no noise

    residual2, _exclude_mask, _track_footprints = _subtract_tracked_contributions(state, residual)

    combined_footprint = (prof_a > 0) | (prof_b > 0)
    assert np.abs(residual2[combined_footprint]).max() < 0.1 * max(true_weight_a, true_weight_b), (
        'jointly fitting overlapping tracks should decompose both signals correctly, '
        'not leave real residual from misattributing the overlap region')


def test_tracked_candidate_never_spawns_a_duplicate_track():
    # while a real cell is being tracked (pre-promotion), the two-stage
    # subtraction should keep explaining its own signal away every frame --
    # if it didn't, leftover residual at the same location could spuriously
    # trigger a second, duplicate track for the same physical source
    onset = 30
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset, f=onset + 15)
    y, x, f = movie.shape

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=default_streaming_fit())
    for ff in range(f):
        realSEUDOfit(movie[:, :, ff], state)

    assert state._next_track_id == 1, (
        f'expected exactly one track ever created for the single new cell, got {state._next_track_id}')


def test_ds_time_default_is_a_noop():
    # ds_time=1 (the default) must apply no averaging at all -- sigma2_ds
    # equals the raw sigma2 exactly, and a single call's avg_frame equals
    # the input frame exactly (np.mean of one array introduces no
    # floating-point change), so behavior is bit-identical to every test
    # above that never mentions ds_time at all.
    y, x = 20, 20
    fit = FitParams(sigma2=0.02)
    state = StreamingState((y, x), fit=fit)
    assert state.sigma2_ds == fit.sigma2

    frame = np.random.default_rng(0).normal(size=(y, x))
    state._frame_buffer.append(frame)
    avg_frame = np.mean(state._frame_buffer, axis=0)
    np.testing.assert_array_equal(avg_frame, frame)


def test_ds_time_smooths_a_stepped_activity_trace():
    # a single known cell with a noise-free step in true activity -- with
    # ds_time=3, the FIRST frame after the step should read back as a
    # partial (averaged-with-the-past) value, not the new true value
    # immediately, and should reach the new true value only once the
    # trailing window has fully slid past the step -- direct evidence the
    # averaging math is doing the right causal thing, not just "some
    # smoothing happens somewhere."
    y, x, f = 30, 30, 20
    prof = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    act = np.full(f, 2.0)
    step_at = 10
    act[step_at:] = 8.0

    movie = np.stack([prof * act[ff] for ff in range(f)], axis=2)

    fit = FitParams(**{**vars(default_streaming_fit()), 'ds_time': 3})
    state = StreamingState((y, x), initial_profiles=prof[:, :, np.newaxis], fit=fit)
    activity = [realSEUDOfit(movie[:, :, ff], state).activity[0] for ff in range(f)]

    # before the step: steady at the true value
    assert np.allclose(activity[step_at - 1], 2.0, atol=1e-6)
    # right at the step: averages 1 new frame at 8.0 with 2 trailing frames at 2.0
    assert np.allclose(activity[step_at], (8.0 + 2.0 + 2.0) / 3.0, atol=1e-6)
    # one frame later: 2 frames at 8.0, 1 trailing frame at 2.0
    assert np.allclose(activity[step_at + 1], (8.0 + 8.0 + 2.0) / 3.0, atol=1e-6)
    # once the window has fully slid past the step: back to steady, at the new value
    assert np.allclose(activity[step_at + 2], 8.0, atol=1e-6)
    assert np.allclose(activity[-1], 8.0, atol=1e-6)


def test_ds_time_improves_detection_in_high_noise():
    # a modest-amplitude cell under noise high enough to make single-frame
    # detection unreliable -- ds_time=3's trailing average should reduce
    # the effective per-pixel noise and detect it, while ds_time=1 (no
    # averaging) fails to promote at all within the same frame budget
    y, x, f = 30, 30, 100
    prof = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    rng = np.random.default_rng(7)
    act = np.zeros(f)
    act[10:] = 1.0

    frames = [prof * act[ff] + 0.25 * rng.normal(size=(y, x)) for ff in range(f)]

    def run(ds_time):
        fit = FitParams(**{**vars(default_streaming_fit()), 'ds_time': ds_time})
        state = StreamingState((y, x), fit=fit)
        promoted_at = None
        for ff in range(f):
            r = realSEUDOfit(frames[ff], state)
            if r.new_cells:
                promoted_at = ff
        return promoted_at, state.profiles.shape[2]

    no_avg_promoted, no_avg_n = run(1)
    smoothed_promoted, smoothed_n = run(3)

    assert smoothed_n == 1, 'ds_time=3 should reliably detect the modest cell under high noise'
    assert no_avg_n == 0 or (no_avg_promoted is not None and smoothed_promoted < no_avg_promoted), (
        'ds_time=3 should detect the cell either faster than ds_time=1, or where ds_time=1 fails entirely')


def test_blobify_radius_defaults_to_mirroring_fit_blob_radius():
    # blobify_radius=None (the default) should build detect_blob identical
    # to one_blob -- the exact coupled behavior this module always had
    # before blobify_radius existed as an independent knob
    y, x = 20, 20
    fit = FitParams(blob_radius=2.5)
    state = StreamingState((y, x), fit=fit)
    np.testing.assert_array_equal(state.detect_blob, state.one_blob)


def test_blobify_radius_override_produces_a_different_kernel():
    # an explicit blobify_radius should genuinely decouple detection
    # smoothing from the fit's own blob dictionary, not silently be ignored
    y, x = 20, 20
    fit = FitParams(blob_radius=3.0)
    detection = DetectionParams(blobify_radius=1.2)
    state = StreamingState((y, x), fit=fit, detection=detection)
    assert state.detect_blob.shape != state.one_blob.shape


def test_should_merge_temp_profiles_eq8():
    # direct check of the paper's Eq. 8 merge condition (realSEUDO sec 3.2),
    # independent of the rest of the tracking machinery
    from seudo.streaming import _should_merge_temp_profiles

    # profile B is a 3x3 block fully contained within a 10x10 profile A --
    # B's every pixel is shared, so U_b=0 <= B_b*0.5 trivially: merge
    mask_a = np.ones((10, 10), dtype=bool)
    bbox_a = (0, 9, 0, 9)
    mask_b = np.ones((3, 3), dtype=bool)
    bbox_b = (4, 6, 4, 6)
    assert _should_merge_temp_profiles(mask_a, bbox_a, mask_b, bbox_b) is True

    # two same-size blocks overlapping by only one corner pixel out of 9 --
    # low common fraction, high unique fraction relative to their own small
    # perimeters: should NOT merge, this is two separate small profiles
    mask_c = np.ones((3, 3), dtype=bool)
    bbox_c = (0, 2, 0, 2)
    mask_d = np.ones((3, 3), dtype=bool)
    bbox_d = (2, 4, 2, 4)
    assert _should_merge_temp_profiles(mask_c, bbox_c, mask_d, bbox_d) is False

    # disjoint bounding boxes: never merge
    mask_e = np.ones((3, 3), dtype=bool)
    bbox_e = (0, 2, 0, 2)
    mask_f = np.ones((3, 3), dtype=bool)
    bbox_f = (20, 22, 20, 22)
    assert _should_merge_temp_profiles(mask_e, bbox_e, mask_f, bbox_f) is False


def test_start_candidate_tracks_absorbs_a_candidate_overlapping_an_existing_track():
    # a raw candidate that's really just leftover, imperfectly-subtracted
    # residual from a track already being followed (Eq. 8 merge condition
    # satisfied) should be silently absorbed, not spawn a duplicate track --
    # this is the mechanism meant to compensate for a tight
    # exclude_radius_known_cells not fully covering a track's footprint
    from seudo.streaming import RawCandidate, _start_candidate_tracks

    y, x = 30, 30
    state = StreamingState((y, x), fit=default_streaming_fit())
    existing_mask = np.ones((10, 10), dtype=bool)
    existing_bbox = (5, 14, 5, 14)
    track_footprints = [(0, existing_mask, existing_bbox)]

    # a small candidate almost entirely inside the existing track's bbox
    overlapping_cand = RawCandidate(
        centroid=(9.0, 9.0), bbox=(8, 10, 8, 10), mask=np.ones((3, 3), dtype=bool))

    n_before = state._next_track_id
    _start_candidate_tracks(state, [overlapping_cand], np.zeros((y, x)), 0, track_footprints)
    assert state._next_track_id == n_before, 'an Eq.8-merging candidate should not start a new track'
    assert len(state.candidate_tracks) == 0


def test_start_candidate_tracks_still_creates_a_track_for_a_disjoint_candidate():
    # sanity check the Eq. 8 filter doesn't over-suppress -- a candidate far
    # from every existing track should still start a genuinely new one
    from seudo.streaming import RawCandidate, _start_candidate_tracks

    y, x = 30, 30
    state = StreamingState((y, x), fit=default_streaming_fit())
    existing_mask = np.ones((5, 5), dtype=bool)
    existing_bbox = (2, 6, 2, 6)
    track_footprints = [(0, existing_mask, existing_bbox)]

    far_cand = RawCandidate(centroid=(25.0, 25.0), bbox=(23, 27, 23, 27), mask=np.ones((5, 5), dtype=bool))

    n_before = state._next_track_id
    _start_candidate_tracks(state, [far_cand], np.zeros((y, x)), 0, track_footprints)
    assert state._next_track_id == n_before + 1
    assert len(state.candidate_tracks) == 1


def test_spatial_denoise_radius_default_is_a_noop():
    y, x = 20, 20
    state = StreamingState((y, x), fit=FitParams())
    assert state.denoise_kernel is None


def test_spatial_denoise_radius_builds_a_sum_normalized_kernel():
    # a smoothing filter must preserve overall intensity (sum to ~1), unlike
    # make_seudo_blob's L2-normalized SEUDO basis-function convention
    y, x = 20, 20
    state = StreamingState((y, x), fit=FitParams(spatial_denoise_radius=1.5))
    assert state.denoise_kernel is not None
    np.testing.assert_allclose(state.denoise_kernel.sum(), 1.0, atol=1e-10)


def test_spatial_denoise_radius_changes_known_cell_activity_under_pixel_noise():
    # a single, sharp, few-pixel-wide noise spike landing inside a known
    # cell's own fit window -- spatial denoising should measurably dilute
    # it before the fit ever sees it, changing the recovered activity
    # relative to fitting the raw (undenoised) frame
    y, x = 30, 30
    prof = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    frame = prof * 3.0
    frame[12, 12] += 20.0  # sharp spike, inside the cell's own padded fit window

    def run(spatial_denoise_radius):
        fit = FitParams(**{**vars(default_streaming_fit()), 'spatial_denoise_radius': spatial_denoise_radius})
        state = StreamingState((y, x), initial_profiles=prof[:, :, np.newaxis], fit=fit)
        return realSEUDOfit(frame, state).activity[0]

    undenoised = run(None)
    denoised = run(1.5)
    assert not np.isclose(undenoised, denoised), (
        'spatial denoising should measurably change the fit when a sharp spike is present')


def test_max_roi_extent_rejects_a_large_diffuse_region():
    # a broad, roughly uniform elevated region (like a real illumination
    # gradient, not a compact cell) spanning much of the frame -- without a
    # cap this crosses threshold as one giant contiguous blob and gets
    # tracked/promoted; max_roi_extent should reject it, while a normal-
    # sized real cell elsewhere is unaffected. Note: the detected extent of
    # even a compact real cell is inherently wider than its own true size,
    # since detection thresholds the blob-convolved + dilated residual, not
    # the raw cell -- the cap (35) is chosen comfortably above that
    # inflated-but-real extent (~27px here at blob_radius=3), well below
    # the gradient region's much larger one.
    y, x, f = 120, 120, 10
    gradient = np.zeros((y, x))
    gradient[:60, :60] = 3.0

    cell = gaussian_patch((y, x), center=(95, 95), sigma=2.0)
    cell[cell < 0.05 * cell.max()] = 0.0

    frames = [gradient + cell * 3.0 for _ in range(f)]
    common = dict(sigma2=0.002, lambda_blob=10.0, blob_radius=3.0, pad_space=5, use_native=False)

    def run(max_roi_extent):
        fit = FitParams(**common)
        detection = DetectionParams(max_roi_extent=max_roi_extent)
        state = StreamingState((y, x), fit=fit, detection=detection)
        max_tracks_ever = 0
        for frame in frames:
            realSEUDOfit(frame, state)
            max_tracks_ever = max(max_tracks_ever, state._next_track_id)
        return state.profiles.shape[2], max_tracks_ever

    n_cells_uncapped, n_tracks_uncapped = run(None)
    n_cells_capped, n_tracks_capped = run(35)

    assert n_tracks_uncapped >= 2, 'sanity check: both the gradient region and the real cell get tracked uncapped'
    assert n_cells_capped == 1, 'only the real cell should survive with the cap on'
    assert n_cells_capped < n_cells_uncapped


def test_rejected_regions_are_quarantined_not_rescanned_every_frame():
    # a noisy (not perfectly static) oversized region -- without permanent
    # quarantine, rejecting it fresh every frame lets different noise-driven
    # sub-regions randomly sneak under the cap on different frames,
    # fragmenting into multiple spurious tracks over time (confirmed this
    # is a REAL effect on real data: capping alone made results measurably
    # worse). With quarantine, the region should be excluded after its
    # first rejection and never contribute a spurious track again, while a
    # real cell elsewhere keeps working normally.
    from seudo.streaming import _quarantine_rejected_regions

    y, x, f = 120, 120, 40
    rng = np.random.default_rng(0)
    cell = gaussian_patch((y, x), center=(95, 95), sigma=2.0)
    cell[cell < 0.05 * cell.max()] = 0.0

    frames = []
    for _ in range(f):
        gradient = 3.0 + rng.normal(0, 0.3, size=(60, 60))
        frame = np.zeros((y, x))
        frame[:60, :60] = gradient
        frame += cell * 3.0
        frames.append(frame)

    common = dict(sigma2=0.002, lambda_blob=10.0, blob_radius=3.0, pad_space=5, use_native=False)
    fit = FitParams(**common)
    detection = DetectionParams(max_roi_extent=35)
    state = StreamingState((y, x), fit=fit, detection=detection)
    for frame in frames:
        realSEUDOfit(frame, state)

    assert state.profiles.shape[2] == 1, 'only the real cell should ever be promoted'
    assert state._next_track_id == 1, (
        'the noisy oversized region should never fragment into spurious tracks once quarantined')
    assert state.rejected_region_mask[:60, :60].sum() > 0, (
        'the quarantine mask should cover (at least part of) the rejected region')


def test_quarantine_rejected_regions_expands_bbox_by_exclude_radius():
    from seudo.streaming import _quarantine_rejected_regions

    y, x = 30, 30
    state = StreamingState((y, x), fit=default_streaming_fit(),
                            detection=DetectionParams(exclude_radius_known_cells=3))
    _quarantine_rejected_regions(state, [(10, 15, 10, 15)])

    assert state.rejected_region_mask[10:16, 10:16].all(), 'the rejected bbox itself must be quarantined'
    assert state.rejected_region_mask[7, 10], 'quarantine should expand by exclude_radius_known_cells'
    assert not state.rejected_region_mask[6, 10], 'expansion should not go further than exclude_radius_known_cells'


def test_eq9_containment_ratios():
    # direct check of the paper's Eq. 9 formula (realSEUDO sec 3.2),
    # independent of the rest of the promotion machinery
    from seudo.streaming import _eq9_containment_ratios

    y, x = 30, 30
    b = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    b[b < 0.05 * b.max()] = 0.0

    # same shape, different amplitude -- A fits entirely inside B and vice versa
    a_same_shape = b * 0.6
    rho_ab, rho_ba = _eq9_containment_ratios(a_same_shape, b)
    np.testing.assert_allclose([rho_ab, rho_ba], [1.0, 1.0], atol=1e-8)

    # A is a strict spatial subset of B's own footprint -- A fits fully
    # inside B (rho_ab=1) but B does NOT fit inside A (rho_ba << 1):
    # exactly the asymmetric signature the paper uses to flag a split
    # candidate, distinct from a clean merge
    ys, xs = np.nonzero(b > 0)
    cy, cx = int(ys.mean()), int(xs.mean())
    mask = np.zeros_like(b, dtype=bool)
    mask[:cy, :cx] = True
    a_partial = np.where(mask, b, 0.0)
    rho_ab, rho_ba = _eq9_containment_ratios(a_partial, b)
    assert rho_ab > 0.99
    assert rho_ba < 0.2

    # no spatial overlap at all
    c = gaussian_patch((y, x), center=(2, 2), sigma=1.0)
    c[c < 0.05 * c.max()] = 0.0
    assert _eq9_containment_ratios(b, c) == (0.0, 0.0)


def test_find_stable_merge_target():
    from seudo.streaming import _find_stable_merge_target

    y, x = 30, 30
    prof = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    state = StreamingState((y, x), initial_profiles=prof[:, :, np.newaxis], fit=default_streaming_fit())

    # a near-duplicate of the existing cell (different amplitude, same shape)
    duplicate = prof * 0.7
    assert _find_stable_merge_target(state, duplicate) == 0

    # two nearby-but-genuinely-distinct cells (same scenario as the joint-
    # fit test) should NOT be flagged as the same cell
    distinct = gaussian_patch((y, x), center=(15, 19), sigma=2.0)
    distinct[distinct < 0.05 * distinct.max()] = 0.0
    assert _find_stable_merge_target(state, distinct) is None


def test_promote_candidate_merges_duplicate_instead_of_creating_new_cell():
    # a candidate track whose built profile turns out to be (nearly) the
    # same shape as an ALREADY-known cell -- Eq. 9 should merge it into the
    # existing cell rather than adding a redundant duplicate profile
    from seudo.streaming import CandidateTrack, _promote_candidate

    y, x = 30, 30
    prof = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    state = StreamingState((y, x), initial_profiles=prof[:, :, np.newaxis], fit=default_streaming_fit())
    n_before = state.profiles.shape[2]

    bbox = (10, 20, 10, 20)
    crop = prof[bbox[0]:bbox[1] + 1, bbox[2]:bbox[3] + 1]
    mask = crop > 0
    track = CandidateTrack(centroid=(15.0, 15.0), consecutive_frames=3, gap=0,
                            history=[(0, bbox, mask, crop * 0.8)] * 3)

    cell_id, is_new = _promote_candidate(state, track, frame_index=0)
    assert is_new is False
    assert cell_id == 0
    assert state.profiles.shape[2] == n_before, 'no duplicate profile should have been added'
